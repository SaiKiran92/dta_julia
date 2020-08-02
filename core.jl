using DataFrames: DataFrame, DataFrame!
using CSV: File
using DataStructures
using Parameters
using Convex, SCS
using Plots
using MathOptInterface
const MOI = MathOptInterface

using EllipsisNotation

#=

Notation:

pᵢ      - inflow proportions
CFᵢ     - cumulative inflows
CFₒ     - cumulative outflows
ECᵢ      - expected costs of vehicles coming in
ECₒ      - expected costs of vehicles going out

DTC       - departure time choices
SR      - splitting rates

=#

ROUND_DIGITS = 13

include("types.jl")
include("crange.jl")
include("utils.jl")
include("models.jl")

ugap, bgap = 3, 4 # used in computing sending and receiving functions
u, w = 1/ugap, 1/bgap
Q, N, δ = 24.0, 80, w/u
demand_level, T, Tm = 1800, 132*ugap, 80*ugap
M = 10 * T
α, β, γ, trgt = 1., 0.5, 2., 60*ugap

linkdata = DataFrame!(File("nguyen.csv"; delim=", "))
linkdata.l *= ugap;
net = Network(linkdata);

nclasses, nsinks, nlinks = 1, 2, numlinks(net);
srcs, snks, mrgs, divs = sources(net), sinks(net), merges(net), diverges(net);
trips = ones(length(srcs), length(snks), nclasses) .* demand_level;

include("initchoices.jl");
include("simulation.jl");
include("costcomputation.jl");
include("equilibrium.jl");

DTC, SR = initchoices();

i = 5; DTC, SR = DTC_list[i], SR_list[i];
CFᵢ, CFₒ, states = simulate();
ECᵢ, ECₒ, rtracs = computecosts();
relgap(ECᵢ, DTC)

DTC, SR = updatechoices!()

relgap(ECᵢ, DTC)


DTC_list = [DTC]
SR_list = [SR]

for i in 1:30
    global DTC, SR, CFᵢ, CFₒ, ECᵢ, ECₒ
    #@show i
    CFᵢ, CFₒ, states = simulate()#(DTC, SR)
    ECᵢ, ECₒ, rtracs = computecosts()#(states)
    @show relgap(ECᵢ, DTC)

    #@assert approxpositive(Fᵢ) # all((Fᵢ .>= 0.) .| (Fᵢ .≈ 0.))
    #@assert approxpositive(Fₒ) # all(Fₒ .>= 0.)
    @assert all(ECᵢ .>= 0.)
    @assert all(ECₒ .>= 0.)

    newDTC, newSR = updatechoices!() #(ECᵢ, DTC, SR)
    push!(DTC_list, DTC)
    push!(SR_list, SR)
    DTC, SR = newDTC, newSR
end
