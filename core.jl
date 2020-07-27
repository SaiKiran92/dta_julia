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

Fᵢ      - inflows
Fₒ      - outflows
TFᵢ     - total (or) aggregated inflows
TFₒ     - total (or) aggregated outflows
CFᵢ     - cumulative inflows
CFₒ     - cumulative outflows
ECᵢ      - expected costs of vehicles coming in
ECₒ      - expected costs of vehicles going out

DTC       - departure time choices
SR      - splitting rates

=#

ROUND_DIGITS = 10

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
trips = Dict((o,d,c) => demand_level for o in srcs for d in snks for c in 1:nclasses);

Fᵢ, Fₒ, CFᵢ, CFₒ, states, ECᵢ, ECₒ = nothing, nothing, nothing, nothing, nothing, nothing, nothing;

include("initchoices.jl");
include("simulation.jl");
include("costcomputation.jl");
include("equilibrium.jl");

DTC, SR = initchoices();
simulate();
computecosts();
relgap(ECᵢ, DTC)

updatechoices!()

relgap(ECᵢ, DTC)



for i in 1:30
    global DTC, SR, ECᵢ, ECₒ
    #@show i
    Fᵢ, Fₒ, states = simulate()#(DTC, SR)
    ECᵢ, ECₒ = computecosts()#(states)
    relgap(ECᵢ, DTC)

    @assert approxpositive(Fᵢ) # all((Fᵢ .>= 0.) .| (Fᵢ .≈ 0.))
    @assert approxpositive(Fₒ) # all(Fₒ .>= 0.)
    @assert approxpositive(ECᵢ) # all(ECᵢ .>= 0.)
    @assert approxpositive(ECₒ) # all(ECₒ .>= 0.)

    updatechoices!() #(ECᵢ, DTC, SR)
end
