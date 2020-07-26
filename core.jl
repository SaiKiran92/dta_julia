using DataFrames: DataFrame, DataFrame!
using CSV: File
using DataStructures
using Parameters
using Convex, SCS
using Plots

using EllipsisNotation

#=

Notation:

Fᵢ      - inflows
Fₒ      - outflows
Nᵢ      - cumulative inflows
Nₒ      - cumulative outflows
Cᵢ      - expected costs of vehicles coming in
Cₒ      - expected costs of vehicles going out

DTC       - departure time choices
SR      - splitting rates

=#

ROUND_DIGITS = 10

include("types.jl")
include("crange.jl")
include("utils.jl")
include("models.jl")

linkdata = DataFrame!(File("nguyen.csv"; delim=", "))

ugap, bgap = 3, 4; # used in computing sending and receiving functions

linkdata.l *= ugap;
net = Network(linkdata);

nclasses, nsinks, nlinks = 1, 2, numlinks(net);
u, w = 1/ugap, 1/bgap;
Q, N, δ = 24.0, 80, w/u;
demand_level, T, Tm = 1800, 132*ugap, 80*ugap;

M = 10 * T
α, β, γ, trgt = 1., 0.5, 2., 60*ugap

srcs, snks, mrgs, divs = sources(net), sinks(net), merges(net), diverges(net);
trips = Dict((o,d,c) => demand_level for o in srcs for d in snks for c in 1:nclasses);

Fᵢ, Fₒ, Nᵢ, Nₒ, Cᵢ, Cₒ, states, revtracs = nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing

include("initchoices.jl")
include("simulation.jl")
include("costcomputation.jl")
include("equilibrium.jl")

DTC, SR = initchoices();

for i in 1:30
    global DTC, SR, Cᵢ, Cₒ
    #@show i
    Fᵢ, Fₒ, states = simulate()#(DTC, SR)
    rstates, Cᵢ, Cₒ = computecosts()#(states)
    relgap(Cᵢ, DTC)

    @assert approxpositive(Fᵢ) # all((Fᵢ .>= 0.) .| (Fᵢ .≈ 0.))
    @assert approxpositive(Fₒ) # all(Fₒ .>= 0.)
    @assert approxpositive(Cᵢ) # all(Cᵢ .>= 0.)
    @assert approxpositive(Cₒ) # all(Cₒ .>= 0.)

    updatechoices!() #(Cᵢ, DTC, SR)
end

DTC, SR = initchoices()

Fᵢ, Fₒ, states = simulate() #(DTC, SR)
rstates, Cᵢ, Cₒ = computecosts() #(states)
relgap(Cᵢ, DTC)
plot(DTC[1,19,1])

updatechoices!()

relgap(Cᵢ, DTC)
