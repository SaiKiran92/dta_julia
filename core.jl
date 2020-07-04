using DataFrames, CSV
using DataStructures
using Parameters

using EllipsisNotation

ROUND_DIGITS = 10

include("types.jl")
include("utils.jl")
include("crange.jl")
include("models.jl")

linkdata = CSV.read("nguyen.csv"; delim=", ")

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

# choice initialization
dtchoices = Dict(k => zeros(T) for k in keys(trips));
for k in keys(trips)
    #dtchoices[k][1:Tm] .= 1/Tm
    dtchoices[k][10] = 1.
end

pathvecs = Dict(snk => dijkstra(net, snk) for snk in snks);
srates = Dict(div => Dict(i => zeros(T, nsinks, nclasses) for i in 1:2) for div in divs);
for div in divs
    for (snkno,snk) in enumerate(snks)
        for cls in 1:nclasses
            if (outneighbors(net, div)[1] == pathvecs[snk][div])
                srates[div][1][:,snkno,cls] .= 1.
            elseif (outneighbors(net, div)[2] == pathvecs[snk][div])
                srates[div][2][:,snkno,cls] .= 1.
            else
                error()
            end
        end
    end
end

# initialization
inflows = zeros(nlinks, T, nsinks, nclasses);
outflows = zeros(nlinks, T, nsinks, nclasses);
for src in sources(net)
    for (snkno,snk) in enumerate(sinks(net))
        for c in 1:nclasses
            inflows[outlinkids(net,src)[1],:,snkno,c] .= dtchoices[src,snk,c] * trips[src,snk,c]
        end
    end
end

srclinks = [outlinkids(net,src)[1] for src in srcs];
states = Matrix{Union{Nothing, LinkState}}(fill(nothing, nlinks, T+1));
for i in 1:nlinks
    l = length(link(net, i))
    states[i,1] = LinkState(-l, 0., (i in srclinks) ? Inf : N*l)
end

nxttracs = ones(nlinks);
for t in 1:T
    nxttracs .*= 0.
    for snk in snks
        i = inlinkids(net,snk)[1]
        lnk = link(net, i)
        l = length(lnk)

        outflows[i,t,:,:] .= (t > l) ? inflows[i,t-l,:,:] : 0.

        nxttracs[i] = 1
    end

    for mrg in mrgs
        # sending and receiving capacities
        ili = inlinkids(net,mrg)
        oli = outlinkids(net,mrg)[1]
        s = [min(Q, svalue(states[li,t])) for li in ili]
        r = min(Q, rvalue(states[oli,t]))

        # calculating flows
        f = [inflows[li, floor(tracker(states[li,t])):(t - length(link(net, li))),..] for li in ili]
        tflows, tracs = mrgflows(f, r, s, decimal.(tracker.(states[ili,t])))

        tflows = round.(tflows, digits=ROUND_DIGITS)

        outflows[ili,t,..] .= tflows
        inflows[oli,t,..] .= sum(tflows, dims=1)[1,..]

        nxttracs[ili] .= tracs
    end

    for div in divs
        # sending and receiving capacities
        ili = inlinkids(net,div)[1]
        oli = outlinkids(net,div)
        s = min(Q, svalue(states[ili,t]))
        r = [min(Q, rvalue(states[li,t])) for li in oli]

        # calculating flows
        sr = [srates[div][ti][t,..] for ti in 1:2]
        tflows, trac = divflows(inflows[ili,floor(tracker(states[ili,t])):(t - length(link(net, ili))),..], sr,r,s,decimal(tracker(states[ili,t])))

        tflows = round.(tflows, digits=ROUND_DIGITS)
        #trac = round(trac, digits=ROUND_DIGITS)

        outflows[ili,t,..] .= sum(tflows,dims=1)[1,..]
        inflows[oli,t,..] .= tflows

        nxttracs[ili] = trac
    end

    # update link states
    for (i,lnk) in enumerate(net.links)
        l = length(lnk)
        @unpack trac, sval, rval = states[i,t]

        sval += try
            sum(inflows[i,t+1-l,..]) - sum(outflows[i,t,..])
        catch BoundsError
             - sum(outflows[i,t,..])
        end

        rval += try
            sum(outflows[i,t+1-Int(round(l/δ)),..]) - sum(inflows[i,t,..])
        catch BoundsError
            - sum(inflows[i,t,..])
        end

        states[i,t+1] = LinkState(trac+nxttracs[i], sval, rval)
    end
end

# cost computation
revtracs = zeros(Int, nlinks, T)
for i in 1:nlinks
    l = length(link(net, i))
    lastt = Int(floor(tracker(states[i,end])))+1
    revtracs[i,lastt:end] .= T+1
    revtracs[i,lastt-1] = T
    for t in (lastt-2):-1:1
        revtracs[i,t] = (t+l-1) + argfilter(x -> x <= t-1, tracker.(states[i,(t+l):revtracs[i,t+1]]))[end]
    end
end

incosts = zeros(nlinks, T, nsinks, nclasses)
outcosts = zeros(nlinks, T, nsinks, nclasses)
incosts[:,end,..] .= M
outcosts[:,end,..] .= M
for (snkid,snk) in enumerate(snks)
    i = inlinkids(net, snk)[1]
    outcosts[i,:,:,1] .= M
    r = collect(1:T)
    outcosts[i,:,snkid,1] .= clamp.(trgt .- r, 0., Inf) * β .+ clamp.(r .- trgt, 0., Inf) * γ # overwriting
end

for t in (T-1):-1:1
    for src in srcs
        i = outlinkids(net, src)[1]
        incosts[i,t,..] .= incost(outcosts[i,..], states[i,..], revtracs[i,t], revtracs[i,t+1])
    end

    for mrg in mrgs
        ili = inlinkids(net,mrg)
        oli = outlinkids(net,mrg)[1]

        incosts[oli,t,..] .= incost(outcosts[oli,..], states[oli,..], revtracs[oli,t], revtracs[oli,t+1])

        outcosts[ili[1],t,..] .= incosts[oli,t,..]
        outcosts[ili[2],t,..] .= incosts[oli,t,..]
    end

    for div in divs
        ili = inlinkids(net,div)[1]
        oli = outlinkids(net,div)

        # incost computations
        for i in oli
            incosts[i,t,..] .= incost(outcosts[i,..], states[i,..], revtracs[i,t], revtracs[i,t+1])
        end

        sr = [srates[div][i][t,..] for i in 1:2]
        outcosts[ili,t,..] .= sr[1] .* incosts[oli[1],t,..] .+ sr[2] .* incosts[oli[2],t,..]
    end
end

# UE computation
