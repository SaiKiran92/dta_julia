using DataFrames, CSV
using DataStructures
using Parameters
#using JuMP, Ipopt
using Convex, SCS
using Plots

using EllipsisNotation

ROUND_DIGITS = 6
ATOL = 0.1^(ROUND_DIGITS)

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

function initchoices()
    # choice initialization
    dtchoices = Dict(k => zeros(T) for k in keys(trips));
    for k in keys(trips)
        #dtchoices[k][1:Tm] .= rand(Tm)
        #dtchoices[k][1:Tm] ./= sum(dtchoices[k][1:Tm])
        #dtchoices[k][1:Tm] .= 1/Tm
        dtchoices[k][1] = 1.
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
    return (dtchoices, srates)
end

inflows, outflows, incosts, outcosts, states, revtracs = nothing, nothing, nothing, nothing, nothing, nothing
function simulate()#(dtchoices, srates)
    global inflows, outflows, states
    # initialization
    inflows = zeros(nlinks, T, nsinks, nclasses);
    outflows = zeros(nlinks, T, nsinks, nclasses);
    for src in sources(net)
        for (snkno,snk) in enumerate(sinks(net))
            for c in 1:nclasses
                inflows[outlinkids(net,src)[1],:,snkno,c] .= round.(dtchoices[src,snk,c] * trips[src,snk,c], digits=ROUND_DIGITS)
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

            sval = round(sval, digits=ROUND_DIGITS)
            rval = round(rval, digits=ROUND_DIGITS)
            #nxttracs[i] = round(nxttracs[i], digits=ROUND_DIGITS)
            if (sval < 0.) || (rval < 0.)
                sval = (sval < 0.) ? 0. : sval
                rval = (rval < 0.) ? 0. : rval
            end

            states[i,t+1] = LinkState(trac+nxttracs[i], sval, rval)
        end
    end

    return (inflows, outflows, states)
end

function computecosts()
    global rstates, tinflows, toutflows, cinflows, coutflows, incosts, outcosts
    # reverse states - for computing costs - only trackers needed
    rstates = zeros(nlinks, T+1)
    tinflows = round.(squeezesum(inflows, dims=(3,4)), digits=ROUND_DIGITS)
    toutflows = round.(squeezesum(outflows, dims=(3,4)), digits=ROUND_DIGITS)
    cinflows = round.(cumsum(tinflows, dims=2), digits=ROUND_DIGITS)
    coutflows = round.(cumsum(toutflows, dims=2), digits=ROUND_DIGITS)
    for i in 1:nlinks
        #println(i)
        l = length(link(net, i))
        maxt = Int(ceil(tracker(states[i,end])+1e-13))
        rstates[i,maxt:end] .= T
        rstates[i,1] = l

        # adjustment for roundoffs
        try
            global lastinidx = argfilter(x -> x > 0., tinflows[i,:])[end]
        catch BoundsError
            # no flow into the link
            rstates[i,1:(maxt-1)] .= collect(1:(maxt-1)) .+ (l-1)
            continue
        end
        lastoutidx = argfilter(x -> x > 0., toutflows[i,:])[end]
        rstates[i,(lastinidx+1):(lastoutidx-l+1)] .= lastoutidx
        rstates[i,(lastoutidx-l+2):(maxt-1)] .= collect((lastoutidx-l+2):(maxt-1)) .+ (l-1)

        for t in 2:lastinidx
            ui = floor(rstates[i,t-1] + (0.1)^(ROUND_DIGITS-1))
            #println(i, " ", t, " ", ui)
            rstates[i,t] = ui + round(argcumval2(coutflows[i,ui:lastoutidx] .- coutflows[i,Int(ui)], max(0.,cinflows[i,t-1] .- coutflows[i,Int(ui)]), 0., :zero_exclude), digits=ROUND_DIGITS)
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

    function costupdate!(i,t)
        if rstates[i,t+1] == T
            incosts[i,t,..] .= M
        else
            r = rstates[i,t]:rstates[i,t+1]
            tmpr = collect(firstidx(r):lastidx(r))
            if size(tmpr) == (0,) # when rstates[i,t] == rstates[i,t+1] == Int(rstates[i,t])
                tmpr = firstidx(r)
                incosts[i,t,..] .= outcosts[i,tmpr,..] .+ tmpr .- t
            else
                x = safedivide.(toutflows[i,r], Ref(tinflows[i,t])) # floating point issues
                incosts[i,t,..] .= sum((outcosts[i,tmpr,..] .+ tmpr .- t) .* (x ./ sum(x)), dims=1)[1,..]
            end
        end
    end

    for t in (T-1):-1:1
        for src in srcs
            i = outlinkids(net, src)[1]
            costupdate!(i,t)
        end

        for mrg in mrgs
            ili = inlinkids(net,mrg)
            oli = outlinkids(net,mrg)[1]

            costupdate!(oli,t)

            outcosts[ili[1],t,..] .= incosts[oli,t,..]
            outcosts[ili[2],t,..] .= incosts[oli,t,..]
        end

        for div in divs
            ili = inlinkids(net,div)[1]
            oli = outlinkids(net,div)

            for i in oli
                costupdate!(i,t)
            end

            sr = [srates[div][i][t,..] for i in 1:2]
            outcosts[ili,t,..] .= sr[1] .* incosts[oli[1],t,..] .+ sr[2] .* incosts[oli[2],t,..]
        end


    end

    return (rstates, incosts, outcosts)
end

# UE computation
function updatechoices!()#(incosts, dtchoices, srates)
    global incosts, dtchoices, srates
    λ = 1e-3
    for src in srcs
        i = outlinkids(net, src)[1]
        for (snkid,snk) in enumerate(snks)
            for clsid in 1:nclasses
                c = incosts[i,1:Tm,snkid,clsid]
                dtc = dtchoices[src,snk,clsid][1:Tm]
                newdtc = proportionalize(solveqp(c, dtc, λ))[:,1]
                dtchoices[src,snk,clsid][1:Tm] .= newdtc
            end
        end
    end

    for div in divs
        ili = inlinkids(net, div)
        for t in 1:T
            for (snkid,snk) in enumerate(snks)
                for clsid in 1:nclasses
                    sr = [srates[div][i][t, snkid,clsid] for i in 1:2]
                    c = incosts[ili,t,snkid,clsid]
                    newsr = proportionalize(solveqp(c, sr, λ))[:,1]
                    for i in 1:2
                        srates[div][i][t,snkid,clsid] = newsr[i]
                    end
                end
            end
        end
    end
end

function relgap(incosts, dtchoices)
    _num, _den = 0., 0.
    for (srcid,src) in enumerate(srcs)
        i = outlinkids(net, src)[1]
        for (snkid,snk) in enumerate(snks)
            for clsid in 1:nclasses
                m = minimum(incosts[i,1:Tm,snkid,clsid])
                a = sum(incosts[i,1:Tm,snkid,clsid] .* dtchoices[src,snk,clsid][1:Tm])
                _num += (a - m)
                _den += m
                if (srcid == 1) && (snkid == 1) && (clsid == 1)
                    @show a, m
                end
            end
        end
    end
    #@show (_num/_den)
    return _num/_den
end

dtchoices, srates = initchoices()
for i in 1:30
    global dtchoices, srates, incosts, outcosts
    #@show i
    inflows, outflows, states = simulate()#(dtchoices, srates)
    rstates, incosts, outcosts = computecosts()#(states)
    relgap(incosts, dtchoices)

    @assert approxpositive(inflows)# all((inflows .>= 0.) .| (inflows .≈ 0.))
    @assert approxpositive(outflows)# all(outflows .>= 0.)
    @assert approxpositive(incosts)# all(incosts .>= 0.)
    @assert approxpositive(outcosts)# all(outcosts .>= 0.)

    updatechoices!()#(incosts, dtchoices, srates)
end

dtchoices, srates = initchoices()

inflows, outflows, states = simulate()#(dtchoices, srates)
rstates, incosts, outcosts = computecosts()#(states)
relgap(incosts, dtchoices)
plot(dtchoices[1,19,1])

updatechoices!()

relgap(incosts, dtchoices)
