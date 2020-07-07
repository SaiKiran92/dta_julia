

# reverse states - for computing costs - only trackers needed
rstates = zeros(nlinks, T+1)
tinflows = round.(squeezesum(inflows, dims=(3,4)), digits=ROUND_DIGITS)
toutflows = round.(squeezesum(outflows, dims=(3,4)), digits=ROUND_DIGITS)
for i in 1:nlinks
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
    rstates[i, (lastinidx+1):(maxt-1)] .= collect((lastinidx+1):(maxt-1)) .+ (l-1)

    for t in 2:lastinidx
        ui = floor(rstates[i,t-1] + 1e-9)
        rstates[i,t] = ui + round(argcumval(toutflows[i,ui:end], tinflows[i,t-1], decimal(rstates[i,t-1]), :zero_exclude), digits=ROUND_DIGITS)
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
    println(i, " ", t)
    if rstates[i,t+1] == T
        incosts[i,t,..] .= M
    else
        r = rstates[i,t]:rstates[i,t+1]
        println(r)
        println(size(outcosts[i,r,..]))
        println(size(collect(firstidx(r):lastidx(r))))
        incosts[i,t,..] .= sum((outcosts[i,r,..] .+ collect(firstidx(r):lastidx(r)) .- t) .* safedivide.(toutflows[i,r], Ref(tinflows[i,t])), dims=1)[1,..]
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

        for i in ili
            costupdate!(i,t)
        end

        sr = [srates[div][i][t,..] for i in 1:2]
        outcosts[ili,t,..] .= sr[1] .* incosts[oli[1],t,..] .+ sr[2] .* incosts[oli[2],t,..]
    end
end
