
function incost(ocosts, st, rtraca, rtracb)
    c = zeros(size(ocosts)[2:end])

    if rtracb == T+1
        c .= M
    else
        t = floor(tracker(st[rtracb]) + 1e-10)
        if rtraca == rtracb
            c .= ocosts[rtraca,..] .+ α * (rtraca - t)
        else
            fracs = tracker.(st[rtraca:rtracb]) .- (t-1)
            fracs[1] = clamp(fracs[1], 0., 1.)
            fracs[end] = clamp(fracs[end], 0., 1.)
            fracs[2:end] .-= fracs[1:(end-1)]

            for (f,rt) in zip(fracs, rtraca:rtracb)
                c .+= f * (α * (rt .- t) .+ ocosts[rt,..])
            end
        end
    end

    return round.(c, digits=ROUND_DIGITS)
end

function computecosts()
    global rstates, tinflows, toutflows, Nᵢ, Nₒ, Cᵢ, Cₒ
    # reverse states - for computing costs - only trackers needed
    rstates = zeros(nlinks, T+1)
    tinflows = round.(squeezesum(inflows, dims=(3,4)), digits=ROUND_DIGITS)
    toutflows = round.(squeezesum(outflows, dims=(3,4)), digits=ROUND_DIGITS)
    Nᵢ = round.(cumsum(tinflows, dims=2), digits=ROUND_DIGITS)
    Nₒ = round.(cumsum(toutflows, dims=2), digits=ROUND_DIGITS)
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
            rstates[i,t] = ui + round(argcumval2(Nₒ[i,ui:lastoutidx] .- Nₒ[i,Int(ui)], max(0.,Nᵢ[i,t-1] .- Nₒ[i,Int(ui)]), 0., :zero_exclude), digits=ROUND_DIGITS)
        end
    end

    Cᵢ = zeros(nlinks, T, nsinks, nclasses)
    Cₒ = zeros(nlinks, T, nsinks, nclasses)
    Cᵢ[:,end,..] .= M
    Cₒ[:,end,..] .= M
    for (snkid,snk) in enumerate(snks)
        i = inlinkids(net, snk)[1]
        Cₒ[i,:,:,1] .= M
        r = collect(1:T)
        Cₒ[i,:,snkid,1] .= clamp.(trgt .- r, 0., Inf) * β .+ clamp.(r .- trgt, 0., Inf) * γ # overwriting
    end

    function costupdate!(i,t)
        if rstates[i,t+1] == T
            Cᵢ[i,t,..] .= M
        else
            r = rstates[i,t]:rstates[i,t+1]
            tmpr = collect(firstidx(r):lastidx(r))
            if size(tmpr) == (0,) # when rstates[i,t] == rstates[i,t+1] == Int(rstates[i,t])
                tmpr = firstidx(r)
                Cᵢ[i,t,..] .= Cₒ[i,tmpr,..] .+ tmpr .- t
            else
                x = safedivide.(toutflows[i,r], Ref(tinflows[i,t])) # floating point issues
                Cᵢ[i,t,..] .= sum((Cₒ[i,tmpr,..] .+ tmpr .- t) .* (x ./ sum(x)), dims=1)[1,..]
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

            Cₒ[ili[1],t,..] .= Cᵢ[oli,t,..]
            Cₒ[ili[2],t,..] .= Cᵢ[oli,t,..]
        end

        for div in divs
            ili = inlinkids(net,div)[1]
            oli = outlinkids(net,div)

            for i in oli
                costupdate!(i,t)
            end

            sr = [SR[div][i][t,..] for i in 1:2]
            Cₒ[ili,t,..] .= sr[1] .* Cᵢ[oli[1],t,..] .+ sr[2] .* Cᵢ[oli[2],t,..]
        end


    end

    return (rstates, Cᵢ, Cₒ)
end
