a, TFᵢ, TFₒ = nothing, nothing, nothing

struct LinkCostState
    trac
end
tracker(s::LinkCostState) = s.trac

function computecosts()
    global a, net, Fᵢ, Fₒ, TFᵢ, TFₒ, CFᵢ, CFₒ, states, ECᵢ, ECₒ

    rstates = Matrix{Union{Nothing, LinkCostState}}(fill(nothing, nlinks, T+1)); #zeros(Int, nlinks, T+1)
    for i in 1:nlinks
        l = length(link(net, i))
        # a link of length 2 - min. time for outflow is 3 => at the beginning of time 1, (-1, 0.)
        rstates[i,(T+1-l):end] .= LinkCostState(l+1)
    end

    """
    for i in 1:nlinks
        τ, newτ = 1, NaN
        l = length(link(net, i))
        t = a[i,τ] = l+1
        while t <= T
            newτ = tracker(states[i,t])
            if newτ > τ
                if newτ > τ+1
                    a[i,(τ+1):(newτ-1)] .= t-1
                end
                τ = newτ
                a[i,τ] = t
            end
            t += 1
        end

        a[i,(newτ+1):end] .= T+1
    end"""

    TFᵢ = squeezesum(Fᵢ, dims=(3,4))
    TFₒ = squeezesum(Fₒ, dims=(3,4))

    ECᵢ, ECₒ = zero(Fᵢ), zero(Fₒ);
    ECᵢ[:,end,..] .= M
    ECₒ[:,end,..] .= M
    for (snkid,snk) in enumerate(snks)
        i = inlinkids(net, snk)[1]
        ECₒ[i,:,:,1] .= M
        r = collect(1:T)
        ECₒ[i,:,snkid,1] .= clamp.(trgt .- r, 0., Inf) * β .+ clamp.(r .- trgt, 0., Inf) * γ # overwriting
    end

    function costupdate!(i, τ)
        y = tracker(rstates[i,τ+1])
        if (τ > 1) & (CFᵢ[i,τ-1] > CFₒ[i,end])
            ECᵢ[i,τ,..] .= M
        else
            
        end
    end

    function costupdate!(i, τ)
        x, y = a[i,τ], a[i,τ+1]
        if y == T+1
            ECᵢ[i,τ,..] .= M
        elseif x == y
            ECᵢ[i,τ,..] .= (x - τ) .+ ECₒ[i,x,..]
        else
            tf = max(0., CFᵢ[i,τ] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.))
            ECᵢ[i,τ,..] .= ((x - τ) .+ ECₒ[i,x,..]) * safedivide(max(0., CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)), tf)
            ECᵢ[i,τ,..] .+= ((y - τ) .+ ECₒ[i,y,..]) * safedivide(max(0., CFₒ[i,τ] - CFᵢ[i,y-1]), tf, 0.)
            #ECᵢ[i,τ,..] .= ((x - τ) .+ ECₒ[i,x,..]) * safedivide(max(0., CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)), TFᵢ[i,τ])
            #ECᵢ[i,τ,..] .+= ((y - τ) .+ ECₒ[i,y,..]) * safedivide(max(0., CFₒ[i,τ] - CFᵢ[i,y-1]), TFᵢ[i,τ], 0.)
            if y > x+1
                ECᵢ[i,τ,..] .+= squeezesum(((((x+1):(y-1)) .- τ) .+ ECₒ[i,(x+1):(y-1),..]) .* expand(safedivide.(TFₒ[i,(x+1):(y-1)], Ref(TFᵢ[i,τ])), dims=(2,3)), dims=1)
            end
        end

        ECᵢ[i,τ,..] .= round.(ECᵢ[i,τ,..], digits=ROUND_DIGITS)
    end

    for τ in T:-1:1
        for src in srcs
            i = outlinkids(net, src)[1]
            costupdate!(i,τ)
        end

        for mrg in mrgs
            ili = inlinkids(net,mrg)
            oli = outlinkids(net,mrg)[1]

            costupdate!(oli,τ)

            ECₒ[ili[1],τ,..] .= ECᵢ[oli,τ,..]
            ECₒ[ili[2],τ,..] .= ECᵢ[oli,τ,..]
        end

        for div in divs
            ili = inlinkids(net,div)[1]
            oli = outlinkids(net,div)

            for li in oli
                costupdate!(li,τ)
            end

            sr = [SR[div][i][τ,..] for i in 1:2]
            ECₒ[ili,τ,..] .= sr[1] .* ECᵢ[oli[1],τ,..] .+ sr[2] .* ECᵢ[oli[2],τ,..]
        end
    end
end
