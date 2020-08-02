function computecosts()
    global ECᵢ, ECₒ, rtracs

    ECᵢ, ECₒ = zero(pᵢ), zero(pᵢ)
    ECᵢ[:,end,..] .= M
    ECₒ[:,end,..] .= M
    for (snkid,snk) in enumerate(snks)
        i = inlinkids(net, snk)[1]
        ECₒ[i,:,:,1] .= M
        r = collect(1:T)
        ECₒ[i,:,snkid,1] .= clamp.(trgt .- r, 0., Inf) * β .+ clamp.(r .- trgt, 0., Inf) * γ # overwriting
    end

    rtracs = zeros(nlinks, T+1)
    for i in 1:nlinks
        l = length(link(net, i))
        rtracs[i,(T-l+1):end] .= T
        for τ in (T-l):-1:2
            if CFᵢ[i,τ] > CFₒ[i,end]
                rtracs[i,τ] = T
            else
                rtracs[i,τ] = rcvalarg(@view(CFₒ[i,:]), CFᵢ[i,τ-1], rtracs[i,τ+1], τ+l-1)
            end
        end
        rtracs[i,1] = l
    end

    function costupdate!(i, τ)
        # goal: set ECᵢ[i,τ,..]
        if rtracs[i,τ+1] == T
            ECᵢ[i,τ,..] .= M
        else
            #@show i, τ
            rng = rtracs[i,τ]:rtracs[i,τ+1]
            inds = indexify(rng)
            #@show i, τ, rng, inds
            if size(rng) == 0.
                ind = lastidx(rng)
                ECᵢ[i,τ,..] .= (ind - τ) .+ ECₒ[i,ind,..]
            else
                tmp = mucsum(@view(CFₒ[i,:]), rng)
                tmp[end] += ϵ;
                tf = (CFᵢ[i,τ] - ((τ == 1) ? 0. : CFᵢ[i,τ-1])) + ϵ
                #@show rng, inds
                #@show tmp, tf, (tmp ./ tf), (inds .- τ)
                ECᵢ[i,τ,..] .= squeezesum((tmp ./ tf) .* ((inds .- τ) .+ ECₒ[i,inds,..]), dims=1)
                #@show i, τ, ECᵢ[i,τ,:,1]
            end
        end
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

    return ECᵢ, ECₒ, rtracs
end
