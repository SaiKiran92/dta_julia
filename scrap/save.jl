
        sval += try
            Nᵢ[i,t+1-l] - Nₒ[i,t]
        catch BoundsError
            -Nₒ[i,t]
        end

        rval += try
            Nₒ[i,t+1-Int(round(l/δ))] - Nᵢ[i,t]
        catch BoundsError
            -Nᵢ[i,t]
        end

        sval = round(sval, digits=ROUND_DIGITS)
        rval = round(rval, digits=ROUND_DIGITS)
        if (sval < 0.) || (rval < 0.)
            sval = (sval < 0.) ? 0. : sval
            rval = (rval < 0.) ? 0. : rval
        end

        """
        function argcut(r::Array, j, i)
            c = (r isa Vector) ? [0.] : zero(r[1,..])
            inti, intj = Int(floor(i))+1, Int(ceil(j))
            if (intj == j)
                intj -= 1
                dec = 1.
            else
                dec = decimal(j)
            end

            @show (i, " ", j)
            @show (inti, " ", intj)

            if inti == intj
                c .= r[inti,..] .* (j - i)
            else
                c .= r[inti,..] .* (inti - i)
                c .+= r[intj,..] .* dec
                if intj > inti + 1
                    c .+= squeezesum(r[(inti+1):(intj-1),..], dims=1)
                end
            end

            return c
        end
        """






        #function computecosts()
        #global rstates, tinflows, toutflows, CFᵢ, CFₒ, Cᵢ, Cₒ
        # reverse states - for computing costs - only trackers needed

        rstates = zeros(nlinks, T+1)
        tinflows = round.(squeezesum(inflows, dims=(3,4)), digits=ROUND_DIGITS)
        toutflows = round.(squeezesum(outflows, dims=(3,4)), digits=ROUND_DIGITS)
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
                rstates[i,t] = ui + round(argcumval2(CFₒ[i,ui:lastoutidx] .- CFₒ[i,Int(ui)], max(0.,CFᵢ[i,t-1] .- CFₒ[i,Int(ui)]), 0., :zero_exclude), digits=ROUND_DIGITS)
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

        function costupdate!(i,τ)
            if rstates[i,τ+1] == T
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

            Cᵢ[i,τ,..] .= round.(Cᵢ[i,τ,..], digits=ROUND_DIGITS)
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

        #return (rstates, Cᵢ, Cₒ)
        #end


        if (srcid == 1) && (snkid == 1) && (clsid == 1)
            @show a, m
        end
