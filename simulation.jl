function simulate()#(DTC, SR)
    global net, Fᵢ, Fₒ, CFᵢ, CFₒ, states

    # initialization
    Fᵢ = zeros(nlinks, T, nsinks, nclasses); # inflows
    Fₒ = zeros(nlinks, T, nsinks, nclasses); # outflows
    CFᵢ = zeros(nlinks, T); # cumulative inflows
    CFₒ = zeros(nlinks, T); # cumulative outflows
    for src in sources(net)
        li = outlinkids(net,src)[1]
        for (snkno,snk) in enumerate(sinks(net))
            for c in 1:nclasses
                Fᵢ[li,:,snkno,c] .= round.(DTC[src,snk,c] * trips[src,snk,c], digits=ROUND_DIGITS)
            end
        end
    end
    CFᵢ .= cumsum(squeezesum(Fᵢ, dims=(3,4)), dims=2)

    srclinkids = [outlinkids(net,src)[1] for src in srcs];
    states = Matrix{Union{Nothing, LinkState}}(fill(nothing, nlinks, T+1));
    for i in 1:nlinks
        l = length(link(net, i))
        # a link of length 2 - min. time for outflow is 3 => at the beginning of time 1, (-1, 0.)
        for j in 1:l
            states[i,j] = LinkState(-l+j, 0., 0., (i in srclinkids) ? Inf : N*l)
        end
    end

    for t in 1:T
        A = ones(Int, nlinks);
        fracs = zeros(nlinks);
        for snk in snks
            i = inlinkids(net,snk)[1]
            lnk = link(net, i)
            l = length(lnk)

            Fₒ[i,t,:,:] .= (t > l) ? Fᵢ[i,t-l,:,:] : 0.

            A[i] = tracker(states[i,t]) + 1
        end

        for mrg in mrgs
            # sending and receiving capacities
            ili = inlinkids(net,mrg)
            oli = outlinkids(net,mrg)[1]
            s = svalue.(states[ili,t])
            r = rvalue(states[oli,t])

            # calculating flows
            trkr = tracker.(states[ili,t])
            ls = length.(link.(Ref(net), ili))
            f = [(t > l) ? Fᵢ[li, trkr[i]:(t - l),..] : zeros(1, size(Fᵢ)[3:end]...)  for (i,(li,l)) in enumerate(zip(ili, ls))]
            fₐ, va = mflows(f, r, s, fracmoved.(states[ili,t]))

            va = round.(va, digits=ROUND_DIGITS)

            for (i,li) in enumerate(ili)
                fₐ[i] .= round.(fₐ[i], digits=ROUND_DIGITS)
                Fₒ[li,t,..] .= fₐ[i]
            end
            Fᵢ[oli,t,..] .= sum(fₐ, dims=1)[1,..]

            for (i,li) in enumerate(ili)
                CFₒ[li,t] = round(sum(Fₒ[li,t,..]) + ((t > 1) ? CFₒ[li,t-1] : 0.), digits=ROUND_DIGITS)
                A[li], fracs[li] = cumvalarg((@view CFᵢ[li,:]), CFₒ[li,t], tracker(states[li,t]), t - ls[i])
            end

"""
            A[ili] .= (trkr .+ Int.(floor.(va)))
            fracs[ili] .= decimal.(va)"""
        end

        for div in divs
            # sending and receiving capacities
            ili = inlinkids(net,div)[1]
            oli = outlinkids(net,div)
            s = svalue(states[ili,t])
            r = rvalue.(states[oli,t])

            # calculating flows
            sr = [SR[div][ti][t,..] for ti in 1:2]
            trkr = tracker(states[ili,t])
            l = length(link(net, ili))
            fₐ, va = dflows(((t > l) ? Fᵢ[ili, trkr:(t - l),..] : zeros(1, size(Fᵢ)[3:end]...)), sr, r, s, fracmoved(states[ili,t]))

            va = round(va, digits=ROUND_DIGITS)

            for (i,li) in enumerate(oli)
                fₐ[i] .= round.(fₐ[i], digits=ROUND_DIGITS)
                Fᵢ[li,t,..] .= fₐ[i]
            end
            Fₒ[ili,t,..] .= sum(fₐ, dims=1)[1,..]

            CFₒ[ili,t] = round(sum(Fₒ[ili,t,..]) + ((t > 1) ? CFₒ[ili,t-1] : 0.), digits=ROUND_DIGITS)
            A[ili], fracs[ili] = cumvalarg((@view CFᵢ[ili,:]), CFₒ[ili,t], tracker(states[ili,t]), t - l)

"""
            A[ili] = (trkr + Int(floor(va)))
            fracs[ili] = decimal(va)"""
        end

        # update cumulative flows
        #CFᵢ[:,t] .= ((t > 1) ? CFᵢ[:,t-1] : 0.) .+ squeezesum(Fᵢ[:,t,..], dims=(2,3))
        CFₒ[:,t] .= ((t > 1) ? CFₒ[:,t-1] : 0.) .+ squeezesum(Fₒ[:,t,..], dims=(2,3))

        # update link states
        for (i,lnk) in enumerate(net.links)
            l = length(lnk)
            sval, rval = - CFₒ[i,t],  -CFᵢ[i,t] + ((i in srclinkids) ? Inf : N*l)
            sval = min(Q, sval + ((t >= l) ? CFᵢ[i,t+1-l] : 0.))
            rval = min(Q, rval + ((t >= Int(round(l/δ))) ? CFₒ[i,t+1-Int(round(l/δ))] : 0.))
            """if t >= l
                sval = min(Q, sval + CFᵢ[i,t+1-l])
                if t >= Int(round(l/δ))
                    rval = min(Q, rval + CFₒ[i,t+1-Int(round(l/δ))])
                end
            end"""
            sval = round(sval, digits=ROUND_DIGITS)
            rval = round(rval, digits=ROUND_DIGITS)

            states[i,t+1] = LinkState(A[i], fracs[i], sval, rval)
        end
    end
end
