#function simulate()#(DTC, SR)
#global net, CFᵢ, CFₒ, states

# initialization
pᵢ = zeros(nlinks, T, nsinks, nclasses); # inflow proportions
CFᵢ = zeros(nlinks, T); # cumulative inflows
CFₒ = zeros(nlinks, T); # cumulative outflows
for (rno,r) in enumerate(sources(net))
    tmp = @view DTC[:,rno,..]
    pᵢ[rno,..] .= tmp./max.(Ref(1e-16), sum(tmp, dims=(2,3)))
    CFᵢ[rno,:] .= cumsum(squeezesum(tmp .* expand(trips[rno,..], dims=1), dims=(2,3)))
end

rlinkids = [outlinkids(net, r)[1] for r in srcs]
states = Matrix{Union{Nothing, LinkState}}(fill(nothing, nlinks, T+1))
for i in 1:nlinks
    l = length(link(net, i))
    for j in 1:l
        states[i,j] = LinkState(j-1-l, 0., (i in rlinkids) ? Inf : Q)
    end
end

for t in 1:T
    tracs = zeros(nlinks)

    for k in snks
        li = inlinkids(net,k)[1]
        l = length(link(net, li))

        CFₒ[li,t] = (t > l) ? CFᵢ[li,t-l] : 0.
        tracs[li] = tracker(states[li,t]) + 1
    end

    for m in mrgs
        #@show t, m
        # sending and receiving capacities
        ili = inlinkids(net,m)
        oli = outlinkids(net,m)[1]
        sval = svalue.(states[ili,t])
        rval = rvalue(states[oli,t])

        ls = length.(link.(Ref(net), ili))

        # cuminflows, props, sval, rval, tracs, maxt
        newcf, p, newtracs = mflows([@view(CFᵢ[ili[1],:]), @view(CFᵢ[ili[2],:])], [@view(pᵢ[ili[1],..]), @view(pᵢ[ili[2],..])], sval, rval, tracker.(states[ili,t]), t.-ls)

        CFᵢ[oli,t] = round(sum(newcf), digits=ROUND_DIGITS)
        CFₒ[ili,t] .= round.(newcf, digits=ROUND_DIGITS)
        pᵢ[oli,t,..] .= p
        tracs[ili] .= newtracs
    end

    for d in divs
        #@show t, d
        # sending and receiving capacities
        ili = inlinkids(net,d)[1]
        oli = outlinkids(net,d)
        sval = svalue(states[ili,t])
        rval = rvalue.(states[oli,t])

        l = length(link(net, ili))

        # cuminflows, props, sr, sval, rval, trac, maxt
        if t > l
            sr = [SR[d][ti][t,..] for ti in 1:2]
            newcf, p, newtrac, newtf = dflows((@view(CFᵢ[ili,:])), @view(pᵢ[ili,..]), sr, sval, rval, tracker(states[ili,t]), t-l)

            CFᵢ[oli,t] .= round.(newtf .+ CFᵢ[oli,t-1], digits=ROUND_DIGITS)
            CFₒ[ili,t] = round(newcf, digits=ROUND_DIGITS)
            for (i,li) in enumerate(oli)
                pᵢ[li,t,..] .= p[i]
            end
            tracs[ili] = newtrac
        end
    end

    for (li,lnk) in enumerate(net.links)
        l = length(lnk)
        if t >= l
            sval = round(min(Q, CFᵢ[li, t+1-l] - CFₒ[li,t]), digits=ROUND_DIGITS)
            rval = if li in rlinkids
                Inf
            else
                round(min(Q, ((t >= Int(round(l/δ))) ? CFₒ[li,t+1-Int(round(l/δ))] : 0.) - CFᵢ[li,t] + N*l), digits=ROUND_DIGITS)
            end

            sval = max.(Ref(0.), sval)
            rval = max(0., rval)
            
            states[li,t+1] = LinkState(tracs[li], sval, rval)
        end
    end
end

#end
