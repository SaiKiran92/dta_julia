
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

function cumvalarg(cumr, v, starti, stopi)
    if starti <= 0
        return (starti+1, 0.)
    end
    i = starti
    while i <= stopi
        if v < cumr[i]
            tmp1 = (i > 1) ? cumr[i-1] : 0.
            tmp2 = (i > 1) ? max(0., cumr[i] - cumr[i-1]) : cumr[i]
            return (i, safedivide((v - tmp1), tmp2))
        end
        i += 1
    end
    return (stopi+1, 0.)
end



function dflows(fᵢ, sr, r, s, m)
    mval = [m*sum(fᵢ[1,..] .* sr[i]) for i in 1:2]
    f = []
    for (k,v) in pairs(sr)
        push!(f, squeezesum(fᵢ .* expand(v, dims=1), dims=(2,3)))
    end

    sva = valarg(sum(f), s, sum(mval))
    rva = min(valarg.(f, r, mval)...)
    va = min(sva, rva)

    fₐ = []
    for (i,(k,v)) in enumerate(pairs(sr))
        push!(fₐ, argcut(fᵢ, va, m) .* v)
    end

    return (fₐ, va)
end



function mflows(cf, p, sval, rval, tracs, maxt)
    cfo = [0., 0.]
    pₐ = zeros(2, size(cf[1])[2:end]...) #[zeros(size(cf[1])[2:end]) for _ in 1:2]
    newtracs = [0., 0.]

    function update!(i, cap)
        cfo[i] = cf[tracs[i]] + cap
        tmparg .= cvalarg(cf, cval, fromt, maxt)
        idxs = indexify(tracs[i]:maxt)
        tmpf = squeezesum(expand(cfo[idxs], dims=(2,3)) .* pₐ[idxs,..], dims=(2,3))
        pₐ[i] .= tmpf ./ sum(tmpf)
    end

    if rval >= sum(sval)
        update!(1, sval[1]) #free
        update!(2, sval[2]) #free
    elseif sval[1] < 0.5*rval
        update!(1, sval[1]) #free
        update!(2, rval-sval[1]) #cong
    elseif sval[2] < 0.5*rval
        update!(1, rval-sval[2]) #free
        update!(2, sval[2]) #cong
    else
        update!(1, 0.5*rval) #cong
        update!(2, 0.5*rval) #cong
    end

    return (cfo, pₐ, newtracs)
end


function mflows(fᵢ, r, s, m)
    fₐ = [zeros(size(fᵢ[1])[2:end]...) for i in 1:2]
    va = [0., 0.]

    function update!(i, cap)
        va[i] = valarg(fᵢ[i], cap, m[i]*sum(fᵢ[i][1,..]))
        fₐ[i] .= argcut(fᵢ[i], va[i], m[i])
    end

    if r >= sum(s)
        update!(1, s[1]) #free
        update!(2, s[2]) #free
    elseif s[1] < 0.5*r
        update!(1, s[1]) #free
        update!(2, r-s[1]) #cong
    elseif s[2] < 0.5*r
        update!(1, r-s[2]) #free
        update!(2, s[2]) #cong
    else
        update!(1, 0.5*r) #cong
        update!(2, 0.5*r) #cong
    end

    return (fₐ, va)
end



for t in 1:T
    for mrg in mrgs
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
            tracs[li], fracs[li] = cumvalarg((@view CFᵢ[li,:]), CFₒ[li,t], tracker(states[li,t]), t - ls[i])
        end

"""
        tracs[ili] .= (trkr .+ Int.(floor.(va)))
        fracs[ili] .= decimal.(va)"""
    end

    for div in divs
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
        tracs[ili], fracs[ili] = cumvalarg((@view CFᵢ[ili,:]), CFₒ[ili,t], tracker(states[ili,t]), t - l)

"""
        tracs[ili] = (trkr + Int(floor(va)))
        fracs[ili] = decimal(va)"""
    end

    # update cumulative flows
    #CFᵢ[:,t] .= ((t > 1) ? CFᵢ[:,t-1] : 0.) .+ squeezesum(Fᵢ[:,t,..], dims=(2,3))
    CFₒ[:,t] .= ((t > 1) ? CFₒ[:,t-1] : 0.) .+ squeezesum(Fₒ[:,t,..], dims=(2,3))
end
#end



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


"""
function cvalarg(r::AbstractArray, v, starti=0, maxi=size(r)[1])
    i = starti
    while i < maxi
        totr = r[i+1] - r[i]
        if v < totr
            return i + ((v - r[i-1])/totr)
        end
        i += 1
    end

    return maxi
end"""