function divflows(infl, sr, r, s, frac)
    sr = [expand(sr[ti], dims=1) for ti in 1:2]
    f = [squeezesum(infl .* sr[ti], dims=(2,3)) for ti in 1:2]

    strac = argcumval(sum(f), s, frac)
    rtrac = min(argcumval.(f, r, frac)...)

    fintrac = min(strac, rtrac)
    tflows = vcat([sum(infl[frac:fintrac,..] .* sr[ti], dims=1) for ti in 1:2]...)
    return (tflows, fintrac - frac)
end

function mrgflows(infl, r, s, fracs)
    tflows = zeros(2, size(infl[1])[2:end]...)
    tracs = [0., 0.]

    function update!(i, cap)
        tmpfrac = argcumval(squeezesum(infl[i], dims=(2,3)), cap, fracs[i])
        tflows[i,..] .= sum(infl[i][fracs[i]:tmpfrac,..], dims=1)[1,..]
        tracs[i] = tmpfrac - fracs[i]
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

    return (tflows, tracs)
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
