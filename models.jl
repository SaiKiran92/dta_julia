function mflows(cf, p, sval, rval, tracs, maxt; digits=ROUND_DIGITS)
    newtracs = [0., 0.]

    function update!(i, f)
        #newcf[i] = cf[i][tracs[i]] + f
        newtracs[i] = cvalarg(cf[i], cf[i][tracs[i]] + f, tracs[i], maxt[i])
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
    newtracs = round.(newtracs, digits=ROUND_DIGITS)

    newcf = round.([cf[i][newtracs[i]] for i in 1:2], digits=ROUND_DIGITS)

    newp = zeros(size(p[1])[2:end]...)
    for i in 1:2
        rng = tracs[i]:newtracs[i]
        if size(rng) > 0.
            newp .+= squeezesum(mucsum(cf[i], rng) .* p[i][indexify(rng),..], dims=1)
        end
    end
    if sum(newp) == 0.
        newp .= 1. /prod(size(newp))
    else
        newp ./= sum(newp)
    end

    return (newcf, newp, newtracs)
end

function dflows(cf, p, sr, sval, rval, trac, maxt; digits=ROUND_DIGITS)
    rng = trac:maxt
    fcut = [(mucsum(cf, rng) .* p[indexify(rng),..]) .* expand(sr[i], dims=1) for i in 1:2] # flow cut
    sva = valarg(sum(fcut), sval, 1)
    rva = min(valarg.(fcut, rval, 1)...)
    va = min(sva, rva)

    tmp = [squeezesum(fcut[i][0.:va,..], dims=1) for i in 1:2]
    newp = [(sum(tmptmp) == 0.) ? ones(size(tmptmp)...)*(1. /prod(size(tmptmp))) : (tmptmp ./ sum(tmptmp)) for tmptmp in tmp]

    newtrac = cvalarg(cf, cf[trac] + sum(sum(tmp)), trac, maxt)
    newtrac = round(newtrac, digits=ROUND_DIGITS)

    newcf = round(cf[newtrac], digits=ROUND_DIGITS)

    return (newcf, newp, newtrac, sum.(tmp))#[sum(tmptmp) for tmptmp in tmp])
end
