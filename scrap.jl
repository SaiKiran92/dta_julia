d, t = 3, 156

ili = inlinkids(net,d)[1]
oli = outlinkids(net,d)
sval = svalue(states[ili,t])
rval = rvalue.(states[oli,t])

l = length(link(net, ili))

sr = [SR[d][ti][t,..] for ti in 1:2]

cf, p, sr, sval, rval, trac, maxt = (@view(CFᵢ[ili,:])), @view(pᵢ[ili,..]), sr, sval, rval, tracker(states[ili,t]), t-l

rng = trac:maxt
fcut = [(mucsum(cf, rng) .* p[indexify(rng),..]) .* expand(sr[i], dims=1) for i in 1:2] # flow cut
sva = valarg(sum(fcut), sval, 1)
rva = min(valarg.(fcut, rval, 1)...)
va = min(sva, rva)

tmp = [squeezesum(fcut[i][0.:va,..], dims=1) for i in 1:2]
newp = [(sum(tmptmp) == 0.) ? ones(size(tmptmp)...)*(1. /prod(size(tmptmp))) : (tmptmp ./ sum(tmptmp)) for tmptmp in tmp]

#newcf = cf[trac] + sum(sum(tmp))
newtrac = cvalarg(cf, cf[trac] + sum(sum(tmp)), trac, maxt)

#newcf = round.(newcf, digits=ROUND_DIGITS)
newtrac = round(newtrac, digits=ROUND_DIGITS)

newcf = round(cf[trac] + sum(cf[trac:newtrac]), digits=ROUND_DIGITS)


newcf, p, newtrac, newtf = dflows((@view(CFᵢ[ili,:])), @view(pᵢ[ili,..]), sr, sval, rval, tracker(states[ili,t]), t-l)

CFᵢ[oli,t] .= round.(newtf .+ CFᵢ[oli,t-1], digits=ROUND_DIGITS)
CFₒ[ili,t] = round(newcf, digits=ROUND_DIGITS)
for (i,li) in enumerate(oli)
    pᵢ[li,t,..] .= p[i]
end
tracs[ili] = newtrac