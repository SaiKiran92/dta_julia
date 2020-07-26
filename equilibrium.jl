
# UE computation
function updatechoices!()#(incosts, dtchoices, srates)
    global incosts, dtchoices, srates
    λ = 1e-3
    for src in srcs
        i = outlinkids(net, src)[1]
        for (snkid,snk) in enumerate(snks)
            for clsid in 1:nclasses
                c = incosts[i,1:Tm,snkid,clsid]
                dtc = dtchoices[src,snk,clsid][1:Tm]
                newdtc = proportionalize(solveqp(c, dtc, λ))[:,1]
                dtchoices[src,snk,clsid][1:Tm] .= newdtc
            end
        end
    end

    for div in divs
        ili = inlinkids(net, div)
        for t in 1:T
            for (snkid,snk) in enumerate(snks)
                for clsid in 1:nclasses
                    sr = [srates[div][i][t, snkid,clsid] for i in 1:2]
                    c = incosts[ili,t,snkid,clsid]
                    newsr = proportionalize(solveqp(c, sr, λ))[:,1]
                    for i in 1:2
                        srates[div][i][t,snkid,clsid] = newsr[i]
                    end
                end
            end
        end
    end
end

function relgap(incosts, dtchoices)
    _num, _den = 0., 0.
    for (srcid,src) in enumerate(srcs)
        i = outlinkids(net, src)[1]
        for (snkid,snk) in enumerate(snks)
            for clsid in 1:nclasses
                m = minimum(incosts[i,1:Tm,snkid,clsid])
                a = sum(incosts[i,1:Tm,snkid,clsid] .* dtchoices[src,snk,clsid][1:Tm])
                _num += (a - m)
                _den += m
                if (srcid == 1) && (snkid == 1) && (clsid == 1)
                    @show a, m
                end
            end
        end
    end
    #@show (_num/_den)
    return _num/_den
end

# equilibrium computation
function optfn!(x, grad, p)
    if length(grad) > 0
        grad .= (x .- p)
    end
    return 0.5*sum((x .- p).^2)
end

function eqcon!(x, grad)
    if length(grad) > 0
        grad .= 1
    end
    return sum(x) - 1
end

function solveqp(c, dtc, λ)
    p = dtc .- λ * c
    x = Variable(length(p))
    problem = minimize(0.5*sumsquares(x - p), sum(x) == 1, x >= 0.)
    solve!(problem, () -> SCS.Optimizer(verbose=false))
    return x.value
end

function proportionalize(x; digits=6)
    x[x .< 0.] .= 0
    x = round.(x, digits=digits)
    x ./= sum(x)
    return x
end
