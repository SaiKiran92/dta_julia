
# UE computation
function updatechoices!()#(ECᵢ, DTC, SR)
    global ECᵢ#, DTC, SR
    newDTC = zero(DTC)
    #newSR = SR
    newSR = Dict(d => Dict(i => zeros(T, nsinks, nclasses) for i in 1:2) for d in divs);
    λ = 1e-4
    for (srcid,src) in enumerate(srcs)
        i = outlinkids(net, src)[1]
        for (snkid,snk) in enumerate(snks)
            for clsid in 1:nclasses
                c = ECᵢ[i,1:Tm,snkid,clsid]
                dtc = DTC[1:Tm, srcid,snkid,clsid]
                newdtc = proportionalize(solveqp(c, dtc, λ))[:,1]
                newDTC[1:Tm, srcid,snkid,clsid] .= newdtc
            end
        end
    end

    for div in divs
        oli = outlinkids(net, div)
        for t in 1:T
            for (snkid,snk) in enumerate(snks)
                for clsid in 1:nclasses
                    sr = [SR[div][i][t, snkid,clsid] for i in 1:2]
                    c = ECᵢ[oli,t,snkid,clsid]
                    newsr = proportionalize(solveqp(c, sr, λ))[:,1]
                    for i in 1:2
                        newSR[div][i][t,snkid,clsid] = newsr[i]
                    end
                end
            end
        end
    end

    return newDTC, newSR
end

function relgap(ECᵢ, DTC)
    _num, _den = 0., 0.
    for (srcid,src) in enumerate(srcs)
        i = outlinkids(net, src)[1]
        for (snkid,snk) in enumerate(snks)
            for clsid in 1:nclasses
                m = minimum(ECᵢ[i,1:Tm,snkid,clsid])
                a = sum(ECᵢ[i,1:Tm,snkid,clsid] .* DTC[1:Tm, srcid,snkid,clsid])
                _num += (a - m)
                _den += m
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
    if (problem.status != MOI.OPTIMAL)
        @show c, dtc, λ
    end
    return x.value
end

function proportionalize(x; digits=6)
    x[x .< 0.] .= 0
    x = round.(x, digits=digits)
    x ./= sum(x)
    return x
end
