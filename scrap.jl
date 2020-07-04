using NLopt

function myfunc(x::Vector, grad::Vector)
    if length(grad) > 0
        grad[1] = 0
        grad[2] = 0.5/sqrt(x[2])
    end
    return sqrt(x[2])
end

function myconstraint(x::Vector, grad::Vector, a, b)
    if length(grad) > 0
        grad[1] = 3a * (a*x[1] + b)^2
        grad[2] = -1
    end
    (a*x[1] + b)^3 - x[2]
end

opt = Opt(:LD_MMA, 2)
opt.lower_bounds = [-Inf, 0.]
opt.xtol_rel = 1e-4

opt.min_objective = myfunc
inequality_constraint!(opt, (x,g) -> myconstraint(x,g,2,0), 1e-8)
inequality_constraint!(opt, (x,g) -> myconstraint(x,g,-1,1), 1e-8)

(minf,minx,ret) = optimize(opt, [1.234, 5.678])
numevals = opt.numevals # the number of function evaluations
println("got $minf at $minx after $numevals iterations (returned $ret)")


p = dtc .- λ * c

opt = Opt(:LD_SLSQP, length(p))
opt.lower_bounds = zeros(length(p))
opt.upper_bounds = ones(length(p))
opt.xtol_rel = 1e-15

myfunc(x, grad) = optfn!(x, grad, p)

opt.min_objective = myfunc
equality_constraint!(opt, eqcon!)

(minf,minx,ret) = optimize(opt, dtc)

using JuMP, Ipopt
m = Model(Ipopt.Optimizer)
#m = Model(solver=NLoptSolver(algorithm=:LD_SLSQP))

@variable(m, x[1:2])
@NLobjective(m, Min, (x[1]-3)^3 + (x[2]-4)^2)
@NLconstraint(m, (x[1]-1)^2 + (x[2]+1)^3 + exp(-x[1]) <= 1)

JuMP.optimize!(m)

c = zeros(size(ocosts)[2:end])

if rtracb == T+1
    c .= M
else
    t = ceil(tracker(st[rtraca]) + 1e-12)
    if rtraca == rtracb
        c .= ocosts[rtraca,..] .+ α * (rtraca - t)
    else
        fracs = t .- tracker.(st[rtraca:rtracb])
        fracs[end] = clamp(fracs[end], 0., Inf)
        fracs[1:(end-1)] .-= fracs[2:end]

        for (f,rt) in zip(fracs, rtraca:rtracb)
            c .+= f * (α * (rt .- t) .+ ocosts[rt,..])
        end
    end
end

t = ceil(tracker(st[rtraca]) + 1e-12)
if rtraca == rtracb
    c .= ocosts[rtraca,..] .+ α * (rtraca - t)
else
    fracs = t .- tracker.(st[rtraca:rtracb])
    fracs[end] = clamp(fracs[end], 0., Inf)
    fracs[1:(end-1)] .-= fracs[2:end]

    for (f,rt) in zip(fracs, rtraca:rtracb)
        c .+= f * (α * (rt .- t) .+ ocosts[rt,..])
    end
end

fracs = t .- tracker.(st[rtraca:rtracb])
fracs[end] = clamp(fracs[end], 0., Inf)
fracs[1:(end-1)] .-= fracs[2:end]

for (f,rt) in zip(fracs, rtraca:rtracb)
    c .+= f * (α * (rt .- t) .+ ocosts[rt,..])
end

revtracs2 = zeros(Int, nlinks, T)
for i in 1:nlinks
    l = length(link(net, i))
    A = Int.(ceil.(tracker.(states[i,:]) .+ 1e-10))

    revtracs2[i,A[end]:end] .= T+1
    revtracs2[i,A[end]-1] = argfilter(x -> x == A[end], A)[1]
    for t in (A[end]-2):-1:1
        revtracs2[i,t] = (t+l-1) + argfilter(x -> x > t, A[(t+l):revtracs2[i,t+1]])[1]
        #revtracs2[i,t] = (t+l-1) + argfilter(x -> x <= t, A[t+l:revtracs2[i,t+1]])[end]
    end
end

    lastt = Int(floor(tracker(states[i,end])))+1
    revtracs2[i,lastt:end] .= T+1

    lastt2 = argfilter(x -> x >= lastt-1, tracker.(states[i,:]))[1]
    revtracs2[i,lastt-1] = lastt2

"""
    lastt2 = argfilter(x -> x >= lastt-1, tracker.(states[i,1:(lastt-1)]))[1]
    revtracs2[i,lastt2:(lastt-1)] .= T
    println(lastt2)
"""
    #revtracs2[i,lastt-1] = T
    for t in (lastt-2):-1:1
        revtracs2[i,t] = (t+l-1) + argfilter(x -> x >= t-1, tracker.(states[i,(t+l):revtracs2[i,t+1]]))[1]
    end
end

c = zeros(size(ocosts)[2:end])

if rtracb == T+1
    c .= M
else
    t = ceil(tracker(st[rtracb]) - 1e-10)
    if rtraca == rtracb
        c .= ocosts[rtraca,..] .+ α * (rtraca - t)
    else
        fracs = tracker.(st[rtraca:rtracb]) .- (t-1)
        fracs[1] = clamp(fracs[1], 0., Inf)
        fracs[2:end] .-= fracs[1:(end-1)]
        #fracs[end] = clamp(fracs[end], 0., Inf)
        #fracs[1:(end-1)] .-= fracs[2:end]

        for (f,rt) in zip(fracs, rtraca:rtracb)
            c .+= f * (α * (rt .- t) .+ ocosts[rt,..])
        end
    end
end

t = floor(tracker(st[rtracb]) + 1e-10)
fracs = tracker.(st[rtraca:rtracb]) .- (t-1)
fracs[1] = clamp(fracs[1], 0., Inf)
fracs[2:end] .-= fracs[1:(end-1)]
#fracs[end] = clamp(fracs[end], 0., Inf)
#fracs[1:(end-1)] .-= fracs[2:end]

for (f,rt) in zip(fracs, rtraca:rtracb)
    c .+= f * (α * (rt .- t) .+ ocosts[rt,..])
end
