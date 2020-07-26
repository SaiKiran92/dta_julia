

# reverse states - for computing costs - only trackers needed
rstates = zeros(nlinks, T+1)
tinflows = round.(squeezesum(inflows, dims=(3,4)), digits=ROUND_DIGITS)
toutflows = round.(squeezesum(outflows, dims=(3,4)), digits=ROUND_DIGITS)
cinflows = round.(cumsum(tinflows, dims=2), digits=ROUND_DIGITS)
coutflows = round.(cumsum(toutflows, dims=2), digits=ROUND_DIGITS)
for i in 1:nlinks
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
    rstates[i, (lastinidx+1):(maxt-1)] .= collect((lastinidx+1):(maxt-1)) .+ (l-1)

    for t in 2:lastinidx
        ui = floor(rstates[i,t-1] + 1e-9)
        rstates[i,t] = ui + round(argcumval2(coutflows[i,ui:lastoutidx], cinflows[i,t-1], 0., :zero_exclude), digits=ROUND_DIGITS)
    end
end


incosts = zeros(nlinks, T, nsinks, nclasses)
outcosts = zeros(nlinks, T, nsinks, nclasses)
incosts[:,end,..] .= M
outcosts[:,end,..] .= M
for (snkid,snk) in enumerate(snks)
    i = inlinkids(net, snk)[1]
    outcosts[i,:,:,1] .= M
    r = collect(1:T)
    outcosts[i,:,snkid,1] .= clamp.(trgt .- r, 0., Inf) * β .+ clamp.(r .- trgt, 0., Inf) * γ # overwriting
end

function costupdate!(i,t)
    println(i, " ", t)
    if rstates[i,t+1] == T
        incosts[i,t,..] .= M
    else
        r = rstates[i,t]:rstates[i,t+1]
        println(r)
        println(size(outcosts[i,r,..]))
        println(size(collect(firstidx(r):lastidx(r))))
        incosts[i,t,..] .= sum((outcosts[i,r,..] .+ collect(firstidx(r):lastidx(r)) .- t) .* safedivide.(toutflows[i,r], Ref(tinflows[i,t])), dims=1)[1,..]
    end
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

        outcosts[ili[1],t,..] .= incosts[oli,t,..]
        outcosts[ili[2],t,..] .= incosts[oli,t,..]
    end

    for div in divs
        ili = inlinkids(net,div)[1]
        oli = outlinkids(net,div)

        for i in ili
            costupdate!(i,t)
        end

        sr = [srates[div][i][t,..] for i in 1:2]
        outcosts[ili,t,..] .= sr[1] .* incosts[oli[1],t,..] .+ sr[2] .* incosts[oli[2],t,..]
    end
end


i, t = 1,2
r = rstates[i,t]:rstates[i,t+1]
tmpr = collect(firstidx(r):lastidx(r))
if size(tmpr) == (0,) # when rstates[i,t] == rstates[i,t+1] == Int(rstates[i,t])
    tmpr = firstidx(r)
    incosts[i,t,..] .= outcosts[i,tmpr,..] .+ tmpr .- t
else
    incosts[i,t,..] .= sum((outcosts[i,r,..] .+ tmpr .- t) .* safedivide.(toutflows[i,r], Ref(tinflows[i,t])), dims=1)[1,..]
end


i = 1
t = 2
ui = floor(rstates[i,t-1] + 1e-9)
cinflows[i,t-1]
#rstates[i,t] = ui + round(argcumval2(coutflows[i,ui:lastoutidx], cinflows[i,t-1], 0., :zero_exclude), digits=ROUND_DIGITS)
ui + round(argcumval2(coutflows[i,ui:lastoutidx], cinflows[i,t-1], 0., :zero_exclude), digits=ROUND_DIGITS)

t = 3
ui = floor(rstates[i,t-1] + 1e-9)
cinflows[i,t-1]
coutflows[i,ui:lastoutidx]
#rstates[i,t] = ui + round(argcumval2(coutflows[i,ui:lastoutidx], cinflows[i,t-1], 0., :zero_exclude), digits=ROUND_DIGITS)
ui + round(argcumval2(coutflows[i,ui:lastoutidx], cinflows[i,t-1], 0., :zero_exclude), digits=ROUND_DIGITS)


@time for i in 1:100
    foo(prod(1:10000), 1)
end

@time foo.(Ref(prod(1:10000)), 1)

function foo(a, b)

end



function argcumval(vec::Vector, val, from=0., z=:zeroinclude)
    if z == :zeroinclude
        f = approxpositive
    else
        f = positive
    end
    cutval = from * vec[1]
    if vec[1] - cutval >= val
        return from + safedivide(val, vec[1])
    else
        val -= (vec[1] - cutval)
        i, l = 2, length(vec)
        while (i <= l) && f(val - vec[i])
            val -= vec[i]
            i += 1
        end

        if (i > l) && f(val)
            return l
        else
            return i-1 + safedivide(val, vec[i])
        end
    end
end

function argcumval2(vec::Vector, val, from=0., z=:zeroinclude)
    if z == :zeroinclude
        f = approxpositive
    else
        f = positive
    end
    cutval = from * vec[1]
    if vec[1] - cutval >= val
        return from + safedivide(val, vec[1])
    else
        val -= (vec[1] - cutval)
        i, l = 2, length(vec)
        while (i <= l) && f(val - (vec[i] - vec[i-1]))
            val -= (vec[i] - vec[i-1])
            i += 1
        end

        if (i > l) && f(val)
            return l
        else
            return i-1 + safedivide(val, (vec[i] - vec[i-1]))
        end
    end
end

argcumval(a::Array{<:AbstractFloat,0}, val) = 0.
argcumval(a::Array, val) = argcumval(squeezesum(a, dims=(2:ndims(a))), val)


positive(x::Real) = (x > 0.)
positive(x::Array) = all(pos.(x))
approxpositive(x::Real) = (x >= 0.) || (x ≈ 0.)
approxpositive(x::Array) = all(approxpositive.(x))

t, div = 8, 3

ili = inlinkids(net,div)[1]
oli = outlinkids(net,div)
s = min(Q, svalue(states[ili,t]))
r = [min(Q, rvalue(states[li,t])) for li in oli]

# calculating flows
sr = [SR[div][ti][t,..] for ti in 1:2]
trkr = tracker(states[ili,t])
l = length(link(net, ili))
#println(size(Fᵢ)[3:end])
fₐ, va = dflows(((t > l) ? Fᵢ[ili, trkr:(t - l),..] : zeros(1, size(Fᵢ)[3:end]...)), sr, r, s, fracmoved(states[ili,t]))

fᵢ = ((t > l) ? Fᵢ[ili, trkr:(t - l),..] : zeros(1, size(Fᵢ)[3:end]...))
m = fracmoved(states[ili,t])

#sr = Dict(k => expand(v, dims=1) for (k,v) in pairs(sr))
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

3288.0 - 5088.0

1644.0 - 3444

t, mrg = 175, 7
ili = inlinkids(net,mrg)
oli = outlinkids(net,mrg)[1]
s = [min(Q, svalue(states[li,t])) for li in ili]
r = min(Q, rvalue(states[oli,t]))

# calculating flows
trkr = tracker.(states[ili,t])
ls = length.(link.(Ref(net), ili))
f = [(t > l) ? Fᵢ[li, trkr[i]:(t - l),..] : zeros(1, size(Fᵢ)[3:end]...)  for (i,(li,l)) in enumerate(zip(ili, ls))]
fₐ, va = mflows(f, r, s, fracmoved.(states[ili,t]))

fᵢ = f
m = fracmoved.(states[ili,t])

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

r = fᵢ[1]
j, i = va[i], m[i]

c = zero(r[1,..])
inti, intj = Int(floor(i))+1, Int(ceil(j))
if (intj == j)
    intj -= 1
    dec = 1.
else
    dec = decimal(j)
end

if inti == intj
    c .= r[inti,..] .* (j - i)
else
    c .= r[inti,..] .* (inti - i)
    c .+= r[intj,..] .* dec
    if intj > inti + 1
        c .+= squeezesum(r[(inti+1):(intj-1),..], dims=1)
    end
end

c .= r[inti,..] .* (inti - i)
c .+= r[intj,..] .* decimal(j)
if intj > inti + 1
    c .+= squeezesum(r[(inti+1):(intj-1),..], dims=1)
end
