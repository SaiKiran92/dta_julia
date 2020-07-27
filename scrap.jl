

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


t, div = 157, 3

ili = inlinkids(net,div)[1]
oli = outlinkids(net,div)
s = svalue(states[ili,t])
r = rvalue.(states[oli,t])

# calculating flows
sr = [SR[div][ti][t,..] for ti in 1:2]
trkr = tracker(states[ili,t])
l = length(link(net, ili))
fₐ, va = dflows(((t > l) ? Fᵢ[ili, trkr:(t - l),..] : zeros(1, size(Fᵢ)[3:end]...)), sr, r, s, fracmoved(states[ili,t]))

fᵢ = ((t > l) ? Fᵢ[ili, trkr:(t - l),..] : zeros(1, size(Fᵢ)[3:end]...))
m = fracmoved(states[ili,t])

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

x, y = a[τ], a[τ+1]
if y == T+1
    ECᵢ[i,τ,..] .= M
elseif x == y
    ECᵢ[i,τ,..] .= (x - τ) .+ ECₒ[i,x,..]
else
    ECᵢ[i,τ,..] .= ((x - τ) .+ ECₒ[i,x,..]) * safedivide(max(0., CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)), TFᵢ[i,τ])
    ECᵢ[i,τ,..] .+= ((y - τ) .+ ECₒ[i,y,..]) * safedivide(max(0., CFₒ[i,τ] - ((τ > 1) ? CFᵢ[i,y-1] : 0.)), TFᵢ[i,τ])
    if y > x+1
        ECᵢ[i,τ,..] .+= squeezesum(((((x+1):(y-1)) .- τ) .+ ECₒ[i,(x+1):(y-1),..]) .* expand(safedivide.(TFₒ[i,(x+1):(y-1)], Ref(TFᵢ[i,τ])), dims=(2,3)), dims=1)
    end
end


function costupdate!(i, τ)

end

x, y = a[i,τ], a[i,τ+1]
if y == T+1
    ECᵢ[i,τ,..] .= M
elseif x == y
    ECᵢ[i,τ,..] .= (x - τ) .+ ECₒ[i,x,..]
else
    ECᵢ[i,τ,..] .= ((x - τ) .+ ECₒ[i,x,..]) * safedivide(max(0., CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)), TFᵢ[i,τ])
    ECᵢ[i,τ,..] .+= ((y - τ) .+ ECₒ[i,y,..]) * safedivide(max(0., CFₒ[i,τ] - ((τ > 1) ? CFᵢ[i,y-1] : 0.)), TFᵢ[i,τ])
    if y > x+1
        ECᵢ[i,τ,..] .+= squeezesum(((((x+1):(y-1)) .- τ) .+ ECₒ[i,(x+1):(y-1),..]) .* expand(safedivide.(TFₒ[i,(x+1):(y-1)], Ref(TFᵢ[i,τ])), dims=(2,3)), dims=1)
    end
end

ECᵢ[i,τ,..] .= ((x - τ) .+ ECₒ[i,x,..]) * safedivide(max(0., CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)), TFᵢ[i,τ])
ECᵢ[i,τ,..] .+= ((y - τ) .+ ECₒ[i,y,..]) * safedivide(max(0., CFₒ[i,τ] - CFᵢ[i,y-1]), TFᵢ[i,τ], 0.)
if y > x+1
    ECᵢ[i,τ,..] .+= squeezesum(((((x+1):(y-1)) .- τ) .+ ECₒ[i,(x+1):(y-1),..]) .* expand(safedivide.(TFₒ[i,(x+1):(y-1)], Ref(TFᵢ[i,τ])), dims=(2,3)), dims=1)
end

ECᵢ[i,τ,..] .= round.(ECᵢ[i,τ,..], digits=ROUND_DIGITS)

d, t = 13, 137

ili = inlinkids(net,d)[1]
oli = outlinkids(net,d)
s = svalue(states[ili,t])
r = rvalue.(states[oli,t])

# calculating flows
sr = [SR[div][ti][t,..] for ti in 1:2]
trkr = tracker(states[ili,t])
l = length(link(net, ili))
fₐ, va = dflows(((t > l) ? Fᵢ[ili, trkr:(t - l),..] : zeros(1, size(Fᵢ)[3:end]...)), sr, r, s, fracmoved(states[ili,t]))

x, y = a[i,τ], a[i,τ+1]
if y == T+1
    ECᵢ[i,τ,..] .= M
elseif x == y
    ECᵢ[i,τ,..] .= (x - τ) .+ ECₒ[i,x,..]
else
    ECᵢ[i,τ,..] .= ((x - τ) .+ ECₒ[i,x,..]) * safedivide(max(0., CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)), TFᵢ[i,τ])
    ECᵢ[i,τ,..] .+= ((y - τ) .+ ECₒ[i,y,..]) * safedivide(max(0., CFₒ[i,τ] - CFᵢ[i,y-1]), TFᵢ[i,τ], 0.)
    if y > x+1
        ECᵢ[i,τ,..] .+= squeezesum(((((x+1):(y-1)) .- τ) .+ ECₒ[i,(x+1):(y-1),..]) .* expand(safedivide.(TFₒ[i,(x+1):(y-1)], Ref(TFᵢ[i,τ])), dims=(2,3)), dims=1)
    end
end

ECᵢ[i,τ,..] .= ((x - τ) .+ ECₒ[i,x,..]) * safedivide(max(0., CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)), TFᵢ[i,τ])
ECᵢ[i,τ,..] .+= ((y - τ) .+ ECₒ[i,y,..]) * safedivide(max(0., CFₒ[i,τ] - CFᵢ[i,y-1]), TFᵢ[i,τ], 0.)
if y > x+1
    ECᵢ[i,τ,..] .+= squeezesum(((((x+1):(y-1)) .- τ) .+ ECₒ[i,(x+1):(y-1),..]) .* expand(safedivide.(TFₒ[i,(x+1):(y-1)], Ref(TFᵢ[i,τ])), dims=(2,3)), dims=1)
end

t, div = 7, 3

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

CFₒ[ili,t] = sum(Fₒ[ili,t,..]) + ((t > 1) ? CFₒ[ili,t-1] : 0.)
A[ili], fracs[ili] = cumvalarg((@view CFᵢ[ili,:]), CFₒ[ili,t], tracker(states[ili,t]), t - l)

tracker(states[ili,t])
t-l


cumr, v, starti, stopi = (@view CFᵢ[ili,:]), CFₒ[ili,t], tracker(states[ili,t]), t - l

if starti <= 0
    return (starti+1, 0.)
end
i = starti
while i <= stopi
    if v < cumr[i]
        tmp1 = (i > 1) ? cumr[i-1] : 0.
        tmp2 = (i > 1) ? max(0., cumr[i] - cumr[i-1]) : 0.
        return (i, safedivide((v - tmp1), tmp2))
    end
    i += 1
end

if v < cumr[i]
    tmp1 = (i > 1) ? cumr[i-1] : 0.
    tmp2 = (i > 1) ? max(0., cumr[i] - cumr[i-1]) : cumr[i]
    return (i, safedivide((v - tmp1), tmp2))
end

tmp1 = (i > 1) ? cumr[i-1] : 0.
tmp2 = (i > 1) ? max(0., cumr[i] - cumr[i-1]) : 0.

i, τ = 3, 389
ili = inlinkids(net,div)[1]
oli = outlinkids(net,div)

for li in oli
    costupdate!(li,τ)
end

i = 1
x, y = a[i,τ], a[i,τ+1]
if y == T+1
    ECᵢ[i,τ,..] .= M
elseif x == y
    ECᵢ[i,τ,..] .= (x - τ) .+ ECₒ[i,x,..]
else
    tf = max(0., CFₒ[i,τ] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.))
    ECᵢ[i,τ,..] .= ((x - τ) .+ ECₒ[i,x,..]) * safedivide(max(0., CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)), tf)
    ECᵢ[i,τ,..] .+= ((y - τ) .+ ECₒ[i,y,..]) * safedivide(max(0., CFₒ[i,τ] - CFᵢ[i,y-1]), tf, 0.)
    if y > x+1
        ECᵢ[i,τ,..] .+= squeezesum(((((x+1):(y-1)) .- τ) .+ ECₒ[i,(x+1):(y-1),..]) .* expand(safedivide.(TFₒ[i,(x+1):(y-1)], Ref(TFᵢ[i,τ])), dims=(2,3)), dims=1)
    end
end


tf = max(0., CFₒ[i,τ] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.))
((x - τ) .+ ECₒ[i,x,..]) * safedivide(max(0., CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)), tf)
((y - τ) .+ ECₒ[i,y,..]) * safedivide(max(0., CFₒ[i,τ] - CFᵢ[i,y-1]), tf, 0.)

τ, d = 140, 3

ili = inlinkids(net,d)[1]
oli = outlinkids(net,d)

i = ili
x, y = a[i,τ], a[i,τ+1]
if y == T+1
    ECᵢ[i,τ,..] .= M
elseif x == y
    ECᵢ[i,τ,..] .= (x - τ) .+ ECₒ[i,x,..]
else
    tf = max(0., CFₒ[i,τ] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.))
    ECᵢ[i,τ,..] .= ((x - τ) .+ ECₒ[i,x,..]) * safedivide(max(0., CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)), tf)
    ECᵢ[i,τ,..] .+= ((y - τ) .+ ECₒ[i,y,..]) * safedivide(max(0., CFₒ[i,τ] - CFᵢ[i,y-1]), tf, 0.)
    #ECᵢ[i,τ,..] .= ((x - τ) .+ ECₒ[i,x,..]) * safedivide(max(0., CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)), TFᵢ[i,τ])
    #ECᵢ[i,τ,..] .+= ((y - τ) .+ ECₒ[i,y,..]) * safedivide(max(0., CFₒ[i,τ] - CFᵢ[i,y-1]), TFᵢ[i,τ], 0.)
    if y > x+1
        ECᵢ[i,τ,..] .+= squeezesum(((((x+1):(y-1)) .- τ) .+ ECₒ[i,(x+1):(y-1),..]) .* expand(safedivide.(TFₒ[i,(x+1):(y-1)], Ref(TFᵢ[i,τ])), dims=(2,3)), dims=1)
    end
end

tf = max(0., CFᵢ[i,τ] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.))
((x - τ) .+ ECₒ[i,x,..]) * safedivide(max(0., CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)), tf)
((y - τ) .+ ECₒ[i,y,..]) * safedivide(max(0., CFₒ[i,τ] - CFᵢ[i,y-1]), tf, 0.)

safedivide(max(0., CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)), tf)
CFₒ[i,x] - ((τ > 1) ? CFᵢ[i,τ-1] : 0.)
