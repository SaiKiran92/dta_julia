
ϵ = 1e-13
safedivide(a, b, v=1.) = (a+(ϵ*v))/(b+ϵ)
decimal(a) = a - floor(a)
argfilter(fn, x) = [i for (i,x) in enumerate(x) if fn(x)]

function getindex(arr::AbstractArray, f::AbstractFloat, rest...)
    if f < 0.
        f = 0.
    end
    i = Int(ceil(f + ϵ))
    a = (i > 1) ? arr[i-1, rest...] : zeros(size(arr)[2:end])
    b = arr[i, rest...]
    return a .+ (b .- a) .* decimal(f)
end

function squeeze(a; dims=nothing)
    if dims === nothing
        return a[((i == 1) ? 1 : Colon() for i in size(a))...]
    else
        return a[((i ∈ dims) ? 1 : Colon() for i in 1:ndims(a))...]
    end
end

function expand(a; dims)
    dims = sort([dims...])
    i = [size(a)...]
    for d in dims
        insert!(i, d, 1)
    end
    reshape(a, i...)
end

squeezesum(a::AbstractArray; dims) = squeeze(sum(a, dims=dims), dims=dims)

function mucsum(a::AbstractArray; dims::Integer=1)
    b = copy(a)
    colonvec(i) = (tmp = Array{Any}(fill(:, ndims(a))); tmp[dims] = i; tmp)
    for i in size(a)[dims]:-1:2
        b[colonvec(i)...] -= b[colonvec(i-1)...]
    end
    b
end

function mucsum(v::AbstractVector, rng::ContinuousRange)
    i, j = firstidx(rng), lastidx(rng)
    rv = mucsum(v[i:j] .- ((i == 1) ? 0. : v[i-1]))

    if i == j
        rv .*= size(rng)
    elseif i < j
        rv[1] *= (1 - decimal(start(rng)))
        rv[end] *= decimal(stop(rng))
    end
    return rv
end

function cvalarg(cvec::AbstractVector, cval, starttrac, maxt)
    starttrac = (starttrac < 0) ? 0. : starttrac
    t = Int(ceil(starttrac + ϵ))
    while t <= maxt
        if cvec[t] > cval
            tmp = (t > 1) ? cvec[t-1] : 0.
            return (t-1) + (cval - tmp)/(cvec[t] - tmp)
        else
            t += 1
        end
    end
    return maxt
end

function rcvalarg(cvec::AbstractVector, cval, endtrac, mint)
    t = Int(floor(endtrac))
    while (t >= mint)
        if (cvec[t] < cval)# && !(cvec[t] ≈ cval)
            return (t+1) - (cvec[t+1] - cval)/(cvec[t+1] - cvec[t])
        else
            t -= 1
        end
    end
    return mint
end

function valarg(r::AbstractArray, v, starti=1, maxi=size(r)[1])
    starti = (starti < 1) ? 1 : starti
    i = starti
    while i <= maxi
        totr = sum(r[i,..])
        if v < totr
            return (i-1) + (v/totr)
        else
            v -= totr
        end
        i += 1
    end

    return maxi
end

argcut(r::AbstractArray, j, i) = squeezesum(r[i:j,..], dims=1)

function projpx(u::Array) # projection onto a probability simplex
    v = sort(u, rev=true)
    
    vsum, lv = 0., 0.
    for j in 1:length(u)
        if v[j] + (1.0/j)*(1 - vsum - v[j]) > 0.
            vsum += v[j]
        else
            lv = (1.0/(j-1))*(1 - vsum) # lagrange variable
            break
        end
    end

    return max.(u .+ lv, Ref(0.))
end