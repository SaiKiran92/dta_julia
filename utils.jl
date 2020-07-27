
ϵ = 1e-13
safedivide(a, b, v=1.) = (a+(ϵ*v))/(b+ϵ)
decimal(a) = a - floor(a)
Base.getindex(d::Dict{Tuple{T,T}}, i, j) where {T<:Integer} = d[(i,j)]
argfilter(fn, x) = [i for (i,x) in enumerate(x) if fn(x)]

function squeeze(a; dims=nothing)
    if dims == nothing
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

function valarg(r::Array, a, b)
    v = a + b
    i, l = 0, size(r)[1]
    while i < l
        totr = sum(r[i+1,..])
        if v < totr
            return i + (v/totr)
        else
            v -= totr
        end
        i += 1
    end

    return l
end

argcut(r::Array, j, i) = squeezesum(r[i:j,..], dims=1)

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
