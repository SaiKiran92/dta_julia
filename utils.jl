
ϵ = 1e-13
safedivide(a, b) = (a+ϵ)/(b+ϵ)
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
    dims = sort([dims...])# .+ (0:length(dims)-1)
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

"""
function argcut(r::Array, j, i)
    c = (r isa Vector) ? [0.] : zero(r[1,..])
    inti, intj = Int(floor(i))+1, Int(ceil(j))
    if (intj == j)
        intj -= 1
        dec = 1.
    else
        dec = decimal(j)
    end

    @show (i, " ", j)
    @show (inti, " ", intj)

    if inti == intj
        c .= r[inti,..] .* (j - i)
    else
        c .= r[inti,..] .* (inti - i)
        c .+= r[intj,..] .* dec
        if intj > inti + 1
            c .+= squeezesum(r[(inti+1):(intj-1),..], dims=1)
        end
    end

    return c
end
"""
