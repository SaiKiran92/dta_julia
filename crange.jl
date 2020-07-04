import Base: getindex, first, last, size, show

struct ContinuousRange{T<:AbstractFloat} <: AbstractRange{T}
    a::T
    b::T

    function ContinuousRange(a::T, b::T) where {T<:AbstractFloat}
        a = clamp(a, 0., Inf)
        b = clamp(b, 0., Inf)
        new{T}(a, b)
    end
end

(::Colon)(a::T, b::T) where {T<:AbstractFloat} = ContinuousRange(a,b)

ϵ = 1e-10;
start(c::ContinuousRange) = c.a
stop(c::ContinuousRange) = c.b
firstidx(c::ContinuousRange) = Int(floor(start(c)+ϵ)+1)
lastidx(c::ContinuousRange) = Int(ceil(stop(c)))
size(c::ContinuousRange) = c.b - c.a
show(io::IO, c::ContinuousRange) = print(io, repr(first(c)), ':', repr(last(c)))

function getindex(v::Array{T}, c::ContinuousRange, rest...) where {T<:AbstractFloat}
    ai,bi = firstidx(c), lastidx(c)
    ai = (ai <= 0) ? 1 : ai
    bi = (bi == 0) ? -1 : bi

    if ai > bi
        if typeof(rest[1]) == EllipsisNotation.Ellipsis
            return zeros(1, size(v)[2:end]...)
        else
            u = [i[1] for i in size.(rest) if length(i)!=0]
            return zeros(1, u...)
        end
    end

    rv = getindex(v, ai:bi, rest...)
    if ai == bi
        rv[1,..] .*= size(c)
    else
        rv[1,..] .*= (firstidx(c) - start(c))
        rv[end,..] .*= (stop(c) - lastidx(c) + 1)
    end
    return rv
end

getindex(a::Array, c::ContinuousRange, rest...) = getindex(convert(Array{Float64}, a), c, rest...)
getindex(a::Array, i::Int, c::ContinuousRange, rest...) = getindex(a[i,..], c, rest...)
