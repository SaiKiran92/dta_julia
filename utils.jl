
ϵ = 1e-13
safedivide(a, b) = (a+ϵ)/(b+ϵ)
positive(x::Real) = (x > 0.)
positive(x::Array) = all(pos.(x))
approxpositive(x::Real) = (x >= 0.) || (x ≈ 0.) #isapprox(x, 0., atol=ATOL)
approxpositive(x::Array) = all(approxpos.(x))
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
    dims = sort([dims...]) .+ (0:length(dims)-1)
    i = [size(a)...]
    for d in dims
        insert!(i, d, 1)
    end
    reshape(a, i...)
end

squeezesum(a::AbstractArray; dims) = squeeze(sum(a, dims=dims), dims=dims)

function dijkstra(net::Network{T}, snk::U) where {T<:Integer, U<:Integer}
    nnodes = numnodes(net)

    costs = fill(typemax(U), nnodes)
    childvec = zeros(T, nnodes)

    P = PriorityQueue{T,U}()
    P[snk] = costs[snk] = 0
    childvec[snk] = snk

    while !isempty(P)
        v = dequeue!(P)
        for u in inneighbors(net, v)
            lidx = idx(net, u, v)
            alt = costs[v] + net.links[lidx].l
            if (costs[u] > alt)
                P[u] = costs[u] = alt
                childvec[u] = v
            end
        end
    end
    childvec
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

argcumval(a::Array{<:AbstractFloat,0}, val) = 0.
argcumval(a::Array, val) = argcumval(squeezesum(a, dims=(2:ndims(a))), val)

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
