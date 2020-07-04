
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

function argcumval(vec::Vector, val, from=0.)
    vec = vec[from:end]

    if vec[1] > val
        return round(from + (val/vec[1]) * (1 - decimal(from)), digits=30)
    end

    tvec = cumsum(vec) .- val
    try
        i = argfilter(x -> (x > 0.), tvec)[1]
        return round(i - tvec[i]/vec[i], digits=30)
    catch BoundsError
        return length(vec)
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

"""
function solveqp(c, dtc, λ)
    p = dtc .- λ * c

    n = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    @variable(n, x[1:length(p)])
    o = @objective(n, Min, 0.5*sum((x .- p).^2))
    @constraint(n, sum(x) == 1.)
    @constraint(n, x .>= 0.)
    @constraint(n, x .<= 1.)

    soln = JuMP.optimize!(n)

    return round.(JuMP.value.(x), digits=7)
end
"""

function proportionalize(x; digits=6)
    x[x .< 0.] .= 0
    x = round.(x, digits=digits)
    x ./= sum(x)
    return x
end
