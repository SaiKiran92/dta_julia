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

function initchoices()
    # choice initialization
    DTC = zeros(T, size(trips)...);
    #DTC[1,..] .= 1.
    DTC[1:Tm,..] .= rand(size(DTC[1:Tm,..])...)
    DTC ./= sum(DTC, dims=1)

    pathvecs = Dict(k => dijkstra(net, k) for k in snks);
    SR = Dict(d => Dict(i => zeros(T, nsinks, nclasses) for i in 1:2) for d in divs);
    for d in divs
        for (kno,k) in enumerate(snks)
            for c in 1:nclasses
                if (outneighbors(net, d)[1] == pathvecs[k][d])
                    SR[d][1][:,kno,c] .= 1.
                elseif (outneighbors(net, d)[2] == pathvecs[k][d])
                    SR[d][2][:,kno,c] .= 1.
                else
                    error()
                end
            end
        end
    end
    return (DTC, SR)
end
