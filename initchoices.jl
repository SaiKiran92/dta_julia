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
    DTC = Dict(k => zeros(T) for k in keys(trips));
    for k in keys(trips)
        DTC[k][1:Tm] .= rand(Tm)
        DTC[k][1:Tm] ./= sum(DTC[k][1:Tm])
        #DTC[k][1:Tm] .= 1/Tm
        #DTC[k][1] = 1.
    end

    pathvecs = Dict(snk => dijkstra(net, snk) for snk in snks);
    SR = Dict(div => Dict(i => zeros(T, nsinks, nclasses) for i in 1:2) for div in divs);
    for div in divs
        for (snkno,snk) in enumerate(snks)
            for cls in 1:nclasses
                if (outneighbors(net, div)[1] == pathvecs[snk][div])
                    SR[div][1][:,snkno,cls] .= 1.
                elseif (outneighbors(net, div)[2] == pathvecs[snk][div])
                    SR[div][2][:,snkno,cls] .= 1.
                else
                    error()
                end
            end
        end
    end
    return (DTC, SR)
end
