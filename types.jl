import Base: length

struct Link{T<:Integer, U<:Integer}
    i::T
    j::T
    l::U
end

length(link::Link) = link.l
upn(link::Link) = link.i
dwn(link::Link) = link.j

@enum NodeType source sink merge diverge
struct Node{T<:Integer}
    id::T
    type::NodeType
end

struct Network{T<:Integer}
    # representations
    ## forward star
    fadjlist::Vector{Vector{T}}
    fpoints::Vector{T}

    ## backward star
    badjlist::Vector{Vector{T}}
    bpoints::Vector{T}

    # data
    ## links
    links::Vector{Link}

    ## nodes
    nodes::Vector{Union{Missing,Node}}
end

numnodes(net::Network) = length(net.nodes)
numlinks(net::Network) = length(net.links)

sources(net::Network) = filter(i -> length(net.badjlist[i]) == 0, 1:numnodes(net))
sinks(net::Network) = filter(i -> length(net.fadjlist[i]) == 0, 1:numnodes(net))
merges(net::Network) = filter(i -> length(net.badjlist[i]) > 1, 1:numnodes(net))
diverges(net::Network) = filter(i -> length(net.fadjlist[i]) > 1, 1:numnodes(net))

haslink(net::Network{T}, i::U, j::U) where {T<:Integer, U<:Integer} = (j in net.fadjlist[i])

idx(net::Network{T}, i::U, j::U) where {T<:Integer, U<:Integer} = haslink(net, i, j) ? net.fpoints[i]-1+searchsortedfirst(net.fadjlist[i],j) : -1
idx(net::Network{T}, l::Link) where {T<:Integer} = idx(net, l.i, l.j)
outneighbors(net::Network, i) = net.fadjlist[i]
inneighbors(net::Network, j) = net.badjlist[j]

outlinkids(net::Network, i) = collect(net.fpoints[i]:(net.fpoints[i+1]-1))
inlinkids(net::Network, j) = idx.(Ref(net), net.badjlist[j], Ref(j))

outlinks(net::Network, i) = net.links[outlinkids(net, i)]
inlinks(net::Network, j) = net.links[inlinkids(net, j)]

link(net::Network, i) = net.links[i]

function Network(linkdata::DataFrame)
    nnodes = max(maximum(linkdata.i), maximum(linkdata.j))
    net = Network(nnodes)

    for row in eachrow(linkdata)
        add_link!(net, Link(row.i, row.j, row.l))
    end

    for i in 1:nnodes
        if length(net.fpoints[i]) == 0
            net.nodes[i] = Node(i, sink)
        elseif length(net.bpoints[i]) == 0
            net.nodes[i] = Node(i, source)
        elseif length(net.fpoints[i]) > 1
            net.nodes[i] = Node(i, diverge)
        else
            net.nodes[i] = Node(i, merge)
        end
    end

    net
end

function Network(nnodes::T) where {T<:Integer}
    fadjlist = [Vector{T}() for _ in one(T):nnodes]
    fpoints = ones(T, nnodes+1)

    badjlist = [Vector{T}() for _ in one(T):nnodes]
    bpoints = ones(T, nnodes+1)

    links = Vector{Link}()

    nodes = fill(missing, nnodes)

    Network{T}(fadjlist, fpoints, badjlist, bpoints, links, nodes)
end

function add_link!(net, l::Link)
    i,j = l.i, l.j

    @inbounds list = net.badjlist[j]
    index = searchsortedfirst(list, i)
    insert!(list, index, i)
    net.bpoints[(j+1):end] .+= 1

    @inbounds list = net.fadjlist[i]
    index = searchsortedfirst(list, j)
    insert!(list, index, j)
    net.fpoints[(i+1):end] .+= 1

    # add to list
    insert!(net.links, net.fpoints[i]-1+index, l)
end

@with_kw struct LinkState
    trac::AbstractFloat # tracker - floating point value
    #frac::AbstractFloat # fraction of the flow already moved from the cohort of 'trac'
    sval::AbstractFloat # sending value
    rval::AbstractFloat # receiving value
end

tracker(ls::LinkState) = ls.trac
#fracmoved(ls::LinkState) = ls.frac
svalue(ls::LinkState) = ls.sval
rvalue(ls::LinkState) = ls.rval
