
"""
    dflows(fᵢ, sr, r, s, m)

Compute flows at a diverging intersection.

# Arguments
- `fᵢ::Array`: slice from inflows from the range relevant to the computation
- `sr::Array`: splitting rates
- `r::Array`: receiving capacities
- `s::Float`: sending capacity
- `m::Float`: flow already moved from the first timestep of fᵢ

"""
function dflows(fᵢ, sr, r, s, m)
    #sr = Dict(k => expand(v, dims=1) for (k,v) in pairs(sr))
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
        push!(fₐ, argcut(fᵢ, va, m) .* v) #[1,..])
    end

    return (fₐ, va)
end

"""
    mflows(fᵢ, sr, r, s, m)

Compute flows at a merging intersection.

# Arguments
- `fᵢ::Array`: slice from inflows from the range relevant to the computation
- `r::Float`: receiving capacity
- `s::Array`: sending capacities
- `m::Float`: flow already moved from the first timestep of fᵢ

"""
function mflows(fᵢ, r, s, m)
    fₐ = [zeros(size(fᵢ[1])[2:end]...) for i in 1:2]
    va = [0., 0.]

    function update!(i, cap)
        va[i] = valarg(fᵢ[i], cap, m[i]*sum(fᵢ[i][1,..]))
        fₐ[i] .= argcut(fᵢ[i], va[i], m[i])
    end

    if r >= sum(s)
        update!(1, s[1]) #free
        update!(2, s[2]) #free
    elseif s[1] < 0.5*r
        update!(1, s[1]) #free
        update!(2, r-s[1]) #cong
    elseif s[2] < 0.5*r
        update!(1, r-s[2]) #free
        update!(2, s[2]) #cong
    else
        update!(1, 0.5*r) #cong
        update!(2, 0.5*r) #cong
    end

    return (fₐ, va)
end
