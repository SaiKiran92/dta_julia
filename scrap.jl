i, τ = 1, 138

if rtracs[i,τ+1] == T
    ECᵢ[i,τ,..] .= M
else
    #@show i, τ
    rng = rtracs[i,τ]:rtracs[i,τ+1]
    inds = indexify(rng)
    #@show i, τ, rng, inds
    if size(rng) == 0.
        ind = lastidx(rng)
        ECᵢ[i,τ,..] .= (ind - τ) .+ ECₒ[i,ind,..]
    else
        tmp = mucsum(@view(CFₒ[i,:]), rng)
        tmp[end] += ϵ;
        tf = (CFᵢ[i,τ] - ((τ == 1) ? 0. : CFᵢ[i,τ-1])) + ϵ
        #@show rng, inds
        #@show tmp, tf, (tmp ./ tf), (inds .- τ)
        ECᵢ[i,τ,..] .= squeezesum((tmp ./ tf) .* ((inds .- τ) .+ ECₒ[i,inds,..]), dims=1)
        #@show i, τ, ECᵢ[i,τ,:,1]
    end
end

v = @view(CFₒ[i,:])

i, j = firstidx(rng), lastidx(rng)
rv = mucsum(v[i:j] .- ((i == 1) ? 0. : v[i-1]))

if i == j
    println("here")
    rv .*= size(rng)
elseif i < j
    println("here2")
    rv[1] *= (1 - decimal(start(rng)))
    rv[end] *= decimal(stop(rng))
end

mucsum(@view(CFₒ[i,:]), rng)