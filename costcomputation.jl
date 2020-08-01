ECᵢ, ECₒ = zero(pᵢ), zero(pᵢ)
ECᵢ[:,end,..] .= M
ECₒ[:,end,..] .= M
for (snkid,snk) in enumerate(snks)
    i = inlinkids(net, snk)[1]
    ECₒ[i,:,:,1] .= M
    r = collect(1:T)
    ECₒ[i,:,snkid,1] .= clamp.(trgt .- r, 0., Inf) * β .+ clamp.(r .- trgt, 0., Inf) * γ # overwriting
end

function costupdate!(i, τ)
    
end

for τ in T:-1:1
    for src in srcs
        i = outlinkids(net, src)[1]
        costupdate!(i,τ)
    end

    for mrg in mrgs
        ili = inlinkids(net,mrg)
        oli = outlinkids(net,mrg)[1]

        costupdate!(oli,τ)

        ECₒ[ili[1],τ,..] .= ECᵢ[oli,τ,..]
        ECₒ[ili[2],τ,..] .= ECᵢ[oli,τ,..]
    end

    for div in divs
        ili = inlinkids(net,div)[1]
        oli = outlinkids(net,div)

        for li in oli
            costupdate!(li,τ)
        end

        sr = [SR[div][i][τ,..] for i in 1:2]
        ECₒ[ili,τ,..] .= sr[1] .* ECᵢ[oli[1],τ,..] .+ sr[2] .* ECᵢ[oli[2],τ,..]
    end
end
#end
