
        sval += try
            Nᵢ[i,t+1-l] - Nₒ[i,t]
        catch BoundsError
            -Nₒ[i,t]
        end

        rval += try
            Nₒ[i,t+1-Int(round(l/δ))] - Nᵢ[i,t]
        catch BoundsError
            -Nᵢ[i,t]
        end

        sval = round(sval, digits=ROUND_DIGITS)
        rval = round(rval, digits=ROUND_DIGITS)
        if (sval < 0.) || (rval < 0.)
            sval = (sval < 0.) ? 0. : sval
            rval = (rval < 0.) ? 0. : rval
        end
