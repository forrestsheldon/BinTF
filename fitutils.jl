using Distributions, StatsBase
using Optim, TaylorSeries

function formcountdict(BC, BCdf, zerobackground)

	genecounts = BCdf[:, BC]
	countdict = countmap(genecounts)
	countdict[0] += zerobackground
	
	return countdict
end

function formcountfreq(countdict)
	counts = 0:maximum(keys(countdict))
	N = sum(values(countdict))
	countfreq = [haskey(countdict, c) ? countdict[c]/N : 0 for c in counts]

	counts, countfreq
end

function MoMinitial(countdict, thresh)
    E1i = 0
    E2i = 0
    nxi = 0
    C1i = 0
    C2i = 0
    nyi = 0
    
    for (c, nc) in countdict
        if c < thresh
            C1i += nc*c
            C2i += nc*c^2
            nyi += nc
        else
            E1i += nc*c
            E2i += nc*c^2
            nxi += nc
        end
    end
    E1i = E1i/nxi
    E2i = E2i/nxi
    C1i = C1i/nyi
    C2i = C2i/nyi

    fi = nxi / (nxi + nyi)
    
    ρμi = E1i
    αi = (E2i - E1i)/E1i^2 - 1

    if C1i != 0
        γi = (1/fi)*(C1i/E1i)^2*(E2i-E1i)/(C2i-C1i-C1i^2)
        νi = E1i/C1i*(C2i-C1i^2-C1i)/(E2i-E1i)
    else
        νi = 0
        γi = 0
    end
    

    [νi, γi, ρμi, αi, fi]
end

function GNB(z, μ, α)
    1/(1-(α*μ)*(z-1))^(1/α)
end

function Gpois(z, λ)
    exp(λ*(z-1))
end


function fitgene_MLE(countdict, lower, upper, initial)
    maxcount = maximum(keys(countdict))
    sumcounts = sum(nc for (c, nc) in countdict)

    function objective(param)
        ν, γ, ρμ, α, f = param

        gex(z) = GNB(z, ρμ, α)
        gexν(z) = GNB(z, ρμ*ν, α)
        gcont(z) = Gpois(gexν(z), f*γ)
        gseq(z) = gcont(z)*(f*gex(z) + (1-f))

        plist = taylor_expand(gseq, order = maxcount).coeffs
        plist[plist .< 0] .= eps(typeof(plist[1]))
        
        sum(-nc*log(plist[c+1]) for (c, nc) in countdict) / sumcounts

    end

    Optim.optimize(objective, lower, upper, initial)
end

# function fitgene_MLE_νγ(countdict, νγ_const, lower, upper, initial)
#     maxcount = maximum(keys(countdict))
#     sumcounts = sum(nc for (c, nc) in countdict)

#     function objective(param)
#         ν, μ, r, f = param
#         gν(z) = Gpois(z, ν)
#         gex(z) = GNB(z, μ, r)
#         gno(z) = Gpois(gex(gν(z)), f*νγ_const/ν)
#         gnoex(z) = gex(z)*gno(z)
#         gseq(z) = f*gnoex(z) + (1-f)*gno(z)

#         plist = taylor_expand(gseq, order = maxcount).coeffs
#         plist[plist .< 0] .= eps(typeof(plist[1]))

#         sum(-nc*log(plist[c+1]) for (c, nc) in countdict) / sumcounts

#     end

#     Optim.optimize(objective, lower, upper, initial)
# end

function fitgene_MLE_reg(countdict, γ_νγ_regparam, lower, upper, initial)
    maxcount = maximum(keys(countdict))
    sumcounts = sum(nc for (c, nc) in countdict)
    γ_med, γ_MAD, νγ_med, νγ_MAD = γ_νγ_regparam

	function objective(param)
        ν, γ, ρμ, α, f = param
        ν = max(1e-20, ν)
        γ = max(1e-20, γ)

        gex(z) = GNB(z, ρμ, α)
        gexν(z) = GNB(z, ρμ*ν, α)
        gcont(z) = Gpois(gexν(z), f*γ)
        gseq(z) = gcont(z)*(f*gex(z) + (1-f))


        plist = taylor_expand(gseq, order = maxcount).coeffs
        plist[plist .< 0] .= eps(typeof(plist[1]))

        
        sum(-nc*log(plist[c+1]) for (c, nc) in countdict) / sumcounts + (ν*γ - νγ_med)^2/(2*(1.48*νγ_MAD)^2)+ (γ - γ_med)^2/(2*(1.48*γ_MAD)^2)

    end

    Optim.optimize(objective, lower, upper, initial)
end

# This is a modified version of fitgene that takes fixed νγ and f parameters and then fits based on ν, μ, r
function fitgene_MLE_regf(celldict, regparam, f_const, lower, upper, initial)
    maxcount = maximum(keys(celldict))
    numcellcounts = sum(nc for (c, nc) in celldict)
    νγ_med, νγ_MAD, hprbl_med, hprbl_MAD = regparam


    function objectivefull(param)
        ν, γ, μ, r = param
        
        gν(z) = Gpois(z, ν)
        gex(z) = GNB(z, μ, r)
        gno(z) = Gpois(gex(gν(z)), f_const*γ)
        gnoex(z) = gex(z)*gno(z)
        gseq(z) = f_const*gnoex(z) + (1-f_const)*gno(z)

        plist = taylor_expand(gseq, order = maxcount).coeffs
        plist[plist .< 0] .= eps(typeof(plist[1]))

        hprbl = asinh((ν-γ)/(2*sqrt(ν*γ)))

        sum(-nc*log(plist[c+1]) for (c, nc) in celldict) / numcellcounts + (ν*γ - νγ_med)^2/(2*(1.48*νγ_MAD)^2)+ (hprbl - hprbl_med)^2/(2*(1.48*hprbl_MAD)^2)


    end

    return Optim.optimize(objectivefull, lower, upper, initial)

end