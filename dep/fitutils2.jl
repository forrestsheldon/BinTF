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
    x1i = 0
    x2i = 0
    nxi = 0
    y1i = 0
    y2i = 0
    nyi = 0
    
    for (c, nc) in countdict
        if c < thresh
            y1i += nc*c
            y2i += nc*c^2
            nyi += nc
        else
            x1i += nc*c
            x2i += nc*c^2
            nxi += nc
        end
    end
    x1i = x1i/nxi
    x2i = x2i/nxi
    y1i = y1i/nyi
    y2i = y2i/nyi

    fi = nxi / (nxi + nyi)
    
    μi = x1i
    ri = 1/(x2i/x1i^2 - 1)
    γi = (1/fi)*(x2i/x1i^2)*y1i^2/(y2i-y1i^2)
    νi = x1i/x2i*(y2i-y1i^2)/y1i
    

    [νi, γi, μi, ri, fi]
end

function GNB(z, μ, r)
    1/(1-(μ/r)*(z-1))^r
end

function Gpois(z, λ)
    exp(λ*(z-1))
end

function Gnoise(z, ν, γ, μ, r)
    gν(z) = Gpois(z, ν)
    gex(z) = GNB(z, μ, r)
    return Gpois(gex(gν(z)), γ)
end

function fitgene_MLE(countdict, lower, upper, initial)
    maxcount = maximum(keys(countdict))
    sumcounts = sum(nc for (c, nc) in countdict)

    function objective(param)
        ν, γ, μ, r, f = param
        gν(z) = Gpois(z, ν)
        gex(z) = GNB(z, μ, r)
        gno(z) = Gpois(gex(gν(z)), f*γ)
        gnoex(z) = gex(z)*gno(z)
        gseq(z) = f*gnoex(z) + (1-f)*gno(z)

        plist = taylor_expand(gseq, order = maxcount).coeffs
        plist[plist .< 0] .= eps(typeof(plist[1]))
        
        sum(-nc*log(plist[c+1]) for (c, nc) in countdict) / sumcounts

    end

    Optim.optimize(objective, lower, upper, initial)
end

function fitgene_MLE_νγ(countdict, νγ_const, lower, upper, initial)
    maxcount = maximum(keys(countdict))
    sumcounts = sum(nc for (c, nc) in countdict)

    function objective(param)
        ν, μ, r, f = param
        gν(z) = Gpois(z, ν)
        gex(z) = GNB(z, μ, r)
        gno(z) = Gpois(gex(gν(z)), f*νγ_const/ν)
        gnoex(z) = gex(z)*gno(z)
        gseq(z) = f*gnoex(z) + (1-f)*gno(z)

        plist = taylor_expand(gseq, order = maxcount).coeffs
        plist[plist .< 0] .= eps(typeof(plist[1]))

        sum(-nc*log(plist[c+1]) for (c, nc) in countdict) / sumcounts

    end

    Optim.optimize(objective, lower, upper, initial)
end

function fitgene_MLE_reg(countdict, regparam, lower, upper, initial)
    maxcount = maximum(keys(countdict))
    sumcounts = sum(nc for (c, nc) in countdict)
    νγ_med, νγ_MAD, hprbl_med, hprbl_MAD = regparam

	function objective(param)
        ν, γ, μ, r, f = param
        ν = max(1e-20, ν)
        γ = max(1e-20, γ)
        gν(z) = Gpois(z, ν)
        gex(z) = GNB(z, μ, r)
        gno(z) = Gpois(gex(gν(z)), f*γ)
        gnoex(z) = gex(z)*gno(z)
        gseq(z) = f*gnoex(z) + (1-f)*gno(z)

        plist = taylor_expand(gseq, order = maxcount).coeffs
        plist[plist .< 0] .= eps(typeof(plist[1]))

        hprbl = asinh((ν-γ)/(2*sqrt(ν*γ)))
        
        sum(-nc*log(plist[c+1]) for (c, nc) in countdict) / sumcounts + (ν*γ - νγ_med)^2/(2*(1.48*νγ_MAD)^2)+ (hprbl - hprbl_med)^2/(2*(1.48*hprbl_MAD)^2)

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