using Distributions, StatsBase
using Optim, TaylorSeries

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

function Gex(z, μ, r)
    1/(1-(μ/r)*(z-1))^r
end

function Gpois(z, λ)
    exp(λ*(z-1))
end

function Gbern(z, ν)
    1 + ν*(z-1)
end

function Gnoise(z, ν, γ, μ, r, gν, gex)
    Gpois(gex(gν(z)), γ)
end

function Gnoise(z, ν, γ, μ, r)
    gν(z) = Gbern(z, ν)
    gex(z) = Gex(z, μ, r)
    return Gpois(gex(gν(z)), γ)
end

function Gnoex(z, ν, γ, μ, r, gν, gex, gno)
    gno(z)*gex(z)
end

function Gnoex(z, ν, γ, μ, r)
    Gnoise(z, ν, γ, μ, r) * Gex(z, μ, r)
end

function Gseq(z, ν, γ, μ, r, f, gν, gex, gno)
    f*gno(z)*gex(z) + (1-f)*gno(z)
end

function Gseq(z, ν, γ, μ, r, f)
    gν(z) = Gbern(z, ν)
    gex(z) = Gex(z, μ, r)
    gno(z) = Gpois(gex(gν(z)), f*γ)
    f*gno(z)*gex(z) + (1-f)*gno(z)
end

function fitgene_MLE(countdict, lower, upper, initial)
    maxcount = maximum(keys(countdict))
    sumcounts = sum(nc for (c, nc) in countdict)

    function objective(param)
        ν, γ, μ, r, f = param
        # gν(z) = Gpois(z, ν)
        gν(z) = 1 + ν*(z-1)
        gex(z) = Gex(z, μ, r)
        gno(z) = Gpois(gex(gν(z)), f*γ)
        gseq(z) = Gseq(z, param..., gν, gex, gno)

        plist = taylor_expand(gseq, order = maxcount).coeffs
        
        sum(-nc*log(plist[c+1]) for (c, nc) in countdict) / sumcounts

    end

    Optim.optimize(objective, lower, upper, initial)
end

