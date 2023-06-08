### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ ea0dc62a-a3e3-11ed-11c7-415e14e982b7
begin
	using Pkg
	Pkg.activate("Project.toml")
end

# ╔═╡ 1fd386df-2836-4cd8-8a24-c3118d221790
begin
	using Distributions, StatsBase
	using Plots, Plots.PlotMeasures, StatsPlots
	using Optim, TaylorSeries
end

# ╔═╡ 0e3a3ccb-07fb-4be6-b6e3-08b1812c0e1f
md"""
This notebook uses simulated data to develop the fitting methods used in the paper. If you are unfamiliar with Julia notebooks, you can use Shift+Enter to execute a cell
"""

# ╔═╡ 89b3fb04-ac91-4845-9738-d18df3c7e7a2
md"""
First we set some parameters to test with
"""

# ╔═╡ e93e30e9-f5f5-45e5-99bf-a1f81908c099
begin
	ν = 1e-2
	γ = 0.5
	μ = 100
	r = 2.
	f = 0.05

	trueparam = (ν, γ, μ, r, f)
	
	numsamples = 1e7

end

# ╔═╡ d5101c3d-cd36-4e1c-a2d3-5261ae997ae0
md"""
We simulate the noise process by generating a random number of cells, each of which contribute a negative binomial draw. These are then from with a poisson distribution.
"""

# ╔═╡ ced52ffe-ee99-44b0-b8b7-193b24298f17
begin
	Burstdist = Poisson(γ)
		
	RNAdist = NegativeBinomial(r, r/(r+μ))
	
	Dmean = μ*γ*ν

	Dlim = 100*Dmean
	dD = Dlim/1000
	Dlist = 0:dD:Dlim

	slim = log((μ+r)/μ)
	ds = 0.95*slim/1000
	slist = 0:ds:(0.95*slim)

	RNA = [sum(rand(RNAdist, rand(Burstdist))) for i in 1:numsamples]
	RNAnoise = [rand(Binomial(R, ν)) for R in RNA]
	
	noisecounts = countmap(RNAnoise)
	noisefreq = Dict{Int,Float64}()
	for (k,v) in noisecounts
		noisefreq[k] = v/numsamples
	end

	bar(noisefreq, label=false, xlabel = "Count", ylabel="Frequency", title="Simulated Count Frequency")
end

# ╔═╡ 1a29ef91-0e94-4209-9b67-018e3ffac41c
md"""
To compare to analytics we define the generating functions we will need.
"""

# ╔═╡ a8e62d89-d929-4d8c-95a4-9a0a9d0cbf39
begin
	function gC(z)
		1/(1-(μ/r)*(z-1))^r
	end

	function gD(z)
		exp(γ*(gC(exp(ν*(z-1))) - 1))
	end

	function gDPoiss(z)
		exp(γ*μ*ν*(z-1))
	end

end

# ╔═╡ 1967c720-5d8f-4691-a184-0e1b410b0f7b
md"""
By expanding the generating functions we can produce histograms
"""

# ╔═╡ c031d021-b2b6-4875-91e8-f5073ef480fd
begin
	countlim = maximum(keys(noisefreq))-4
	
	ana = taylor_expand(gD, order = countlim).coeffs
	emp = [haskey(noisefreq, k) ? noisefreq[k] : 0 for k in 0:countlim]
	poisson = taylor_expand(gDPoiss, order = countlim).coeffs
	
	groupedbar(0:countlim, [emp ana poisson], label=false, xlabel="RNA count", ylabel="P(D)", xlim = [-0.7, countlim-4.5], xticks=0:10, title="ν=$(ν), γ=$(γ), μ=$(μ), r=$(r)", size=(600, 350), tickfontsize=12)
	groupedbar!(0:countlim, [emp ana poisson], label=["Empirical" "Exact" "Poisson"], yaxis=:log10, inset=(1, bbox(0.31, 0.25, 0.69, 0.72, :bottom, :left)), subplot = 2, xlabel = "RNA count", xlim = (-1, countlim-4.5), xticks=0:2:10, ylim =(1e-11, 1), tickfontsize=12, legendfontsize=10, legend=:bottomleft)

	# savefig("./Plots/ExactHist_Poisson.png")
end

# ╔═╡ 0cc8fd42-a66b-4a55-bf03-4847163fb855
md"""
This seems to be working well. Let's make some functions to generate sample data.
"""

# ╔═╡ 372374aa-9ecf-41da-92a0-6e1886a92118
begin
	# generate N samples from sequencing data
	function sampleseq(N, ν, γ, μ, r, f)
		burstdist = Poisson(f*γ)
		RNAdist = NegativeBinomial(r, r/(r+μ))
		

		samples = []
		for i in 1:N
			RNA = sum(rand(RNAdist, rand(burstdist)))
			noise = rand(Binomial(RNA, ν))
			if rand() < f
				push!(samples, rand(RNAdist) + noise)
			else
				push!(samples, noise)
			end
		end
		return samples
	end

	numsimcounts= 10000
	simcounts = sampleseq(numsimcounts, trueparam...)

	# Generate samples from expression
	function sampleexp(N, μ, r)

		RNAdist = NegativeBinomial(r, r/(r+μ))
		samples = []
		
		for i in 1:N
			push!(samples, rand(RNAdist))
		end
		return samples
	end

	# Generate noise samples
	function samplenoise(N, ν, γ, μ, r, f)
		burstdist = Poisson(f*γ)
		RNAdist = NegativeBinomial(r, r/(r+μ))
		

		samples = []
		for i in 1:N
			RNA = sum(rand(RNAdist, rand(burstdist)))
			push!(samples, rand(Binomial(RNA, ν)))
		end
		return samples
	end

	simexpcounts = sampleexp(numsimcounts, μ, r)
	simnoisecounts = samplenoise(numsimcounts, trueparam...)

	histogram(simcounts, ylim = (0, 30))
end

# ╔═╡ 4922aef3-599d-421a-9634-9258aa161a00
md"""
### Method of Moments
"""

# ╔═╡ 74befa96-96bd-46d8-ac1b-6a6252fb0137
md"""
We begin with a method of moments test for the negative binomial distribution. The moment conditions for these are

$$\langle X\rangle = \mu,\quad \langle X^2\rangle = \mu^2 \frac{r+1}{r}$$
from which we get

$$\mu = \langle X\rangle $$
$$ r = \frac{\mu^2}{\langle X^2\rangle - \mu^2} = \frac{\langle X\rangle^2}{\langle X^2\rangle - \langle X\rangle^2}$$
"""

# ╔═╡ b5484c36-a167-4c47-9f62-5c23ec436641
md"""
For the noise distribution,

$$\langle Y\rangle = f\gamma\mu\nu, \quad \langle Y^2\rangle = (f\gamma\mu\nu)^2 + f\gamma (\mu\nu)^2 \frac{r+1}{r}$$

$$\gamma = \frac{1}{f}\frac{\langle X^2\rangle}{\langle X\rangle^2}\frac{\langle Y\rangle^2}{\langle Y^2\rangle - \langle Y\rangle^2}$$
$$\nu = \frac{\langle Y\rangle}{f\gamma\mu} = \frac{\langle X\rangle}{\langle X^2\rangle}\frac{\langle Y^2\rangle - \langle Y\rangle^2}{\langle Y\rangle}$$
"""

# ╔═╡ 99fab631-16a2-42ce-ad18-e3a60c985667
[trueparam...]

# ╔═╡ dd3f1203-c9a8-4bff-b800-96958884b530
begin

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

	testthresh = 3
	simcountdict = countmap(simcounts)
	initial = MoMinitial(simcountdict, testthresh)
end

# ╔═╡ 63797ab2-646e-4be0-9352-c17368f3bc6b
md"""
### Maximum Likelihood Inference
"""

# ╔═╡ 9946ac22-aedd-43fa-b732-01a057f714ca
begin
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

	gν(z) = Gpois(z, ν)
	gex(z) = Gex(z, μ, r)
	gno(z) = Gpois(gex(gν(z)), f*γ)	
	
	taylor_expand(z->Gseq(z, ν, γ, μ, r, f, gν, gex, gno), order = 1000)
end

# ╔═╡ 0a6cd7d5-90ba-4232-9ead-c3357a3e0cc2
begin
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
	
	lower = [0., 0., 0., 0., 0.]
	upper = [1., 1e4, 1e4, 10, 1.]
	
	res = fitgene_MLE(simcountdict, lower, upper, initial)
end

# ╔═╡ f56dc1b7-b6cd-4d67-bc45-19450c604d97
abs.(res.minimizer .- trueparam)./trueparam

# ╔═╡ Cell order:
# ╠═ea0dc62a-a3e3-11ed-11c7-415e14e982b7
# ╠═1fd386df-2836-4cd8-8a24-c3118d221790
# ╟─0e3a3ccb-07fb-4be6-b6e3-08b1812c0e1f
# ╟─89b3fb04-ac91-4845-9738-d18df3c7e7a2
# ╠═e93e30e9-f5f5-45e5-99bf-a1f81908c099
# ╟─d5101c3d-cd36-4e1c-a2d3-5261ae997ae0
# ╠═ced52ffe-ee99-44b0-b8b7-193b24298f17
# ╟─1a29ef91-0e94-4209-9b67-018e3ffac41c
# ╠═a8e62d89-d929-4d8c-95a4-9a0a9d0cbf39
# ╠═1967c720-5d8f-4691-a184-0e1b410b0f7b
# ╠═c031d021-b2b6-4875-91e8-f5073ef480fd
# ╟─0cc8fd42-a66b-4a55-bf03-4847163fb855
# ╠═372374aa-9ecf-41da-92a0-6e1886a92118
# ╟─4922aef3-599d-421a-9634-9258aa161a00
# ╟─74befa96-96bd-46d8-ac1b-6a6252fb0137
# ╟─b5484c36-a167-4c47-9f62-5c23ec436641
# ╠═99fab631-16a2-42ce-ad18-e3a60c985667
# ╠═dd3f1203-c9a8-4bff-b800-96958884b530
# ╟─63797ab2-646e-4be0-9352-c17368f3bc6b
# ╠═9946ac22-aedd-43fa-b732-01a057f714ca
# ╠═0a6cd7d5-90ba-4232-9ead-c3357a3e0cc2
# ╠═f56dc1b7-b6cd-4d67-bc45-19450c604d97
