### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ ea0dc62a-a3e3-11ed-11c7-415e14e982b7
begin
	using Pkg
	Pkg.activate("../Project.toml")
end

# ╔═╡ 1fd386df-2836-4cd8-8a24-c3118d221790
begin
	using Roots
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
	λ = 0.02
	ν = 5e-3
	γ = λ/ν
	ρ = 0.5
	μ = 100
	α = 0.8
	f = 0.1

	trueparam = (ν, γ, ρ, μ, α, f)
	
	numsamples = 1e7
end

# ╔═╡ d5101c3d-cd36-4e1c-a2d3-5261ae997ae0
md"""
We simulate the noise process by generating a random number of cells, each of which contribute a negative binomial draw. These are then from with a poisson distribution.
"""

# ╔═╡ ced52ffe-ee99-44b0-b8b7-193b24298f17
begin
	Burstdist = Poisson(γ)

	# takes arguments r, p. α = 1/r and p=1/(1+αμ)
	RNAdist = NegativeBinomial(1/α, 1/(1+α*μ))
	
	Cmean = μ*γ*ν*ρ

	Clim = 100*Cmean
	dC = Clim/1000
	Dlist = 0:dC:Clim

	slim = log((α*μ+1)/(α*μ))
	ds = 0.95*slim/1000
	slist = 0:ds:(0.95*slim)

	RNA = [sum(rand(RNAdist, rand(Burstdist))) for i in 1:numsamples]
	RNAnoise = [rand(Binomial(R, ρ*ν)) for R in RNA]
	
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
	function gExp(z)
		1/(1-(α*μ)*(z-1))^(1/α)
	end

	function gCont(z)
		exp(γ*(gExp((1+ν*ρ*(z-1))) - 1))
	end

	function gCPoiss(z)
		exp(γ*μ*ν*ρ*(z-1))
	end
end

# ╔═╡ 1967c720-5d8f-4691-a184-0e1b410b0f7b
md"""
By expanding the generating functions we can produce histograms
"""

# ╔═╡ c031d021-b2b6-4875-91e8-f5073ef480fd
begin
	countlim = maximum(keys(noisefreq))-3
	
	ana = taylor_expand(gCont, order = countlim).coeffs
	emp = [haskey(noisefreq, k) ? noisefreq[k] : 0 for k in 0:countlim]
	poisson = taylor_expand(gCPoiss, order = countlim).coeffs
	
	groupedbar(0:countlim, [emp ana poisson], label=false, xlabel="RNA Count", ylabel="P(D)", xlim = [-0.7, countlim-0.5], xticks=0:10, title="ν=$(ν), γ=$(γ), ρ=$(ρ), μ=$(μ), α=$(α)", size=(600, 350), tickfontsize=12, labelfontsize=14, lw=0)
	
	groupedbar!(0:countlim, [emp ana poisson], label=["Empirical" "Exact" "Poisson"], yaxis=:log10, inset=(1, bbox(0.4, 0.23, 0.6, 0.72, :bottom, :left)), subplot = 2, xlabel = "RNA Count", xlim = (-1, countlim-0.5), xticks=0:2:10, ylim =(1e-11, 1), tickfontsize=10, legendfontsize=10, legend=:bottomleft, lw=0)

	# savefig("../Plots/ExactHist_Poisson.png")
end

# ╔═╡ 1dfb0517-8529-4715-a0d7-018cead56ecf
ana

# ╔═╡ 0cc8fd42-a66b-4a55-bf03-4847163fb855
md"""
This seems to be working well. Let's make some functions to generate sample data.
"""

# ╔═╡ 372374aa-9ecf-41da-92a0-6e1886a92118
begin
	# generate N samples from sequencing data
	function sampleseq(N, ν, γ, ρ, μ, α, f)
		burstdist = Poisson(f*γ)
		RNAdist = NegativeBinomial(1/α, 1/(1+α*μ))
		

		samples = []
		for i in 1:N
			RNA = sum(rand(RNAdist, rand(burstdist)))
			contaminant = rand(Binomial(RNA, ν*ρ))
			if rand() < f
				E = rand(RNAdist)
				push!(samples, rand(Binomial(E, ρ)) + contaminant)
			else
				push!(samples, contaminant)
			end
		end
		return samples
	end

	numsimcounts= 10000
	simcounts = sampleseq(numsimcounts, trueparam...)

	# Generate samples from expression
	function sampleexp(N, μ, α, ρ)

		RNAdist = NegativeBinomial(1/α, 1/(1+α*μ*ρ))
		samples = []
		
		for i in 1:N
			push!(samples, rand(RNAdist))
		end
		return samples
	end

	# Generate noise samples
	function samplenoise(N, ν, γ, ρ, μ, α, f)
		burstdist = Poisson(f*γ)
		RNAdist = NegativeBinomial(1/α, 1/(1+α*μ))
		

		samples = []
		for i in 1:N
			RNA = sum(rand(RNAdist, rand(burstdist)))
			push!(samples, rand(Binomial(RNA, ν*ρ)))
		end
		return samples
	end

	simexpcounts = sampleexp(numsimcounts, μ, α, ρ)
	simnoisecounts = samplenoise(numsimcounts, trueparam...)

	histogram(simcounts, ylim=(0, Int(0.02*numsimcounts)))
end

# ╔═╡ 4922aef3-599d-421a-9634-9258aa161a00
md"""
### Method of Moments
"""

# ╔═╡ 74befa96-96bd-46d8-ac1b-6a6252fb0137
md"""
We begin with a method of moments test for the negative binomial distribution. The moment conditions for these are

$$\langle E\rangle = \mu,\quad \langle E^2\rangle = \mu + \mu^2(1+\alpha)$$
from which we get

$$\mu = \langle E\rangle $$
$$\alpha = \frac{\langle E \rangle^2 - \mu}{\mu^2} -1 = \frac{\langle E\rangle^2 - \langle E\rangle}{\langle E\rangle^2} -1$$

Note that for sequenced counts, we will take $\mu \to \rho \mu$.
"""

# ╔═╡ b5484c36-a167-4c47-9f62-5c23ec436641
md"""
For the contaminant distribution,

$$\langle C\rangle = f\gamma\nu\rho\mu, \quad \langle C^2\rangle = f\gamma\nu\rho\mu + (f\gamma\nu\rho\mu)^2 -f\gamma\mu(\rho\nu)^2+ f\gamma (\rho\nu)^2\langle E^2\rangle$$
$$\langle C^2\rangle = f\gamma\nu\rho\mu + (f\gamma\nu\rho\mu)^2 + f\gamma \nu^2(\rho\mu)^2(1+\alpha)$$
$$\langle C^2\rangle = \langle C\rangle + \langle C\rangle^2 + f\gamma \nu^2(\rho\mu)^2(1+\alpha)$$

These give us the simple relationships

$$f\gamma\nu = \frac{\langle C\rangle}{\langle E\rangle}$$
$$f\gamma\nu^2 = \frac{\langle C^2\rangle - \langle C\rangle - \langle C\rangle^2}{\langle E^2\rangle - \langle E\rangle}$$


from which we get,

$$\gamma = \frac{1}{f}\frac{\langle C\rangle^2}{\langle E\rangle^2}\frac{\langle E^2\rangle - \langle E\rangle}{\langle C^2\rangle - \langle C\rangle^2 - \langle C\rangle}$$
$$\nu = \frac{\langle E\rangle}{\langle C\rangle}\frac{\langle C^2\rangle - \langle C\rangle^2 - \langle C\rangle}{\langle E^2\rangle - \langle E\rangle}$$
"""

# ╔═╡ 99fab631-16a2-42ce-ad18-e3a60c985667
begin
	reducedparam = [ν, γ, ρ*μ, α, f]
	poissonparam = [ν*γ, ρ*μ, α, f]
end

# ╔═╡ 31f2e4b6-5e1b-4277-86e7-55f3c03d1450
begin
	function MoMtest(countdict, thresh)
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

		γi = (1/fi)*(C1i/E1i)^2*(E2i-E1i)/(C2i-C1i-C1i^2)
		νi = E1i/C1i*(C2i-C1i^2-C1i)/(E2i-E1i)
		

		[νi, γi, ρμi, αi, fi]
	end
	
	simcountdict = countmap(simcounts)
	initial = MoMtest(simcountdict, 5)
	initialpois = [initial[1]*initial[2], initial[3:end]...]
	# Not bad!
end

# ╔═╡ 63797ab2-646e-4be0-9352-c17368f3bc6b
md"""
### Maximum Likelihood Inference
"""

# ╔═╡ 9946ac22-aedd-43fa-b732-01a057f714ca
begin
	function GNB(z, μ, α)
		1/(1-(α*μ)*(z-1))^(1/α)
	end
	
	function Gpois(z, λ)
		exp(λ*(z-1))
	end

	gE(z) = GNB(z, ρ*μ, α)
	gEb(z) = GNB(z, μ*ν*ρ, α)
	gC(z) = Gpois(gEb(z), f*γ)	
	
	taylor_expand(z->gC(z)*(f*gE(z) + 1-f), order = 1000)
end

# ╔═╡ 0a6cd7d5-90ba-4232-9ead-c3357a3e0cc2
begin
	function fitgene_MLE(countdict, lower, upper, initial)
		maxcount = maximum(keys(countdict))
		sumcounts = sum(nc for (c, nc) in countdict)
	
		function objective(param)
			ν, γ, ρμ, α, f = param
			gE(z) = GNB(z, ρμ, α)
			gEν(z) = GNB(z, ρμ*ν, α)
			gC(z) = Gpois(gEν(z), f*γ)
			gseq(z) = gC(z)*(f*gE(z) + (1-f))
	
			plist = taylor_expand(gseq, order = maxcount).coeffs
			# for numerical errors at very small probabilities
			plist[plist .<= 0] .= eps(typeof(plist[1]))
			
			sum(-nc*log(plist[c+1]) for (c, nc) in countdict) / sumcounts
	
		end
	
		Optim.optimize(objective, lower, upper, initial)
	end
	
	lower = [0., 0., 0., 0., 0.]
	upper = [1., 1e3, 1e5, 100, 1.]
	
	res = fitgene_MLE(simcountdict, lower, upper, initial)
end

# ╔═╡ c70a53df-c408-41c4-a0a3-683c1d65c891
νf, γf, ρμf, αf, ff = res.minimizer

# ╔═╡ f56dc1b7-b6cd-4d67-bc45-19450c604d97
abs.(res.minimizer .- reducedparam)./reducedparam

# ╔═╡ 48b37d5a-cf46-4eaa-9b2e-14b2115ef0c7
begin
	gEf(z) = GNB(z, ρμf, αf)
	gEνf(z) = GNB(z, ρμf*ν, αf)
	gCf(z) = Gpois(gEνf(z), ff*γf)
	gseqf(z) = gCf(z)*(f*gEf(z) + (1-f))

	maxcounts = maximum(keys(simcountdict))
	
	
	pcont = taylor_expand(gCf, order = maxcounts).coeffs
	pcontex = taylor_expand(z->gCf(z)*gEf(z), order = maxcounts).coeffs
	pseqlist = taylor_expand(gseqf, order=maxcounts).coeffs

	plot(0:maxcounts, pseqlist, ylim = (0, 0.02))
	plot!(0:maxcounts, [get(simcountdict, c, 0)/numsimcounts for c in 0:maxcounts])
end

# ╔═╡ b0672e37-530f-4da1-ad73-c24e1170bc43
begin
	xlim = (0, 140)
	xticks = [0, 100, 200]
	ylim = (0, 0.004)
	yticks = [0.0, 0.002, 0.004]
	
	plot(0:maxcounts, [get(simcountdict, c, 0)/numsimcounts for c in 0:maxcounts], color = 1, lw=3, label="Counts")

	plot!(0:maxcounts, (1-ff).*pcont, lw=3, label="Contaminant")
	plot!(0:maxcounts, ff.*pcontex, lw=3, label="Exp+Cont", c=5)
	
	# plot!(plotcounts, pseqlist.+4e-5, lw=3, label="Total")
	
	plot!(xlim=xlim, xticks = xticks, ylim=ylim, yticks=yticks, size=(600, 300), tickfontsize=16, legendfontsize=12)

	# savefig("../Plots/Simulated_ZoomCounts.png")

end

# ╔═╡ f5424f56-6d78-41f1-a895-897377546423
md"""
### Generating function fits
"""

# ╔═╡ 4c820b67-23af-4a04-8e4d-06f3d9e70402
# begin
# 	Gemp(z) = sum(nc*z^c for (c, nc) in simcountdict)/numsimcounts

# 	function objectiveg(param)
# 		ν, γ, ρμ, α, f = param
# 		dz = 0.01
# 		gE(z) = GNB(z, ρμ, α)
# 		gEν(z) = GNB(z, ρμ*ν, α)
# 		gC(z) = Gpois(gEν(z), f*γ)
# 		gseq(z) = gC(z)*(f*gE(z) + (1-f))
	
# 		sum((gseq(z)-Gemp(z))^2*dz for z in dz:dz:1)
# 	end
	
# 	objectiveg(initial)
# end

# ╔═╡ 9a692ce9-8c37-4426-ae60-a4109eb43a5d
# resg = Optim.optimize(objectiveg, lower, upper, initial)

# ╔═╡ 318925af-739a-4893-a660-76e239256aff
# resg.minimizer

# ╔═╡ 76649067-a0ce-4aff-b4e8-5c1ea7274d9d
# abs.(resg.minimizer .- reducedparam)./reducedparam

# ╔═╡ 6807e4e4-5887-4a13-bfd0-3ea0eb70e7f3
# begin
# 	νg, γg, ρμg, αg, fg = resg.minimizer
# 	gEg(z) = GNB(z, ρμg, αg)
# 	gEνg(z) = GNB(z, ρμg*νg, αg)
# 	gCg(z) = Gpois(gEνg(z), fg*γg)
# 	gseqg(z) = gCg(z)*(fg*gEg(z) + (1-fg))
# end

# ╔═╡ 2f9a15f4-d1ab-4baf-a422-92d24837c487
# begin
	
# 	plot(0:0.01:1, [Gemp(z) for z in 0:0.01:1])
# 	plot!(0:0.01:1, [gseqg(z) for z in 0:0.01:1])
# end

# ╔═╡ fdb01519-b8d3-475f-8b28-60e84813d329
# begin
# 	pseqglist = taylor_expand(gseqg, order=maxcounts).coeffs
	
# 	plot(0:maxcounts, pseqglist, ylim = (0, 0.02))
# 	plot!(0:maxcounts, [get(simcountdict, c, 0)/numsimcounts for c in 0:maxcounts])
# end

# ╔═╡ beffd53e-d905-4448-b089-b03a373b1d35
# begin
# 	zvalues = [0.]
# 	nz = 20
# 	Δg = (Gemp(1) - Gemp(0))/nz
	
# 	for n in 1:nz
# 		push!(zvalues, fzero(z->Gemp(z) - (Gemp(0)+n*Δg), 0, 1))
# 	end
# end

# ╔═╡ 0620ffd5-c094-4649-a220-aee05bc9b30e
# begin
# 	function objectivegz(param)
# 		ν, γ, ρμ, α, f = param
# 		gE(z) = GNB(z, ρμ, α)
# 		gEν(z) = GNB(z, ρμ*ν, α)
# 		gC(z) = Gpois(gEν(z), f*γ)
# 		gseq(z) = gC(z)*(f*gE(z) + (1-f))
	
# 		Gemp(z) = sum(nc*z^c for (c, nc) in simcountdict)/numsimcounts
		
# 		sum(((gseq(z1)-Gemp(z1))^2+(gseq(z2)-Gemp(z2))^2)/2*(z2-z1)  for (z1,z2) in zip(zvalues[1:end-1], zvalues[2:end]))
# 	end
	
# 	objectivegz(initial)
# end

# ╔═╡ eafcf8a8-ba24-4cb6-bdbf-0352c21a008f
# resgz = Optim.optimize(objectivegz, lower, upper, initial)

# ╔═╡ 0b2b4594-8d73-42d2-9d93-d6487d82e048
# resgz.minimizer

# ╔═╡ b47265bf-d6f3-409d-a97b-9f535c94225b
# abs.(resgz.minimizer .- reducedparam)./reducedparam

# ╔═╡ 3a67dfd7-ef4d-4b24-b640-3c62bcdc942e
# begin
# 	νgz, γgz, ρμgz, αgz, fgz = resg.minimizer
# 	gEgz(z) = GNB(z, ρμgz, αgz)
# 	gEνgz(z) = GNB(z, ρμgz*νgz, αgz)
# 	gCgz(z) = Gpois(gEνgz(z), fgz*γgz)
# 	gseqgz(z) = gCgz(z)*(fgz*gEgz(z) + (1-fgz))
# end

# ╔═╡ 660f0812-63c7-41cc-a576-7a296d3c8289
# begin
# 	pseqgzlist = taylor_expand(gseqgz, order=maxcounts).coeffs
	
# 	plot(0:maxcounts, pseqgzlist, ylim = (0, 0.02))
# 	plot!(0:maxcounts, [get(simcountdict, c, 0)/numsimcounts for c in 0:maxcounts])
# end

# ╔═╡ 0361d0bc-5ec1-4f73-9567-605cd0d364bc
# begin
# 	lowerpois = [0., 0., 0., 0.]
# 	upperpois = [1., 1e5, 100, 1.]
# 	function objectivegzpois(param)
# 		λ, ρμ, α, f = param
# 		gE(z) = GNB(z, ρμ, α)
# 		gC(z) = Gpois(z, λ*f*ρμ)
# 		gseq(z) = gC(z)*(f*gE(z) + (1-f))
	
# 		Gemp(z) = sum(nc*z^c for (c, nc) in simcountdict)/numsimcounts
		
# 		sum(((gseq(z1)-Gemp(z1))^2+(gseq(z2)-Gemp(z2))^2)/2*(z2-z1)  for (z1,z2) in zip(zvalues[1:end-1], zvalues[2:end]))
# 	end
	
# 	objectivegzpois(initialpois)
# end

# ╔═╡ 6cf630e3-8ebe-404a-acf3-42de1f7f23a5
# resgzpois = Optim.optimize(objectivegzpois, lowerpois, upperpois, initialpois)

# ╔═╡ 0e290a6c-31be-49fd-ad2b-03b58c89d918
# resgzpois.minimizer

# ╔═╡ 70f02cec-4caa-4bbf-bdf0-c46e7ced0b86
# abs.(resgzpois.minimizer .- poissonparam)./poissonparam

# ╔═╡ a6faa73a-df7c-4d8b-8f24-10ca8c38b8ae
# begin
# 	λgzp, ρμgzp, αgzp, fgzp = resgzpois.minimizer
# 	gEgzp(z) = GNB(z, ρμgzp, αgzp)
# 	gCgzp(z) = Gpois(z, λgzp*fgzp*ρμgzp)
# 	gseqgzp(z) = gCgzp(z)*(fgzp*gEgzp(z) + (1-fgzp))
# end

# ╔═╡ c080f143-fc5e-4c30-822d-1021d7036d4b
# begin
# 	pseqgzplist = taylor_expand(gseqgzp, order=maxcounts).coeffs
	
# 	plot(0:maxcounts, pseqgzplist, ylim = (0, 0.02))
# 	plot!(0:maxcounts, [get(simcountdict, c, 0)/numsimcounts for c in 0:maxcounts])
# end

# ╔═╡ 276d9af3-9cce-483e-9b8b-087131acfc8a
md"""
### Convergence
"""

# ╔═╡ 424f59c4-b537-485f-b6fc-6ebc9ee0eb6f
begin
	
	numcountvec = [1000, 2000, 4000, 6000, 8000, 10000]
	errν = [[] for n in numcountvec]
	errγ = [[] for n in numcountvec]
	errμ = [[] for n in numcountvec]
	errα = [[] for n in numcountvec]
	errf = [[] for n in numcountvec]
	errνγ = [[] for n in numcountvec]
	
	for (nidx, numsimcounts) in enumerate(numcountvec)
		for _ in 1:100
			simcounts = sampleseq(numsimcounts, trueparam...)
			simcountdict = countmap(simcounts)
			initial = MoMtest(simcountdict, 5)
			res = fitgene_MLE(simcountdict, lower, upper, initial)
			νf, γf, ρμf, αf, ff = res.minimizer
	
			push!(errν[nidx], (νf-ν)/ν)
			push!(errγ[nidx], (γf-γ)/γ)
			push!(errμ[nidx], (ρμf-ρ*μ)/(ρ*μ))
			push!(errα[nidx], (αf-α)/α)
			push!(errf[nidx], (ff-f)/f)
			push!(errνγ[nidx], (νf*γf - ν*γ)/(ν*γ))
		end
	end

end

# ╔═╡ ca77d91d-8ff2-4a4b-b330-96889e7affad
begin
	errvar = errf
	scatter()
	for (idx, nc) in enumerate(numcountvec)
		scatter!([numcountvec[idx] for _ in errν[idx]], errvar[idx], markerstrokewidth=0, color = 1, alpha=0.2, label = false)
	end
	plot!(numcountvec, [mean(e) for e in errvar], label = "Mean")
	
end

# ╔═╡ 4ba17ac7-aeeb-4ea7-85bc-a02219570bb3
histogram(errvar[1])

# ╔═╡ 6e61189c-f823-480b-afb7-5392a2f0ecab
begin
	
	plot(numcountvec, [var(e) for e in errvar])
	scatter!(numcountvec, [var(e) for e in errvar], markerstrokewidth=0, color=1)

	plot!(yaxis = :log10, xaxis=:log10)
end

# ╔═╡ 586b4f00-4b85-4ca0-97a1-8eaed9b1fb95
var(errvar[1]), var(errvar[end])

# ╔═╡ 21f8e18f-50cb-442a-9cb9-4088063f1a98
(log(var(errvar[end])) - log(var(errvar[1])))/(log(numcountvec[end]) - log(numcountvec[1]))

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
# ╠═1dfb0517-8529-4715-a0d7-018cead56ecf
# ╟─0cc8fd42-a66b-4a55-bf03-4847163fb855
# ╠═372374aa-9ecf-41da-92a0-6e1886a92118
# ╟─4922aef3-599d-421a-9634-9258aa161a00
# ╟─74befa96-96bd-46d8-ac1b-6a6252fb0137
# ╟─b5484c36-a167-4c47-9f62-5c23ec436641
# ╠═99fab631-16a2-42ce-ad18-e3a60c985667
# ╠═31f2e4b6-5e1b-4277-86e7-55f3c03d1450
# ╟─63797ab2-646e-4be0-9352-c17368f3bc6b
# ╠═9946ac22-aedd-43fa-b732-01a057f714ca
# ╠═0a6cd7d5-90ba-4232-9ead-c3357a3e0cc2
# ╠═c70a53df-c408-41c4-a0a3-683c1d65c891
# ╠═f56dc1b7-b6cd-4d67-bc45-19450c604d97
# ╠═48b37d5a-cf46-4eaa-9b2e-14b2115ef0c7
# ╠═b0672e37-530f-4da1-ad73-c24e1170bc43
# ╟─f5424f56-6d78-41f1-a895-897377546423
# ╠═4c820b67-23af-4a04-8e4d-06f3d9e70402
# ╠═9a692ce9-8c37-4426-ae60-a4109eb43a5d
# ╠═318925af-739a-4893-a660-76e239256aff
# ╠═76649067-a0ce-4aff-b4e8-5c1ea7274d9d
# ╠═6807e4e4-5887-4a13-bfd0-3ea0eb70e7f3
# ╠═2f9a15f4-d1ab-4baf-a422-92d24837c487
# ╠═fdb01519-b8d3-475f-8b28-60e84813d329
# ╠═beffd53e-d905-4448-b089-b03a373b1d35
# ╠═0620ffd5-c094-4649-a220-aee05bc9b30e
# ╠═eafcf8a8-ba24-4cb6-bdbf-0352c21a008f
# ╠═0b2b4594-8d73-42d2-9d93-d6487d82e048
# ╠═b47265bf-d6f3-409d-a97b-9f535c94225b
# ╠═3a67dfd7-ef4d-4b24-b640-3c62bcdc942e
# ╠═660f0812-63c7-41cc-a576-7a296d3c8289
# ╠═0361d0bc-5ec1-4f73-9567-605cd0d364bc
# ╠═6cf630e3-8ebe-404a-acf3-42de1f7f23a5
# ╠═0e290a6c-31be-49fd-ad2b-03b58c89d918
# ╠═70f02cec-4caa-4bbf-bdf0-c46e7ced0b86
# ╠═a6faa73a-df7c-4d8b-8f24-10ca8c38b8ae
# ╠═c080f143-fc5e-4c30-822d-1021d7036d4b
# ╟─276d9af3-9cce-483e-9b8b-087131acfc8a
# ╠═424f59c4-b537-485f-b6fc-6ebc9ee0eb6f
# ╠═ca77d91d-8ff2-4a4b-b330-96889e7affad
# ╠═4ba17ac7-aeeb-4ea7-85bc-a02219570bb3
# ╠═6e61189c-f823-480b-afb7-5392a2f0ecab
# ╠═586b4f00-4b85-4ca0-97a1-8eaed9b1fb95
# ╠═21f8e18f-50cb-442a-9cb9-4088063f1a98
