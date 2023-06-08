### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ dd3be7a0-cd7e-11ed-0c04-77087f8f9e51
begin
	using Pkg
	Pkg.activate("Project.toml")
end

# ╔═╡ be5dc067-4398-4511-a741-7c8115656aab
begin
	using StatsBase
	using TaylorSeries
	using CSV, DataFrames
	using Plots, Plots.PlotMeasures, StatsPlots
	using Printf
end

# ╔═╡ 9edb8c81-c6b8-4422-8ab6-ab58a3810e0e
include("fitutils.jl")

# ╔═╡ a55ade2e-de74-40f7-9b7a-2cd50a94c16a
begin
	rep = "R11"
	marker = "NG"
end

# ╔═╡ 275c2e81-026c-4a7c-9a4d-303ebbb08f56
md"""
### Load Data
"""

# ╔═╡ f7a000a2-ad36-4506-94e9-28cc451ad463
celldf = DataFrame(CSV.File("../RandomDesigns/ScreenData/NKv/bb20230124104525_NKv_R11R12R13_counts/bb20230124104525_NKv_$(rep)_$(marker)_eTF_count.txt", delim="\t"))

# ╔═╡ c4f26fa4-e9c8-4dcd-aa0f-88ecae320885
begin

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
end

# ╔═╡ d2827474-d4a4-425b-bdc8-6ded3177449b
begin
	TFlist = names(celldf)[2:end]
end

# ╔═╡ 09257638-4766-44a3-80d1-7d5332e62afe
length(TFlist)

# ╔═╡ 6ce3fb4a-bfde-4796-8e54-631e7d825e5f
begin
	TFkept = []

	# Filter for genes with at least 200 counts greater than 10
	for TF in TFlist
		if sum(celldf[:, TF] .> 5) > 100
			push!(TFkept, TF)
		end
	end
	length(TFkept)
end

# ╔═╡ 46f9cce2-9b21-49c1-92df-c964c1181c81
begin
	testTF = TFkept[1]
	testcountdict = formcountdict(testTF, celldf, 0)
	testmax = maximum(keys(testcountdict))
	testnumcounts = sum(values(testcountdict))
	testcounts, testcountfreq = formcountfreq(testcountdict)
end

# ╔═╡ 414c2324-2556-4eca-ab5a-cd94ab09e8fb
plot(testcounts, testcountfreq, ylim=(0, 0.02), label=false, xlabel="Counts", ylabel="Frequency", title=testTF)

# ╔═╡ f117944e-b6a5-40a9-ba8b-28c0578afbeb
testcountfreq[2], sum(testcountfreq[2:end]), testcountfreq[2]/sum(testcountfreq[2:end])

# ╔═╡ 39fac9de-58d8-47c3-8025-f237be4f8e56
# begin
# 	testlower = [0., 0., 0., 0., 0.]
# 	testupper = [1., 1e4, 1e4, 10, min(1, 2*sum(testcountfreq[2:end]))]
# 	testthresh = 10
	
# 	testres = fitgene_MLE(testcountdict, testlower, testupper, MoMinitial(testcountdict, testthresh))
# end

# ╔═╡ 3cb2d27e-37ae-43ea-ac38-842f5b20eeb8
# begin
# 	νfit, γfit, μfit, rfit, ffit = testres.minimizer
	
# 	pnoexlist=taylor_expand(z->ffit*Gnoex(z, νfit, ffit*γfit, μfit, rfit), order=testmax).coeffs
# 	peak = max(maximum(pnoexlist), 0.001/4)

# 	pnolist = taylor_expand(z->(1-ffit)*Gnoise(z, νfit, ffit*γfit, μfit, rfit), order=testmax).coeffs

# 	plot(testcounts, testcountfreq, ylim=(0, 0.001), xlim=(0, 0.75*testmax), label="Counts", xlabel="$(testTF) Count\nν=$(@sprintf("%.3g", νfit)) , γ=$(@sprintf("%.1f", γfit)) , μ=$(@sprintf("%.1f", μfit)) , r=$(@sprintf("%.3g", rfit)), f=$(@sprintf("%.3g", ffit))", ylabel="Fraction of Droplets")
	
# 	plot!(0:testmax, pnolist, label="Noise", lw=2)

# 	plot!(0:testmax, pnoexlist, label="Expression", lw=2)
	
# 	plot!(0:testmax, taylor_expand(z->Gseq(z, testres.minimizer...), order=testmax).coeffs .+ peak*0, label="Mixture", lw=2)
# 	plot!(0:testmax, taylor_expand(z->Gseq(z, MoMinitial(testcountdict, testthresh)...), order=testmax).coeffs, label="Initial", lw=2)
# end

# ╔═╡ 3a33e1f4-ecde-4be6-9eb3-8f83ad82321c
# testcounts[findfirst(pnoexlist .> pnolist)]

# ╔═╡ bbd1ddee-8717-4ce5-b855-3f168a33193c
# begin
# 	fitres = []
# 	thresholds = []
	
# 	thresh = 10
# 	lower = [0., 0., 0., 0., 0.]
# 	upper = [1., 1e4, 1e4, 10, 1.]
	
# 	for (TFidx, TF) in enumerate(TFkept)
# 		countdict = formcountdict(TF, celldf, 0)
# 		countmax = maximum(keys(countdict))
# 		numcounts = sum(values(countdict))
# 		counts, countfreq = formcountfreq(countdict)

# 		results = fitgene_MLE(countdict, lower, upper, MoMinitial(countdict, thresh))
# 		push!(fitres, results)
# 		println(TF)

# 		νfit, γfit, μfit, rfit, ffit = results.minimizer
	
# 		pnoexlist=taylor_expand(z->ffit*Gnoex(z, νfit, ffit*γfit, μfit, rfit), order=countmax).coeffs
# 		peak = max(maximum(pnoexlist), 0.0015/4)

# 		pnolist = taylor_expand(z->(1-ffit)*Gnoise(z, νfit, ffit*γfit, μfit, rfit), order=countmax).coeffs

# 		push!(thresholds, counts[findfirst(pnoexlist .> pnolist)])


# 		# plot(counts, countfreq, ylim=(0, 3*peak), xlim=(0, 0.75*countmax), label="Counts", xlabel="$(TF) Count\nν=$(@sprintf("%.2g", νfit)) , γ=$(@sprintf("%.1f", γfit)) , μ=$(@sprintf("%.1f", μfit)) , r=$(@sprintf("%.3g", rfit)), f=$(@sprintf("%.2g", ffit))", ylabel="Fraction of Droplets")
	
# 		# plot!(0:countmax, pnolist, label="Noise", lw=2)

# 		# plot!(0:countmax, pnoexlist, label="Expression", lw=2)
	
# 		# plot!(0:countmax, taylor_expand(z->Gseq(z, results.minimizer...), order=countmax).coeffs .+ peak*2e-2, label="Mixture", lw=2)

# 		# savefig("/Users/forrestsheldon/Dropbox/Projects/GeneBinarisation/Plots/GeneFits/$(TF)_MixtureFit_$(thresh).png")
# 	end
	
# end

# ╔═╡ 1808e0b7-43d0-410b-9f7a-755764fb1584
median([res.minimizer[1]*res.minimizer[2] for res in fitres])

# ╔═╡ 7809c70d-eefd-4c83-ab42-265d0972c5e0
begin
	histogram([res.minimizer[1]*res.minimizer[2] for res in fitres], label = false, size=(300,200), tickfontsize=12)#, xlabel="νγ", ylabel = "Number of TFs")
	# savefig("Plots/Histnugamma.png")
end

# ╔═╡ 7615e8b4-8216-465e-92e5-60ca30833a2f
begin
	scatter([res.minimizer[1] for res in fitres], [res.minimizer[2] for res in fitres], legend=false, tickfontsize=24, markerstrokewidth=0, markersize=8, size=(900, 500))
	# savefig("Plots/fitparam.png")
end

# ╔═╡ a5e0deb8-6e89-4a8e-8044-f851e181c3e0
begin
	scatter([res.minimizer[1] for res in fitres], [res.minimizer[2] for res in fitres], legend=false, tickfontsize=24, markerstrokewidth=0, markersize=8, size=(900, 500), xscale=:log10, yscale=:log10)
	x1 = 10 .^collect(-4:0.1:-2)
	plot!(x1,  0.27 ./ x1, lw=4)
	# savefig("Plots/fitparamcorr.png")
end

# ╔═╡ dc6c9bbb-9766-4df9-aeb2-0f640e14e555
begin
	histogram(thresholds, label = false, size=(600,400), xaxis=:identity, tickfontsize=24)#, xlabel="Threshold", ylabel="Number of TFs")
	# savefig("Plots/Thresholds.png")
end

# ╔═╡ a23becdf-a966-4f5e-9fbd-129ce1bd799b
thresholds

# ╔═╡ Cell order:
# ╠═dd3be7a0-cd7e-11ed-0c04-77087f8f9e51
# ╠═be5dc067-4398-4511-a741-7c8115656aab
# ╠═9edb8c81-c6b8-4422-8ab6-ab58a3810e0e
# ╠═a55ade2e-de74-40f7-9b7a-2cd50a94c16a
# ╟─275c2e81-026c-4a7c-9a4d-303ebbb08f56
# ╠═f7a000a2-ad36-4506-94e9-28cc451ad463
# ╠═c4f26fa4-e9c8-4dcd-aa0f-88ecae320885
# ╠═d2827474-d4a4-425b-bdc8-6ded3177449b
# ╠═09257638-4766-44a3-80d1-7d5332e62afe
# ╠═6ce3fb4a-bfde-4796-8e54-631e7d825e5f
# ╠═46f9cce2-9b21-49c1-92df-c964c1181c81
# ╠═414c2324-2556-4eca-ab5a-cd94ab09e8fb
# ╠═f117944e-b6a5-40a9-ba8b-28c0578afbeb
# ╠═39fac9de-58d8-47c3-8025-f237be4f8e56
# ╠═3cb2d27e-37ae-43ea-ac38-842f5b20eeb8
# ╠═3a33e1f4-ecde-4be6-9eb3-8f83ad82321c
# ╠═bbd1ddee-8717-4ce5-b855-3f168a33193c
# ╠═1808e0b7-43d0-410b-9f7a-755764fb1584
# ╠═7809c70d-eefd-4c83-ab42-265d0972c5e0
# ╠═7615e8b4-8216-465e-92e5-60ca30833a2f
# ╠═a5e0deb8-6e89-4a8e-8044-f851e181c3e0
# ╠═dc6c9bbb-9766-4df9-aeb2-0f640e14e555
# ╠═a23becdf-a966-4f5e-9fbd-129ce1bd799b
