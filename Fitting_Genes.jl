### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 9f56ed72-a60a-11ed-0a9c-078cd1e09182
begin
	using Pkg
	Pkg.activate("Project.toml")
end

# ╔═╡ 04c7ebc3-d911-4f6d-ac1d-c21262506972
begin
	using StatsBase
	using TaylorSeries
	using CSV, DataFrames
	using Plots, Plots.PlotMeasures, StatsPlots
	using Printf
end

# ╔═╡ 7eb6af35-b356-4b9d-a203-22913e239eab
include("fitutils.jl")

# ╔═╡ adffddb0-d063-490a-a17d-4587fa350749
md"""
### Load Data
"""

# ╔═╡ 02d9a29a-b282-48c2-a883-c06295f576b9
begin
	replist = ["R01", "R02"]
	condlist = ["PreAmp", "PostAmp"]

	replicate = replist[1]
	condition = condlist[1]
	runident = "$(condition)_$(replicate)"

	
	datapath = "./Data/$(condition)_$(join(replist))"
	celldf = DataFrame(CSV.File(joinpath(datapath, "$(runident)_CellCounts.tsv"), delim="\t"))
end

# ╔═╡ 70400124-5b17-4ad6-bddb-4a2715d0d7d2
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

# ╔═╡ 7b028b95-79ee-4133-a986-64be7579c1df
begin
	TFlist = names(celldf)[2:end]
end

# ╔═╡ 127cd7d8-241d-40ff-a9f6-5e2e2f307d56
begin
	TFkept = []

	# Filter for genes with at least 200 counts greater than 10
	for TF in TFlist
		if sum(celldf[:, TF] .> 5) > 100
			push!(TFkept, TF)
		end
	end
	println("$(length(TFkept)) of $(length(TFlist)) Plasmids kept for fitting")
end

# ╔═╡ 1044643c-884f-4017-a7c8-55a129374cb9
md"""
### Test Gene
"""

# ╔═╡ fdcea879-6ae0-4e74-b931-934e551f2fed
begin
	testTF = TFkept[3]
	testcountdict = formcountdict(testTF, celldf, 0)
	testmax = maximum(keys(testcountdict))
	testnumcounts = sum(values(testcountdict))
	testcounts, testcountfreq = formcountfreq(testcountdict)
end

# ╔═╡ 725a2c87-b9c3-418c-8fc1-48c448713400
plot(testcounts, testcountfreq, ylim=(0, 0.003), label=false, xlabel="Counts", ylabel="Frequency")

# ╔═╡ 31d069e8-3838-4686-9f05-8ae423f70c47
begin
	testlower = [0., 0., 0., 0., 0.]
	testupper = [1., 1e4, 1e4, 10, 1.]
	testthresh = 10
	
	testres = fitgene_MLE(testcountdict, testlower, testupper, MoMinitial(testcountdict, testthresh))
end

# ╔═╡ 426cd89c-25de-4c59-98cd-a0a976adf186
begin
	νfit, γfit, μfit, rfit, ffit = testres.minimizer
	
	pnoexlist=taylor_expand(z->ffit*Gnoex(z, νfit, ffit*γfit, μfit, rfit), order=testmax).coeffs
	peak = max(maximum(pnoexlist), 0.0015/4)

	plot(testcounts, testcountfreq, ylim=(0, 0.003), xlim=(0, 0.75*testmax), label="Counts", xlabel="$(testTF) Count\nν=$(@sprintf("%.3g", νfit)) , γ=$(@sprintf("%.1f", γfit)) , μ=$(@sprintf("%.1f", μfit)) , r=$(@sprintf("%.3g", rfit)), f=$(@sprintf("%.3g", ffit))", ylabel="Fraction of Droplets")
	
	plot!(0:testmax, taylor_expand(z->(1-ffit)*Gnoise(z, νfit, ffit*γfit, μfit, rfit), order=testmax).coeffs, label="Noise", lw=2)

	plot!(0:testmax, pnoexlist, label="Expression", lw=2)
	
	plot!(0:testmax, taylor_expand(z->Gseq(z, testres.minimizer...), order=testmax).coeffs .+ peak*2e-2, label="Mixture", lw=2)
	plot!(0:testmax, taylor_expand(z->Gseq(z, MoMinitial(testcountdict, testthresh)...), order=testmax).coeffs, label="Initial", lw=2)
end

# ╔═╡ af7132a5-7f3a-4c84-a093-c403387016d7
# begin
# # Make similar plots without labels for placing in panel figures
# 	νfit, γfit, μfit, rfit, ffit = testres.minimizer
	
# 	pnoexlist=taylor_expand(z->ffit*Gnoex(z, νfit, ffit*γfit, μfit, rfit), order=testmax).coeffs
# 	peak = max(maximum(pnoexlist), 0.0015/4)

# 	plot(testcounts, testcountfreq, ylim=(0, 0.007), xlim=(0, 102), label="Counts", tickfontsize=12, legendfontsize=12, size=(500, 300), xticks=[0, 25, 50, 75, 100], yticks=[0., 0.002, 0.004, 0.006])#, xlabel="$(testTF) Count\nν=$(@sprintf("%.3g", νfit)) , γ=$(@sprintf("%.1f", γfit)) , μ=$(@sprintf("%.1f", μfit)) , r=$(@sprintf("%.3g", rfit)), f=$(@sprintf("%.3g", ffit))", ylabel="Fraction of Droplets")
	
# 	plot!(0:testmax, taylor_expand(z->(1-ffit)*Gnoise(z, νfit, ffit*γfit, μfit, rfit), order=testmax).coeffs, label="Noise", lw=2)

# 	plot!(0:testmax, pnoexlist, label="Expression", lw=2)
	
# 	plot!(0:testmax, taylor_expand(z->Gseq(z, testres.minimizer...), order=testmax).coeffs, label="Mixture", lw=2)
# 	plot!(0:testmax, taylor_expand(z->Gseq(z, MoMinitial(testcountdict, testthresh)...), order=testmax).coeffs, label="Initial", lw=2, legend=false)

# 	savefig("Plots/GeneFit_$(runident)_Plasmid_$(testTF)_components.png")
# end

# ╔═╡ cbd23588-9a0f-4b99-ad6e-0c81506e25d6
md"""
### Fit All Genes
"""

# ╔═╡ ff816c47-ec1e-4969-ba98-585232eeca05
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


# 		plot(counts, countfreq, ylim=(0, 3*peak), xlim=(0, 0.75*countmax), label="Counts", xlabel="$(TF) Count\nν=$(@sprintf("%.2g", νfit)) , γ=$(@sprintf("%.1f", γfit)) , μ=$(@sprintf("%.1f", μfit)) , r=$(@sprintf("%.3g", rfit)), f=$(@sprintf("%.2g", ffit))", ylabel="Fraction of Droplets")
	
# 		plot!(0:countmax, pnolist, label="Noise", lw=2)

# 		plot!(0:countmax, pnoexlist, label="Expression", lw=2)
	
# 		plot!(0:countmax, taylor_expand(z->Gseq(z, results.minimizer...), order=countmax).coeffs .+ peak*2e-2, label="Mixture", lw=2)

# 		savefig("/Users/forrestsheldon/Dropbox/Projects/GeneBinarisation/Plots/GeneFits/$(TF)_MixtureFit_$(thresh).png")
# 	end
	
# end

# ╔═╡ 7864f60b-c61f-4e03-9952-15dd110e3304
md"""
### Examining Results of the fit
"""

# ╔═╡ 6fa4ddf4-77a1-496d-b20e-de76694e8435
begin
	histogram([res.minimizer[1] for res in fitres], nbins=8, label=false, size=(300,200),xticks = [0, 0.002, 0.004], yticks = [0, 5, 10],tickfontsize=12)#, xlabel="ν", ylabel = "Number of TFs")
	# savefig("Plots/Histnu.png")
end

# ╔═╡ d805a9ca-d12e-485f-9206-3996605002d5
begin
	histogram([res.minimizer[2] for res in fitres], label = false, size=(300,200), tickfontsize=12)#, xlabel="γ", ylabel = "Number of TFs")
	# savefig("Plots/Histgamma.png")
end

# ╔═╡ 540d79c6-07d6-48f5-97e8-bb95e6e61c61
begin
	histogram([res.minimizer[1]*res.minimizer[2] for res in fitres], label = false, xlim = (0.013,  0.033), size=(300,200), tickfontsize=12)#, xlabel="νγ", ylabel = "Number of TFs")
	# savefig("Plots/Histnugamma.png")
end

# ╔═╡ 36d0fd95-ee29-4cd8-8399-d957868fe4ce
median([res.minimizer[1]*res.minimizer[2] for res in fitres])

# ╔═╡ 01c4d4d0-24ee-4cc2-a81a-65cda8583b4f
md"""
### Errors
"""

# ╔═╡ 92a2c98e-3ab6-4e74-b8d6-f0965f636a22
begin
	fpall = []
	fnall = []
	for (erridx, TF) in enumerate(TFkept)
		
		thresherr = thresholds[erridx]
		
		countdicterr = formcountdict(TF, celldf, 0)
		countmaxerr = maximum(keys(countdicterr))
		
		νerr, γerr, μerr, rerr, ferr = fitres[erridx].minimizer
	
		pnoexlisterr=taylor_expand(z->ferr*Gnoex(z, νerr, ferr*γerr, μerr, rerr), order=countmaxerr).coeffs
	
		pnolisterr = taylor_expand(z->(1-ferr)*Gnoise(z, νerr, ferr*γerr, μerr, rerr), order=countmaxerr).coeffs
	
		psumerr = pnoexlisterr .+ pnolisterr
	
		fplist = []
		fnlist=[]
	
		for (c, nc) in countdicterr
			if c < thresherr
				for i in 1:nc
					push!(fnlist, pnoexlisterr[c+1]/psumerr[c+1])
				end
			else
				for i in 1:nc
					push!(fplist, pnolisterr[c+1]/psumerr[c+1])
				end
			end
		end

		push!(fpall, fplist)
		push!(fnall, fnlist)
	end
end

# ╔═╡ b380cb6d-7532-4159-aecd-a38a457d6cf1
begin
	fpflatten = round.(vcat(fpall...), digits=2)
	histogram(fpflatten, nbins=50, legend=false)
end

# ╔═╡ 7e50d391-3fa6-448b-8547-939022c01168
histogram(round.(vcat(fnall...), digits=2), nbins=50, legend=false)

# ╔═╡ 446147dd-3af5-459d-ad58-88345fb16dc3
md"""
### False Positive rates with Threshold 1
"""

# ╔═╡ 0fb0dd35-06c0-43c4-a321-60d0fe39e4e1
begin

	fprates = []
	for (TFidx, TF) in enumerate(TFkept)
		TF = TFkept[TFidx]
		thresherr = thresholds[TFidx]
		
		countdicterr = formcountdict(TF, celldf, 0)
		countmaxerr = maximum(keys(countdicterr))
	
		fpcounts = 0
		tpcounts = 0
	
		for (c, nc) in countdicterr
			if c >= 1 && c < thresherr
				fpcounts += nc
			elseif c >= thresherr
				tpcounts += nc
			end
		end
	
		push!(fprates, fpcounts/(fpcounts + tpcounts))
	end
	
end

# ╔═╡ 97dd83cf-c914-4334-994a-bf28f352af9c
begin
	histogram(fprates, xlim = (0, 1.03), legend=false, tickfontsize=25, size=(600, 400), xticks=[0, 0.2, 0.4, 0.6, 0.8, 1.0])
	# savefig("Plots/FalsePositive.png")
end

# ╔═╡ f8d5e9c8-7f2a-486d-b6d7-dfed2450c96b
begin
	histogram(thresholds, label = false, size=(600,400), xaxis=:identity, tickfontsize=24, yticks = [0, 10, 20])#, xlabel="Threshold", ylabel="Number of TFs")
	# savefig("Plots/Thresholds.png")
end

# ╔═╡ Cell order:
# ╠═9f56ed72-a60a-11ed-0a9c-078cd1e09182
# ╠═04c7ebc3-d911-4f6d-ac1d-c21262506972
# ╠═7eb6af35-b356-4b9d-a203-22913e239eab
# ╟─adffddb0-d063-490a-a17d-4587fa350749
# ╠═02d9a29a-b282-48c2-a883-c06295f576b9
# ╠═70400124-5b17-4ad6-bddb-4a2715d0d7d2
# ╠═7b028b95-79ee-4133-a986-64be7579c1df
# ╠═127cd7d8-241d-40ff-a9f6-5e2e2f307d56
# ╟─1044643c-884f-4017-a7c8-55a129374cb9
# ╠═fdcea879-6ae0-4e74-b931-934e551f2fed
# ╠═725a2c87-b9c3-418c-8fc1-48c448713400
# ╠═31d069e8-3838-4686-9f05-8ae423f70c47
# ╠═426cd89c-25de-4c59-98cd-a0a976adf186
# ╠═af7132a5-7f3a-4c84-a093-c403387016d7
# ╟─cbd23588-9a0f-4b99-ad6e-0c81506e25d6
# ╠═ff816c47-ec1e-4969-ba98-585232eeca05
# ╟─7864f60b-c61f-4e03-9952-15dd110e3304
# ╠═6fa4ddf4-77a1-496d-b20e-de76694e8435
# ╠═d805a9ca-d12e-485f-9206-3996605002d5
# ╠═540d79c6-07d6-48f5-97e8-bb95e6e61c61
# ╠═36d0fd95-ee29-4cd8-8399-d957868fe4ce
# ╟─01c4d4d0-24ee-4cc2-a81a-65cda8583b4f
# ╠═92a2c98e-3ab6-4e74-b8d6-f0965f636a22
# ╠═b380cb6d-7532-4159-aecd-a38a457d6cf1
# ╠═7e50d391-3fa6-448b-8547-939022c01168
# ╟─446147dd-3af5-459d-ad58-88345fb16dc3
# ╠═0fb0dd35-06c0-43c4-a321-60d0fe39e4e1
# ╠═97dd83cf-c914-4334-994a-bf28f352af9c
# ╠═f8d5e9c8-7f2a-486d-b6d7-dfed2450c96b
