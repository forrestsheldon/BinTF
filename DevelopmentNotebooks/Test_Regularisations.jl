### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 7878da9c-fe22-11ed-35ff-afa72bda41b4
begin
	using Pkg
	Pkg.activate("../Project.toml")
end

# ╔═╡ 3be08f88-58e2-47e1-80d4-70bf3e31f824
begin
	using DataFrames, CSV, JSON
	using TaylorSeries, Optim, StatsBase
	using Plots
	
	include("../FitUtilities.jl")
end

# ╔═╡ e699d95c-ad8f-4c06-840c-a70a8861b47a
begin
	#######################################################
	# Set paths and parameters for screen
	#######################################################
	replist = readlines("../replicates.txt")
	condlist = readlines("../condition.txt")
	datapathprefix = "../Data"
	
	rep = replist[1]
	cond = condlist[1]
	contdist = "Full"
	
	datapath = joinpath(datapathprefix, "$(cond)_"*join(replist))
	outputpath = joinpath(datapath, "GeneFitData", "allfits")
	
	countpath = joinpath(datapath, "$(cond)_$(rep)_CellCounts.tsv")
	parampath = joinpath(datapath, "GeneFitData", "contaminantfits", "$(cond)_$(rep)_ContaminantParameters_$(contdist).txt")
	noBCcells = parse(Int, readline(joinpath(datapath, "$(cond)_$(rep)_noBCcells.txt")))

	contfitpath = joinpath(datapath, "GeneFitData", "contaminantfits","$(cond)_$(rep)_Parameters_$(contdist).tsv")
	selectedpath = joinpath(datapath, "GeneFitData", "contaminantfits", "$(cond)_$(rep)_SelectedFits_$(contdist).json")
end

# ╔═╡ 3894a63c-a2c6-4fbd-b6ce-ed959e9266c9
begin
	allparamdf = load_tsv(contfitpath)
	selected = JSON.parse(readline(selectedpath))
	paramdf = allparamdf[:,selected]
	
	ρμvec = paramdf[3,:]
	αvec = paramdf[4,:]
	# fvec = paramdf[5,:]
	α_mean, α_var = mean(αvec), var(αvec)
	α_med = median(αvec)
	kα, θα = α_mean^2/α_var, α_var/α_mean
end

# ╔═╡ f6f1cd10-aa02-4b6f-895c-32ec22fa0803
begin
	#######################################################
	# Data Processing Loops
	#######################################################
	# Load the Noise Constant set by the last fit
	contparam = [parse(Float64, s) for s in readlines(parampath)]

	# Load Data
	println("Loading data for $(cond) $(rep)")
	allcountsdf = load_tsv(countpath)
	BClist = names(allcountsdf)[2:end]
end

# ╔═╡ 839fd959-027b-4d73-976f-9bd95fe33cb8
begin
	# weighting by sample size makes almost no difference in the final
	# estimates of the paramters
	# Just the median will produce a good estimate
	numcounts_selected = sum.(eachcol(allcountsdf[:, selected]))
	
	paramvec = [paramdf[4, plasmid] for plasmid in selected]
	
	medcount = median(paramvec, FrequencyWeights(numcounts_selected))
	med = median(paramvec)
	
end

# ╔═╡ 63da7372-6406-44d2-bfa9-fde2de738340
begin
	scatter(1 ./numcounts_selected, paramvec)
	plot!([0, maximum(1 ./numcounts_selected)], [medcount, medcount])
	plot!([0, maximum(1 ./numcounts_selected)], [med, med])
end

# ╔═╡ da17eb63-fad8-46e5-8075-0118fd373c8c
begin
	αmed = median([paramdf[4, plasmid] for plasmid in selected])
	νmed = median([paramdf[1, plasmid] for plasmid in selected])
	γmed = median([paramdf[2, plasmid] for plasmid in selected])
	νγmed = median([paramdf[1, plasmid]*paramdf[2, plasmid] for plasmid in selected])
	(νmed*γmed - νγmed)/νγmed
end

# ╔═╡ fdd77086-3c6b-48f3-9d6d-715aa8a8215a
begin
	function fμMOM(Savg, νγ)
		Savg/(1+νγ)
	end
	
	function μMOM(Savg, S2avg, ν, γ, α)
		fμ = fμMOM(Savg, ν*γ)
		((S2avg - Savg^2)/fμ + fμ - (1+ν*γ))/((1+α)*(1+ν^2*γ))
	end

	fμvec = [fμMOM(mean(allcountsdf[!, plasmid]), νγmed) for plasmid in selected]
	μvec = [μMOM(mean(allcountsdf[!, plasmid]), mean(allcountsdf[!, plasmid].^2), νγmed/γmed, γmed, αmed) for plasmid in selected]
	fvec = fμvec ./ μvec
end

# ╔═╡ 95e34375-4c15-48fe-839a-fd0717f0802a
begin
	# logPα(α) = (kα-1)*log(α) - α/θα
	logPα(α) = -α/(100*α_med)

	f_mean = mean(fvec)
	αf = 1.5
	βf = αf*(1-f_mean)/f_mean

	logPf(f) = (αf-1)*log(f) + (βf-1)*log(1-f)
end

# ╔═╡ befb945e-0534-4eea-9e20-fe39a6706a6a
scatter(μvec, [paramdf[3, plasmid] for plasmid in selected])

# ╔═╡ 5c0e1a79-f404-4810-96be-3507adede964
scatter(fvec, [paramdf[5, plasmid] for plasmid in selected])

# ╔═╡ 8c7f3c76-d946-4d6c-98db-7dcaebace718
(fvec .- [paramdf[5, plasmid] for plasmid in selected])./fvec

# ╔═╡ f0661533-954f-44aa-b50c-5b3947589644
begin
	BC = BClist[10]
	BC = selected[10]
	BCcountdict = countmap(allcountsdf[!, BC])
    BCcountdict[0] += noBCcells

	countvec = collect(keys(BCcountdict))
	sort!(countvec)
	numcountvec = [BCcountdict[c] for c in countvec]
	fw = FrequencyWeights(numcountvec)
	qthresh = max(1, quantile(countvec, fw, 1-f_mean))
end

# ╔═╡ 6b45c684-8798-4441-b471-5ddfc24dabf6
plot(countvec, numcountvec, ylim = (0, 10))

# ╔═╡ bfd165fc-bed5-4b9d-92d8-4e652ce390bf
begin
    νi, γi, ρμi, αi, fi = MoMFull(BCcountdict, qthresh)
	γ_med, γ_MAD, νγ_med, νγ_MAD = contparam
    
end

# ╔═╡ bcb24809-8e51-444d-920a-4c165fd33b44
begin
	maxcount = maximum(keys(BCcountdict))
	sumcounts = sum(nc for (c, nc) in BCcountdict)

	# function objective(param)
 #        ν, γ, ρμ, α, f = param
 #        ν = max(1e-20, ν)
 #        γ = max(1e-20, γ)
	# 	ρμ = max(1e-20, ρμ)
	# 	α = max(1e-20, α)
	# 	f = min(1, max(0, f))
	# 	gex, gcont, gseq = genfuncs(ν, γ, ρμ, α, f)
    
 #        plist = taylor_expand(gseq, order = maxcount).coeffs
	# 	plist[plist .< 0] .= eps(typeof(plist[1]))
        
 #        sum(-nc*log(plist[c+1]) for (c, nc) in BCcountdict) / sumcounts + (ν*γ - νγ_med)^2/(2*(1.48*νγ_MAD)^2) + (γ - γ_med)^2/(2*(1.48*γ_MAD)^2)

 #    end

	# objective(initial)
end

# ╔═╡ 3f4c400d-190b-45c7-a419-4710e683c8c9
begin
	# fitting μ and f with fixed ν, γ, α
	function objective(param)
		ρμ, f = param
		gex, gcont, gseq = genfuncs(νγ_med/γ_med, γ_med, ρμ, α_med, f)
		plist = taylor_expand(gseq, order = maxcount).coeffs
		plist[plist .< 0] .= eps(typeof(plist[1]))
	
		sum(-nc*log(plist[c+1]) for (c, nc) in BCcountdict) / sumcounts
	end
	lower = [0., 0.]
	upper = [1e3, 1.]
	initial = [ρμi, fi]
	
	objective([ρμi, fi])
end

# ╔═╡ da1aaf34-d9f4-4548-ac9c-f948a9329713
# begin
#   # Fitting μ, f and α with an exponential regularisation for α
# 	function objective(param)
# 		ρμ, α, f= param
# 		gex, gcont, gseq = genfuncs(νγ_med/γ_med, γ_med, ρμ, α, f)
# 		plist = taylor_expand(gseq, order = maxcount).coeffs
# 		plist[plist .< 0] .= eps(typeof(plist[1]))
	
# 		sum(-nc*log(plist[c+1]) for (c, nc) in BCcountdict) / sumcounts - logPα(α)
# 	end
# 	lower = [0., 0., 0.]
# 	upper = [1e3, 10., 1.]
# 	initial = [ρμi, α_med, fi]
# 	objective([ρμi, α_med, fi])
# end

# ╔═╡ a8e79472-807c-4959-a9f1-593c7926c6d0
fitresults = Optim.optimize(objective, lower, upper, initial)

# ╔═╡ 336442e1-a717-4e8b-9da6-e9aec532dc47
begin
	# begin
	# 	fitresults = fitgene_MLE_reg(celldict, regparam, lower, upper, initial)
	ρμ, f = fitresults.minimizer
	ν = νγ_med/γ_med
	γ = γ_med
	α = α_med
	# end
end

# ╔═╡ 8b3b50f5-57b3-4a16-9bc3-a55eeefc3257
begin
	# mixture fits
	cellcountvec, cellcountfreq = formcountfreq(BCcountdict)

	
	gex, gcont, gseq = genfuncs(ν, γ, ρμ, α, f)
	gcontex(z) = gex(z)*gcont(z)
	pcoex = taylor_expand(gcontex, order = cellcountvec[end]).coeffs
	pseq = taylor_expand(gseq, order = cellcountvec[end]).coeffs
	pco = taylor_expand(gcont, order = cellcountvec[end]).coeffs
	
	threshold = cellcountvec[findfirst(f.*pcoex[1:length(cellcountvec)] .> (1-f).*pco[1:length(cellcountvec)])]
	
end

# ╔═╡ 6f2eacee-50ad-44ec-bae2-e5d361b73b3f
begin
	plot(cellcountvec, cellcountfreq, label = "Counts",ylim=(0, 0.005))
	plot!(cellcountvec, pseq, label = "Pseq")
	plot!(cellcountvec, f.*pcoex, label = "Pex")
	plot!(cellcountvec, (1-f).*pco, label="Pno")
end

# ╔═╡ f207c6ed-6b7f-4d57-b361-83bdcee0984f
begin
	fμ = fμMOM(mean(allcountsdf[!, BC]), νγmed)
	μ = μMOM(mean(allcountsdf[!, BC]), mean(allcountsdf[!, BC].^2), νγmed/γmed, γmed, αmed)
	gex2, gcont2, gseq2 = genfuncs(νγmed/γmed, γmed, μ, αmed, fμ/μ)

	gcontex2(z) = gex2(z)*gcont2(z)
	pcoex2 = taylor_expand(gcontex2, order = cellcountvec[end]).coeffs
	pseq2 = taylor_expand(gseq2, order = cellcountvec[end]).coeffs
	pco2 = taylor_expand(gcont2, order = cellcountvec[end]).coeffs
	
	threshold2 = cellcountvec[findfirst(fμ/μ.*pcoex2[1:length(cellcountvec)] .> (1-fμ/μ).*pco2[1:length(cellcountvec)])]
	
end

# ╔═╡ a41a8618-1358-4bca-abd9-6ed543204984
begin
	plot(cellcountvec, cellcountfreq, label = "Counts",ylim=(0, 0.005))
	plot!(cellcountvec, pseq2, label = "Pseq")
	plot!(cellcountvec, fμ/μ.*pcoex2, label = "Pex")
	plot!(cellcountvec, (1-fμ/μ).*pco2, label="Pno")
end

# ╔═╡ ccddf90a-2dca-4806-ad42-b4e25cbfbfdd
abs(fμ/μ - f)/(fμ/μ)

# ╔═╡ Cell order:
# ╠═7878da9c-fe22-11ed-35ff-afa72bda41b4
# ╠═3be08f88-58e2-47e1-80d4-70bf3e31f824
# ╠═e699d95c-ad8f-4c06-840c-a70a8861b47a
# ╠═3894a63c-a2c6-4fbd-b6ce-ed959e9266c9
# ╠═95e34375-4c15-48fe-839a-fd0717f0802a
# ╠═f6f1cd10-aa02-4b6f-895c-32ec22fa0803
# ╠═839fd959-027b-4d73-976f-9bd95fe33cb8
# ╠═63da7372-6406-44d2-bfa9-fde2de738340
# ╠═da17eb63-fad8-46e5-8075-0118fd373c8c
# ╠═fdd77086-3c6b-48f3-9d6d-715aa8a8215a
# ╠═befb945e-0534-4eea-9e20-fe39a6706a6a
# ╠═5c0e1a79-f404-4810-96be-3507adede964
# ╠═8c7f3c76-d946-4d6c-98db-7dcaebace718
# ╠═f0661533-954f-44aa-b50c-5b3947589644
# ╠═6b45c684-8798-4441-b471-5ddfc24dabf6
# ╠═bfd165fc-bed5-4b9d-92d8-4e652ce390bf
# ╠═bcb24809-8e51-444d-920a-4c165fd33b44
# ╠═3f4c400d-190b-45c7-a419-4710e683c8c9
# ╠═da1aaf34-d9f4-4548-ac9c-f948a9329713
# ╠═a8e79472-807c-4959-a9f1-593c7926c6d0
# ╠═336442e1-a717-4e8b-9da6-e9aec532dc47
# ╠═8b3b50f5-57b3-4a16-9bc3-a55eeefc3257
# ╠═6f2eacee-50ad-44ec-bae2-e5d361b73b3f
# ╠═f207c6ed-6b7f-4d57-b361-83bdcee0984f
# ╠═a41a8618-1358-4bca-abd9-6ed543204984
# ╠═ccddf90a-2dca-4806-ad42-b4e25cbfbfdd
