### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 9f56ed72-a60a-11ed-0a9c-078cd1e09182
begin
	using Pkg
	Pkg.activate("../Project.toml")
end

# ╔═╡ 04c7ebc3-d911-4f6d-ac1d-c21262506972
begin
	using StatsBase
	using TaylorSeries
	using DataFrames, CSV, JSON
	using Plots, Plots.PlotMeasures, StatsPlots
	using Printf
	using Colors
end

# ╔═╡ 7eb6af35-b356-4b9d-a203-22913e239eab
include("../FitUtilities.jl")

# ╔═╡ adffddb0-d063-490a-a17d-4587fa350749
md"""
### Load Data
"""

# ╔═╡ ab46ed63-338f-4297-94d1-0bda0b9a9704
begin
	#######################################################
	# Set paths and parameters for screen
	#######################################################
	replist = ["R01", "R02"]
	condlist = ["withoutAmp", "withAmp"]
	contdist = "Full"
	
	datapathprefix = "../Data"
	
	rep = replist[1]
	cond = condlist[1]
	
	datapath(cond, rep) = joinpath(datapathprefix, "$(cond)_"*join(replist))
	outputpath(cond, rep) = joinpath(datapath(cond, rep), "GeneFitData", "allfits")
	
	countpath(cond, rep) = joinpath(datapath(cond, rep), "$(cond)_$(rep)_CellCounts.tsv")
	parampath(cond, rep) = joinpath(datapath(cond, rep), "GeneFitData", "contaminantfits", "$(cond)_$(rep)_ContaminantParameters_$(contdist).txt")
	noBCcells(cond, rep) = parse(Int, readline(joinpath(datapath(cond, rep), "$(cond)_$(rep)_noBCcells.txt")))
	ddPCR(cond, rep) = parse(Float64, readline(joinpath(datapath(cond, rep), "$(rep)_ddPCR.txt")))

	contfitpath(cond, rep) = joinpath(datapath(cond, rep), "GeneFitData", "contaminantfits","$(cond)_$(rep)_Parameters_$(contdist).tsv")
	selectedpath(cond, rep) = joinpath(datapath(cond, rep), "GeneFitData", "contaminantfits", "$(cond)_$(rep)_SelectedFits_$(contdist).json")

	labelpath(cond, rep) = joinpath(datapath(cond, rep), "$(cond)_$(rep)_Labeled_$(contdist).tsv")
	finalfitpath(cond, rep) = joinpath(datapath(cond, rep), "GeneFitData", "allfits","$(cond)_$(rep)_Parameters_$(contdist).tsv")
	
end

# ╔═╡ 33361e24-b008-4996-b36e-90ffbd1c9874
md"""
### Contaminant Fit
"""

# ╔═╡ 1c1e08d4-dc8d-46bb-a275-b8146b1f75c0
begin
	contparamdfwo1 = load_tsv(contfitpath("withoutAmp", "R01"))
	selectedwo1 = JSON.parse(readline(selectedpath("withoutAmp", "R01")))
	contparamdfwo1 = contparamdfwo1[:,selectedwo1]

	contparamdfwo2 = load_tsv(contfitpath("withoutAmp", "R02"))
	selectedwo2 = JSON.parse(readline(selectedpath("withoutAmp", "R02")))
	contparamdfwo2 = contparamdfwo2[:,selectedwo2]

	contparamdfw1 = load_tsv(contfitpath("withAmp", "R01"))
	selectedw1 = JSON.parse(readline(selectedpath("withAmp", "R01")))
	contparamdfw1 = contparamdfw1[:,selectedw1]

	contparamdfw2 = load_tsv(contfitpath("withAmp", "R02"))
	selectedw2 = JSON.parse(readline(selectedpath("withAmp", "R02")))
	contparamdfw2 = contparamdfw2[:,selectedw2]
end

# ╔═╡ 87d16f2c-b2a9-44e6-a4ec-88eda65a8739
begin
	νwo1 = Vector(contparamdfwo1[1,:])
	γwo1 = Vector(contparamdfwo1[2,:])
	νγmedwo1 = median(νwo1.*γwo1)

	νwo2 = Vector(contparamdfwo2[1,:])
	γwo2 = Vector(contparamdfwo2[2,:])
	νγmedwo2 = median(νwo2.*γwo2)
	
	νw1 = Vector(contparamdfw1[1,:])
	γw1 = Vector(contparamdfw1[2,:])
	νγmedw1 = median(νw1.*γw1)

	νw2 = Vector(contparamdfw2[1,:])
	γw2 = Vector(contparamdfw2[2,:])
	νγmedw2 = median(νw2.*γw2)
end

# ╔═╡ dc7489d6-35f5-46f5-a978-e16abbf8072d
begin
	colorpalette = theme_palette(:auto)
	bluegrad = cgrad([:white, colorpalette[1], :black])
	myblues = cgrad([get(bluegrad, 0.1), colorpalette[1], get(bluegrad, 0.8)])
end

# ╔═╡ 95ee2eb2-a3f1-4e0a-99d9-3913d6334e71
begin
	orangegrad = cgrad([:white, colorpalette[2], :black])
	myoranges = cgrad([get(orangegrad, 0.1), colorpalette[2], get(orangegrad, 0.8)])
end

# ╔═╡ 29302f0f-4ccb-4905-a8d4-fc79934809df
begin
	scatter(νwo1, γwo1, markerstrokecolor=:auto, label="Replicate 1")
	linebounds = 10^-4:10^-5:0.0047
	plot!(linebounds, νγmedwo1./linebounds, color=1, lw=2, label=false)

	scatter!(νwo2, γwo2, color=2, markerstrokecolor=:auto, label="Replicate 2", tickfontsize=16, xticks=[0, 0.002, 0.004], yticks=[0, 100, 200], xlim=[0, 0.005], ylim = [0, 200], size=(600, 350))
	plot!(linebounds, νγmedwo2./linebounds, color=2, lw=2, label=false, legendfontsize=12)
	# savefig("../Plots/PreScatter.png")
end

# ╔═╡ b13cb5c8-e27b-47ef-b523-a95c73b4a23f
md"""
To look for evidence of convergence, we pull the count numbers for each gene
"""

# ╔═╡ 68822eff-a59b-486b-9ac3-5c9ae7af734d
begin
	wo1countdf = load_tsv(countpath("withoutAmp", "R01"))
	wo2countdf = load_tsv(countpath("withoutAmp", "R02"))

	w1countdf = load_tsv(countpath("withAmp", "R01"))
	w2countdf = load_tsv(countpath("withAmp", "R02"))
end

# ╔═╡ aed79458-b324-4662-a1c3-7368768a4047
begin
	# wo1numcounts = sum.(eachcol(wo1countdf[!,selectedwo1]))
	# wo2numcounts = sum.(eachcol(wo2countdf[!,selectedwo2]))

	# w1numcounts = sum.(eachcol(w1countdf[!,selectedw1]))
	# w2numcounts = sum.(eachcol(w2countdf[!,selectedw2]))

	wo1numcounts = sum.(c .!= 0 for c in eachcol(wo1countdf[!,selectedwo1]))
	wo2numcounts = sum.(c .!= 0 for c in eachcol(wo2countdf[!,selectedwo2]))

	w1numcounts = sum.(c .!= 0 for c in eachcol(w1countdf[!,selectedw1]))
	w2numcounts = sum.(c .!= 0 for c in eachcol(w2countdf[!,selectedw2]))

end

# ╔═╡ 972a908d-5fa8-406e-8832-4f30af6ffdb6
begin
	wo1perm = sortperm(wo1numcounts)
	wo2perm = sortperm(wo2numcounts)
	
	scatter(νwo1[wo1perm], γwo1[wo1perm], zcolor=log.(wo1numcounts[wo1perm]), colorbar=false, c=myblues, markerstrokewidth=0, xlim=(0, 0.005))

	scatter!(νwo2[wo2perm], γwo2[wo2perm], zcolor=log.(wo2numcounts[wo2perm]), colorbar=false, c=myoranges, markerstrokewidth=0, alpha=0.7)
	
	plot!(size=(800, 500), ylims=(0, 200), xlim=(0, 0.0047), legend=false, tickfontsize=16)

end

# ╔═╡ e2f55db3-c8d5-4f69-b540-59e7ca551abd
begin
	scatter(νw1, γw1, markerstrokecolor=:auto)
	
	scatter!(νw2, γw2, markerstrokecolor=:auto)
end

# ╔═╡ 653b2c74-2dca-4988-865c-02d99b1ee49f
begin
	w1perm = sortperm(w1numcounts)
	w2perm = sortperm(w2numcounts)
	
	scatter(νw1[w1perm], γw1[w1perm], zcolor=log.(w1numcounts[w1perm]), colorbar=false, c=myblues, markerstrokewidth=0)

	scatter!(νw2[w2perm], γw2[w2perm], zcolor=log.(w2numcounts[w2perm]), colorbar=false, c=myoranges, markerstrokewidth=0, alpha=0.7)
	
	# plot!(size=(800, 500), ylims=(0, 200), xlim=(0, 0.0047), legend=false, tickfontsize=16)

end

# ╔═╡ 293a0eb4-d300-4e67-8388-dd489e2dd427
md"""
Assemble these into a single figure
"""

# ╔═╡ d7d7ebda-ffe7-4a58-9362-7df922c5c2fd
begin
	scatter(νwo1, γwo1, label="Replicate 1", color=1, markerstrokewidth=0)
	wo1linebounds = [1.1*10^-4, 7*10^-3]
	plot!(wo1linebounds, νγmedwo1./wo1linebounds, color = 1, lw=2, label=:none)
	# scatter!([NaN], [NaN], markerstrokewidth=0, color=colorpalette[1])  # Add dummy series


	scatter!(νwo2, γwo2, label="Replicate 2", color=2, markerstrokewidth=0)
	wo2linebounds = [10^-4, 3*10^-3]
	plot!(wo2linebounds, νγmedwo2./wo2linebounds, color = 2, lw=2, label=:none)
	# # scatter!([NaN], [NaN], markerstrokewidth=0, color=colorpalette[2])  # Add dummy series
	
	scatter!(νw1, γw1, label=:none, color=1, markerstrokewidth=0)
	w1linebounds = [6*10^(-3), 2.3*10^-2]
	plot!(w1linebounds, νγmedw1./w1linebounds, color = 1, lw=2, label=:none)

	scatter!(νw2, γw2, markerstrokecolor=:auto, label=:none, tickfontsize=15, color=2, markerstrokewidth=0)
	w2linebounds = [6*10^(-3), 1.8*10^-2]
	plot!(w2linebounds, νγmedw2./w2linebounds, color = 2, lw=2, label=:none, size=(600, 350), legendfontsize=12)

	plot!(xscale=:log10, yscale=:log10, colorbar=:false, xlim=(8*10^-5, 3*10^-2), ylim=(3.5, 2.5*10^2))
	# savefig("../Plots/LogScatter_noshading.png")
end

# ╔═╡ f36c3230-92cd-4c9e-b01c-590bbe5776bd
begin
	scatter(νwo1, γwo1, markerstrokecolor=:auto, label="Replicate 1")
	plot!(linebounds, νγmedwo1./linebounds, color=1, lw=2, label=false)

	scatter!(νwo2, γwo2, color=2, markerstrokecolor=:auto, label="Replicate 2", tickfontsize=16, xticks=[0, 0.002, 0.004], yticks=[0, 100, 200], xlim=[0, 0.005], ylim = [0, 200], size=(600, 350))
	plot!(linebounds, νγmedwo2./linebounds, color=2, lw=2, label=false, legendfontsize=12, legend=false)

	scatter!(νwo1, γwo1, inset=(1, bbox(0.25, 0.3, 0.75, 0.7, :bottom, :left)), subplot = 2, label="Replicate 1", color=1, markerstrokewidth=0, xscale=:log10, yscale=:log10)
	plot!(wo1linebounds, νγmedwo1./wo1linebounds, subplot=2, color = 1, lw=2, label=:none)

	scatter!(νwo2, γwo2, subplot=2, label="Replicate 2", color=2, markerstrokewidth=0)
	plot!(wo2linebounds, νγmedwo2./wo2linebounds, subplot=2, color = 2, lw=2, label=:none)

	scatter!(νw1, γw1, subplot=2, label=:none, color=1, markerstrokewidth=0)
	plot!(w1linebounds, νγmedw1./w1linebounds, subplot=2, color = 1, lw=2, label=:none)

	scatter!(νw2, γw2, subplot=2, markerstrokecolor=:auto, label=:none, tickfontsize=15, color=2, markerstrokewidth=0)
	plot!(w2linebounds, νγmedw2./w2linebounds, subplot=2, color = 2, lw=2, label=:none)

	plot!(subplot=2, tickfontsize=16, legendfontsize=12, size=(600, 400))

	# savefig("../Plots/LogScatter_inset.png")
end

# ╔═╡ 10cdb8f2-95c2-429f-b7e4-a39a99f6db65
begin
	scatter(νwo1[wo1perm], γwo1[wo1perm], label=nothing, zcolor=log.(wo1numcounts[wo1perm]), color=myblues, markerstrokewidth=0)
	wo1lineboundsc = [1.1*10^-4, 4*10^-3]
	plot!(wo1lineboundsc, νγmedwo1./wo1lineboundsc, color = 1, lw=2, label=:none)
	# scatter!([NaN], [NaN], markerstrokewidth=0, color=colorpalette[1])  # Add dummy series


	scatter!(νwo2[wo2perm], γwo2[wo2perm], label=nothing, zcolor=log.(wo2numcounts[wo2perm]), color=myoranges, markerstrokewidth=0)
	wo2lineboundsc = [10^-4, 3*10^-3]
	plot!(wo2lineboundsc, νγmedwo2./wo2lineboundsc, color = 2, lw=2, label=:none)
	# # scatter!([NaN], [NaN], markerstrokewidth=0, color=colorpalette[2])  # Add dummy series
	
	scatter!(νw1[w1perm], γw1[w1perm], label=:none, zcolor=log.(w1numcounts[w1perm]), color=myblues, markerstrokewidth=0)
	w1lineboundsc = [6*10^(-3), 2*10^-2]
	plot!(w1lineboundsc, νγmedw1./w1lineboundsc, color = 1, lw=2, label=:none)

	scatter!(νw2[w2perm], γw2[w2perm], markerstrokecolor=:auto, label=:none, tickfontsize=15, zcolor=log.(w2numcounts[w2perm]), color=myoranges, markerstrokewidth=0)
	w2lineboundsc = [6*10^(-3), 1.6*10^-2]
	plot!(w2lineboundsc, νγmedw2./w2lineboundsc, color = 2, lw=2, label=:none, size=(600, 350), legendfontsize=12)

	plot!(xscale=:log10, yscale=:log10, colorbar=:false, xlim=(8*10^-5, 3*10^-2), ylim=(6, 2.5*10^2))
	# savefig("../Plots/LogScatter.png")
end

# ╔═╡ bdd3a1f1-8e69-4559-8f35-4c8182340425
begin
	
	histogram(γw1, lw=0, alpha=0.5, label="with Amp.", color=2)
	
	histogram!(γwo1, lw=0, alpha=0.6, label="w/o Amp.", color=1)
	
	plot!(lw=0, tickfontsize=16, legendfontsize=14, size=(600,400))

	# savefig("../Plots/Histograms_gamma1.png")
end

# ╔═╡ 277256d1-b2e4-4d52-8efd-b74b477a1fb4
begin
	histogram(γw2, lw=0, alpha=0.8)
	histogram!(γwo2, lw=0, alpha=0.5)
	plot!(lw=0, tickfontsize=16)
end

# ╔═╡ 8aa03f0c-a67a-42c1-a913-8a7a8f591b2f
begin
	histogram(vcat(γw1,γw2), lw=0, alpha=0.6, label="with Amp.")
	histogram!(vcat(γwo1, γwo2), lw=0, alpha=0.5, label="w/o Amp.")
	plot!(lw=0, tickfontsize=16, legendfontsize=14, size=(600,300))
	# savefig("../Plots/Histograms_gamma12.png")
end

# ╔═╡ 7f274c2e-dc27-4cb4-92dc-4a649ca3e3e3
md"""
To look for evidence of concentration, we look at the finite size scaling plots of the parameters
"""

# ╔═╡ ff34c6db-35f8-4e36-9f6c-f75227b39c5c
begin
	scatter(1 ./wo1numcounts[wo1perm], γwo1[wo1perm], markerstrokewidth=0,  xlabel="Inverse Total Counts 1/N", ylabel="γ")
	plot!(1 ./wo1numcounts[wo1perm], [median(γwo1) for _ in wo1perm])
	plot!(1 ./wo1numcounts[wo1perm], [median(γwo1, FrequencyWeights(wo1numcounts)) for _ in wo1perm])
end

# ╔═╡ ece369cf-5820-42c5-b4ea-5c7fa9f6cd01
median(γwo1, FrequencyWeights(wo1numcounts))

# ╔═╡ 7332a6d2-ed2c-4619-85fc-c1aa26f228e6
begin
	scatter(1 ./wo2numcounts[wo2perm], γwo2[wo2perm], markerstrokewidth=0,  xlabel="Inverse Total Counts 1/N", ylabel="γ")
	plot!(1 ./wo2numcounts[wo2perm], [median(γwo2) for _ in wo2perm])
	plot!(1 ./wo2numcounts[wo2perm], [median(γwo2, FrequencyWeights(wo2numcounts)) for _ in wo2perm])
	
end

# ╔═╡ 5ef3aa92-7638-464c-8714-1fe78f6aaa3c
median(γwo2, FrequencyWeights(wo2numcounts))

# ╔═╡ 01ea329b-4781-401b-aa8d-a7ee347a5b81
begin
	scatter(1 ./w1numcounts[w1perm], γw1[w1perm], markerstrokewidth=0,  xlabel="Inverse Total Counts 1/N", ylabel="γ")
	plot!(1 ./w1numcounts[w1perm], [median(γw1) for _ in w1perm])
	plot!(1 ./w1numcounts[w1perm], [median(γw1, FrequencyWeights(w1numcounts)) for _ in w1perm])
end

# ╔═╡ d2858b3d-e9b1-49d5-8310-136fb8cba12f
median(γw1, FrequencyWeights(w1numcounts))

# ╔═╡ a55d2c29-284e-418c-a683-c793d3db8444
begin
	scatter(1 ./w2numcounts[w2perm], γw2[w2perm], markerstrokewidth=0,  xlabel="Inverse Total Counts 1/N", ylabel="γ")
	plot!(1 ./w2numcounts[w2perm], [median(γw2) for _ in w2perm])
	plot!(1 ./w2numcounts[w2perm], [median(γw2, FrequencyWeights(w2numcounts)) for _ in w2perm])
end

# ╔═╡ 43c4454f-b901-4a53-a9fc-6fe0e32520db
median(γw2, FrequencyWeights(w2numcounts))

# ╔═╡ b8fd4822-e4d0-4cf8-9efa-589f9205ae28
begin
	αwo1 = Vector(contparamdfwo1[4,:])
	αwo2 = Vector(contparamdfwo2[4,:])

	αw1 = Vector(contparamdfw1[4,:])
	αw2 = Vector(contparamdfw2[4,:])
	
	histogram(αwo1, nbins = 40, lw=0, alpha = 0.9, label="Replicate 1")
	histogram!(αwo2, nbins = 20, lw = 0, alpha = 0.7, label="Replicate 2")
	plot!(tickfontsize = 16, legendfontsize=14, xlabel="α value", ylabel="Number of Barcodes", yticks = [0, 5, 10])

end

# ╔═╡ f030de7f-2af4-4dff-bca5-caff7ac3a0fc
begin
	histogram(αw1, nbins = 40, lw=0, alpha = 0.9, label="Replicate 1")
	histogram!(αw2, nbins = 20, lw = 0, alpha = 0.7, label="Replicate 2")
	plot!(tickfontsize = 16, legendfontsize=14, xlabel="α value", ylabel="Number of Barcodes", yticks = [0, 5, 10])
end

# ╔═╡ a29d9e64-6325-442e-9d9a-2e26730f818c
begin
	scatter(1 ./wo1numcounts[wo1perm], αwo1[wo1perm], markerstrokewidth=0, label="Replicate 1")#,  xlabel="Inverse Total Counts 1/N", ylabel="α")
	scatter!(1 ./wo2numcounts[wo2perm], αwo2[wo2perm], markerstrokewidth=0, alpha=0.7, label="Replicate 2")
	

	plot!(tickfontsize=16, legendfontsize=14, xticks=[0.0004,0.0012, 0.002], labelfontsize=16, xlabel="1/Sum of Counts", ylabel="α", size=(600, 400))

	# savefig("../Plots/alphaplot_withoutAmp.png")
	# scatter!(1 ./wo1numcounts[wo1perm], αwo2[wo1perm], markerstrokewidth=0)
	# plot!(1 ./wo1numcounts[wo1perm], [median(αwo1) for _ in wo1perm])
	# plot!(1 ./wo1numcounts[wo1perm], [median(αwo1, FrequencyWeights(wo1numcounts)) for _ in wo1perm])
end

# ╔═╡ 3ca5df27-ae5e-423d-9228-b67389132215
begin
	scatter(1 ./w1numcounts[w1perm], αw1[w1perm], markerstrokewidth=0, label="Replicate 1")#,  xlabel="Inverse Total Counts 1/N", ylabel="α")
	scatter!(1 ./w2numcounts[w2perm], αw2[w2perm], markerstrokewidth=0, alpha=0.7, label="Replicate 2")
	

	plot!(tickfontsize=16, legendfontsize=14, xticks=[0.0002, 0.0004, 0.0006], xlim=(0.000125, 0.00063), labelfontsize=16, xlabel="1/Sum of Counts", ylabel="α", ylim=(0, 7.1))

	# savefig("../Plots/alphaplot_withAmp.png")
	# scatter!(1 ./wo1numcounts[wo1perm], αwo2[wo1perm], markerstrokewidth=0)
	# plot!(1 ./wo1numcounts[wo1perm], [median(αwo1) for _ in wo1perm])
	# plot!(1 ./wo1numcounts[wo1perm], [median(αwo1, FrequencyWeights(wo1numcounts)) for _ in wo1perm])
end

# ╔═╡ c3828bf3-8401-477e-899a-b6ac375552ff
md"""
### Labeling Histograms
"""

# ╔═╡ ab30e5f7-a0ea-4e5e-8124-7a1052eef7e5
begin
	wo1labeldf = load_tsv(labelpath("withoutAmp", "R01"))
	wo2labeldf = load_tsv(labelpath("withoutAmp", "R02"))
	
	w1labeldf = load_tsv(labelpath("withAmp", "R01"))
	w2labeldf = load_tsv(labelpath("withAmp", "R02"))
end

# ╔═╡ 523a26af-9b92-45a3-9f2f-8b6f50b5a790
begin
	wo1noBC = noBCcells("withoutAmp", "R01")
	wo2noBC = noBCcells("withoutAmp", "R02")

	w1noBC = noBCcells("withAmp", "R01")
	w2noBC = noBCcells("withAmp", "R02")
end

# ╔═╡ 8a6b0fe0-9be0-4b3d-80c7-ac4bae5ac822
begin
	wo1ddPCR = ddPCR("withoutAmp", "R01")
	wo2ddPCR = ddPCR("withoutAmp", "R02")

	w1ddPCR = ddPCR("withAmp", "R01")
	w2ddPCR = ddPCR("withAmp", "R02")
end

# ╔═╡ c6757638-8730-457e-9dd5-0c12c78ba4cf
begin
	wo1uniqueinit = vcat([sum(val != 0 for val in row) for row in eachrow(wo1countdf[:, 2:end])], zeros(Int, wo1noBC))
	wo2uniqueinit = vcat([sum(val != 0 for val in row) for row in eachrow(wo2countdf[:, 2:end])], zeros(Int, wo2noBC))
	
	w1uniqueinit = vcat([sum(val != 0 for val in row) for row in eachrow(w1countdf[:, 2:end])], zeros(Int, w1noBC))
	w2uniqueinit = vcat([sum(val != 0 for val in row) for row in eachrow(w2countdf[:, 2:end])], zeros(Int, w2noBC))
end

# ╔═╡ c484975f-8bc9-413b-b28a-35947491cf60
begin
	wo1minit = mean(wo1uniqueinit)
	wo1hist = fit(Histogram, wo1uniqueinit, nbins=80)
	bar(wo1hist.edges[1][1:end-1], wo1hist.weights, legend=false, linecolor=:match, xlim=(-0.7, wo1hist.edges[1][end]+2), ylim=(0, 1.05*maximum(wo1hist.weights)), size=(600, 200), tickfontsize=16, xticks=[0, 20, 40, 60], yticks=[0, 800, 1600])
	vline!([wo1minit], line=(2, 4))
	vline!([wo1ddPCR], line=(:black, 4))

	# Fill the region between a and b
	wo1fill_between_x = [wo1ddPCR, wo1minit, wo1minit, wo1ddPCR]
	wo1fill_between_y = [0, 0, 2*maximum(wo1hist.weights), 2*maximum(wo1hist.weights)]
	plot!(wo1fill_between_x, wo1fill_between_y, fill=(0, 0.4, 1))
	# savefig("../Plots/CopyHist_PreLabel_withoutAmp_R01.png")
end

# ╔═╡ d09e3965-d0fc-499e-b303-804ec35e7023
begin
	wo2minit = mean(wo2uniqueinit)
	wo2hist = fit(Histogram, wo2uniqueinit, nbins=80)
	bar(wo2hist.edges[1][1:end-1], wo2hist.weights, legend=false, linecolor=:match, xlim=(-0.7, wo2hist.edges[1][end]+2), ylim=(0, 1.05*maximum(wo2hist.weights)), size=(600, 200), tickfontsize=16, xticks=[0, 20, 40, 60], yticks=[0, 800, 1600])
	vline!([wo2minit], line=(2, 4))
	vline!([wo2ddPCR], line=(:black, 4))

	# Fill the region between a and b
	wo2fill_between_x = [wo2ddPCR, wo2minit, wo2minit, wo2ddPCR]
	wo2fill_between_y = [0, 0, 2*maximum(wo2hist.weights), 2*maximum(wo2hist.weights)]
	plot!(wo2fill_between_x, wo2fill_between_y, fill=(0, 0.4, 1))
	# savefig("../Plots/CopyHist_PreLabel_withoutAmp_R02.png")
end


# ╔═╡ 071c0c0a-bc39-46bf-86bc-a88338e5e060
begin
	w1minit = mean(w1uniqueinit)
	w1hist = fit(Histogram, w1uniqueinit, nbins=80)
	bar(w1hist.edges[1][1:end-1], w1hist.weights, legend=false, linecolor=:match, xlim=(-0.7, w1hist.edges[1][end]+2), ylim=(0, 1.05*maximum(w1hist.weights)), size=(600, 200), tickfontsize=16, xticks=[0, 20, 40, 60], yticks=[0, 800, 1600])
	vline!([w1minit], line=(2, 4))
	vline!([w1ddPCR], line=(:black, 4))

	# Fill the region between a and b
	w1fill_between_x = [w1ddPCR, w1minit, w1minit, w1ddPCR]
	w1fill_between_y = [0, 0, 2*maximum(w1hist.weights), 2*maximum(w1hist.weights)]
	plot!(w1fill_between_x, w1fill_between_y, fill=(0, 0.4, 1))
	# savefig("../Plots/CopyHist_PreLabel_withAmp_R01.png")
end


# ╔═╡ ac85e63f-df7c-4af1-9331-348c571bef7a
begin
	w2minit = mean(w2uniqueinit)
	w2hist = fit(Histogram, w2uniqueinit, nbins=80)
	bar(w2hist.edges[1][1:end-1], w2hist.weights, legend=false, linecolor=:match, xlim=(-0.7, w2hist.edges[1][end]+2), ylim=(0, 1.05*maximum(w2hist.weights)), size=(600, 200), tickfontsize=16, xticks=[0, 20, 40, 60], yticks=[0, 800, 1600])
	vline!([w2minit], line=(2, 4))
	vline!([w2ddPCR], line=(:black, 4))

	# Fill the region between a and b
	w2fill_between_x = [w2ddPCR, w2minit, w2minit, w2ddPCR]
	w2fill_between_y = [0, 0, 2*maximum(w2hist.weights), 2*maximum(w2hist.weights)]
	plot!(w2fill_between_x, w2fill_between_y, fill=(0, 0.4, 1))
	# savefig("../Plots/CopyHist_PreLabel_withAmp_R02.png")
end


# ╔═╡ a52dde76-615a-4a19-832d-9e3ce635bc75
begin
	wo1uniquefinal = vcat([sum(val != 0 for val in row) for row in eachrow(wo1labeldf[:, 2:end])], zeros(Int, wo1noBC))
	wo2uniquefinal = vcat([sum(val != 0 for val in row) for row in eachrow(wo2labeldf[:, 2:end])], zeros(Int, wo2noBC))
	
	w1uniquefinal = vcat([sum(val != 0 for val in row) for row in eachrow(w1labeldf[:, 2:end])], zeros(Int, w1noBC))
	w2uniquefinal = vcat([sum(val != 0 for val in row) for row in eachrow(w2labeldf[:, 2:end])], zeros(Int, w2noBC))
end

# ╔═╡ 65d6643b-fea1-4816-8a5c-da485a3e9d5d
begin
	wo1mfinal = mean(wo1uniquefinal)
	wo1histfinal = fit(Histogram, wo1uniquefinal, nbins=80)
	bar(wo1histfinal.edges[1][1:end-1], wo1histfinal.weights, legend=false, linecolor=:match, xlim=(-0.7, 22), ylim=(0, 2600), size=(300, 200), tickfontsize=16, xticks=[0, 10, 20], yticks=[0, 1000, 2000])
	vline!([wo1mfinal], line=(2, 4))
	vline!([wo1ddPCR], line=(:black, 4))

	# Fill the region between a and b
	wo1fill_between_xfinal = [wo1ddPCR, wo1mfinal, wo1mfinal, wo1ddPCR]
	wo1fill_between_yfinal = [0, 0, 2*maximum(wo1histfinal.weights), 2*maximum(wo1histfinal.weights)]
	plot!(wo1fill_between_xfinal, wo1fill_between_yfinal, fill=(0, 0.4, 1))
	# savefig("../Plots/CopyHist_PostLabel_withoutAmp_R01.png")
end

# ╔═╡ 16738e01-a06f-4330-ab9f-80aea36fa226
wo1mfinal

# ╔═╡ ecbe3dc0-af16-452e-8f6d-0a11b0ca22f4
begin
	wo2mfinal = mean(wo2uniquefinal)
	wo2histfinal = fit(Histogram, wo2uniquefinal, nbins=80)
	bar(wo2histfinal.edges[1][1:end-1], wo2histfinal.weights, legend=false, linecolor=:match, xlim=(-0.7, 22), ylim=(0, 2600), size=(300, 200), tickfontsize=16, xticks=[0, 10, 20], yticks=[0, 1000, 2000])
	vline!([wo2mfinal], line=(2, 4))
	vline!([wo2ddPCR], line=(:black, 4))

	# Fill the region between a and b
	wo2fill_between_xfinal = [wo2ddPCR, wo2mfinal, wo2mfinal, wo2ddPCR]
	wo2fill_between_yfinal = [0, 0, 2*maximum(wo2histfinal.weights), 2*maximum(wo2histfinal.weights)]
	plot!(wo2fill_between_xfinal, wo2fill_between_yfinal, fill=(0, 0.4, 1))
	# savefig("../Plots/CopyHist_PostLabel_withoutAmp_R02.png")
end


# ╔═╡ 53532863-e43f-4fec-8f03-05a31c0ebe89
wo2mfinal

# ╔═╡ 6ddefc81-f4aa-4102-bcca-345371e1491e
begin
	w1mfinal = mean(w1uniquefinal)
	w1histfinal = fit(Histogram, w1uniquefinal, nbins=80)
	bar(w1histfinal.edges[1][1:end-1], w1histfinal.weights, legend=false, linecolor=:match, xlim=(-0.7, 22), ylim=(0, 2600), size=(300, 200), tickfontsize=16, xticks=[0, 10, 20], yticks=[0, 1000, 2000])
	vline!([w1mfinal], line=(2, 4))
	vline!([w1ddPCR], line=(:black, 4))

	# Fill the region between a and b
	w1fill_between_xfinal = [w1ddPCR, w1mfinal, w1mfinal, w1ddPCR]
	w1fill_between_yfinal = [0, 0, 2*maximum(w1histfinal.weights), 2*maximum(w1histfinal.weights)]
	plot!(w1fill_between_xfinal, w1fill_between_yfinal, fill=(0, 0.4, 1))
	# savefig("../Plots/CopyHist_PostLabel_withAmp_R01.png")
end

# ╔═╡ cc5792d3-54f1-4f54-bd63-c290cfe00174
w1mfinal

# ╔═╡ af6006d3-614b-4e2b-b9ce-0dbfb0b3ec48
begin
	w2mfinal = mean(w2uniquefinal)
	w2histfinal = fit(Histogram, w2uniquefinal, nbins=80)
	bar(w2histfinal.edges[1][1:end-1], w2histfinal.weights, legend=false, linecolor=:match, xlim=(-0.7, 22), ylim=(0, 2600), size=(300, 200), tickfontsize=16, xticks=[0, 10, 20], yticks=[0, 1000, 2000])
	vline!([w2mfinal], line=(2, 4))
	vline!([w2ddPCR], line=(:black, 4))

	# Fill the region between a and b
	w2fill_between_xfinal = [w2ddPCR, w2mfinal, w2mfinal, w2ddPCR]
	w2fill_between_yfinal = [0, 0, 2*maximum(w2histfinal.weights), 2*maximum(w2histfinal.weights)]
	plot!(w2fill_between_xfinal, w2fill_between_yfinal, fill=(0, 0.4, 1))
	# savefig("../Plots/CopyHist_PostLabel_withAmp_R02.png")
end

# ╔═╡ c1bc847d-b29b-4865-9a49-281adaa19062
w2mfinal

# ╔═╡ 6377a680-8d8e-4675-858f-427bb75db022
md"""
### Count Figures
"""

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

# ╔═╡ 1044643c-884f-4017-a7c8-55a129374cb9
md"""
### Test Gene
"""

# ╔═╡ fdcea879-6ae0-4e74-b931-934e551f2fed
begin
	testTF = selectedw1[20]
	testcountdict = formcountdict(testTF, wo1countdf, noBCcells("withoutAmp","R01"))
	
	testmax = maximum(keys(testcountdict))
	testnumcounts = sum(values(testcountdict))
	testcounts, testcountfreq = formcountfreq(testcountdict)
end

# ╔═╡ 725a2c87-b9c3-418c-8fc1-48c448713400
plot(testcounts, testcountfreq, ylim=(0, 0.003), label=false, xlabel="Counts", ylabel="Frequency")

# ╔═╡ dc1c2c00-1a17-461e-8dfb-3ee6a68c9fb8
testTF

# ╔═╡ e310b53e-de2a-4295-9b69-6b7657e29951
begin
	plotTF = "Plasmid_47"
	plotcond = "withoutAmp"
	plotrep = "R01"
	plotdf = load_tsv(countpath(plotcond, plotrep))
	plotfit = load_tsv(contfitpath(plotcond, plotrep))
	finalparamdf = load_tsv(finalfitpath(plotcond, plotrep))
end

# ╔═╡ dce8b0ee-6bf9-41da-8096-30cf02139cbc
begin
	plotcountdict = formcountdict(plotTF, plotdf, noBCcells(plotcond, plotrep))
	plotmax = maximum(keys(plotcountdict))
	plotnumcounts = sum(values(plotcountdict))
	plotcounts, plotcountfreq = formcountfreq(plotcountdict)
end

# ╔═╡ 4f3b95e9-ca2a-41c8-a0e8-4a149642e4b8
begin
	plot(plotcounts, plotcountfreq, ylim=(0, 1), xlim=(0, plotmax+10), label=false, color = 1, lw=10, tickfontsize=16)
	plot!(size=(600, 300))
	# savefig("../Plots/AllCounts_$(plotcond)_$(plotrep)_$(plotTF).png")
end

# ╔═╡ 41929758-8a1a-4e5d-8078-36e715fbb2a9
begin
	t = 7
	xlim = (0, 140)
	xticks = [0, 100, 200]
	ylim = (0, 0.003)
	yticks = [0.0, 0.001, 0.002, 0.003]
	plot(plotcounts[1:t], plotcountfreq[1:t], fillrange=0, fillalpha=0.3, fillcolor=4, label=false)
	plot!(plotcounts, plotcountfreq, label=false, color = 1, lw=4, tickfontsize=16)  
	
	plot!(xlim=xlim,xticks = xticks, ylim=ylim, yticks=yticks, size=(600, 300))
	# savefig("../Plots/CountsZoomed_$(plotcond)_$(plotrep)_$(plotTF).png")
end

# ╔═╡ eb493270-6eec-4883-9237-c0b8862a332b
plotcountfreq[1], sum(plotcountfreq[2:t]), sum(plotcountfreq[t+1:end])

# ╔═╡ 1fced61e-e061-42e5-b153-0d64f7d3d4b8
begin
	ν, γ, ρμ, α, f, threshold = plotfit[:, plotTF]
	
	gex(z) = GNB(z, ρμ, α)
	gexν(z) = GNB(z, ρμ*ν, α)
	gcont(z) = Gpois(gexν(z), f*γ)
	gseq(z) = gcont(z)*(f*gex(z) + (1-f))

	pnoexlist = taylor_expand(z->gcont(z)*gex(z), order = plotcounts[end]).coeffs
    pseqlist = taylor_expand(gseq, order = plotcounts[end]).coeffs
    pnolistlong = taylor_expand(gcont, order = plotcounts[end]).coeffs
end

# ╔═╡ 3db6e14f-aca3-44ed-bd3c-12319468b1ed
begin
	plot(plotcounts, plotcountfreq, color = 1, lw=3, label="Counts")

	plot!(plotcounts, (1-f).*pnolistlong, lw=3, label="Contaminant")
	plot!(plotcounts, f.*pnoexlist, lw=3, label="Exp+Cont", c=5)
	
	# plot!(plotcounts, pseqlist.+4e-5, lw=3, label="Total")
	
	plot!(xlim=xlim, xticks = xticks, ylim=ylim, yticks=yticks, size=(600, 300), tickfontsize=16, legendfontsize=12)
	# savefig("../Plots/Fit_Pre_$(plotcond)_$(plotrep)_$(plotTF).png")
end

# ╔═╡ 3f9816ad-7661-45e4-88fe-6af36353720b
md"""
### Final Fits
"""

# ╔═╡ aa008f08-94dd-4bfb-9648-e6b247a9b7dc
begin
	fitparam = finalparamdf[1:5, plotTF]
	ff = fitparam[end]
	gexf, gcontf, gseqf = genfuncs(fitparam...)

	pnoexf = taylor_expand(z->gcontf(z)*gexf(z), order = plotcounts[end]).coeffs
    pseqf = taylor_expand(gseqf, order = plotcounts[end]).coeffs
    pnof = taylor_expand(gcontf, order = plotcounts[end]).coeffs
end

# ╔═╡ 1e595a25-dd80-4878-8cad-5666b9f5a817
begin
	plot(plotcounts, plotcountfreq, color = 1, lw=3, label="Counts")

	plot!(plotcounts, (1-ff).*pnof, lw=3, label="Contaminant")
	plot!(plotcounts, ff.*pnoexf, lw=3, label="Exp+Cont", c=5)
	
	# plot!(plotcounts, pseqlist.+4e-5, lw=3, label="Total")
	
	plot!(xlim=xlim, xticks = xticks, ylim=ylim, yticks=yticks, size=(600, 300), tickfontsize=16, legendfontsize=12)
	# savefig("../Plots/Fit_Post_$(plotcond)_$(plotrep)_$(plotTF).png")
end

# ╔═╡ ebdf1dec-89e3-41ef-bf85-69e112c420fe
begin
	BClistlow = [n for n in names(finalparamdf) if n ∉ union(selectedw1, selectedw2, selectedwo1, selectedwo2)]
	
	testTFlow = BClistlow[5]	
end

# ╔═╡ efea7b8a-ab20-4aef-8156-e3f664a3dab7
begin
	plotcountdictl = formcountdict(testTFlow, plotdf, noBCcells(plotcond, plotrep))
	plotmaxl = maximum(keys(plotcountdictl))
	plotnumcountsl = sum(values(plotcountdictl))
	plotcountsl, plotcountfreql = formcountfreq(plotcountdictl)
end

# ╔═╡ 928e4c83-67a7-4c1c-8c82-0552f3348d82
begin
	fitparaml = finalparamdf[1:5, testTFlow]
	fl = fitparaml[end]
	gexl, gcontl, gseql = genfuncs(fitparaml...)

	pnoexl = taylor_expand(z->gcontl(z)*gexl(z), order = plotcountsl[end]).coeffs
    pseql = taylor_expand(gseql, order = plotcountsl[end]).coeffs
    pnol = taylor_expand(gcontl, order = plotcountsl[end]).coeffs
end

# ╔═╡ 78084e5d-1b00-4740-ada8-322cb53a3bdb
begin
	xliml = (0, 42)
	xticksl = [0, 20]
	yliml = (0, 0.01)
	yticksl = [0, 0.025, 0.05]
	plot(plotcountsl, plotcountfreql, color = 1, lw=3, label="Counts")

	plot!(plotcountsl, (1-fl).*pnol, lw=3, label="Contaminant")
	plot!(plotcountsl, fl.*pnoexl, lw=3, label="Exp+Cont", c=5)
	
	# plot!(plotcounts, pseqlist.+4e-5, lw=3, label="Total")
	
	plot!(xlim=xliml, xticks = xticksl, ylim=yliml, yticks=yticksl, size=(300, 300), tickfontsize=16, legendfontsize=12, legend = false)
	# savefig("../Plots/Fit_PostLow_$(plotcond)_$(plotrep)_$(testTFlow).png")
end

# ╔═╡ d8f00d51-b809-4b37-a05c-8450998c58a9
md"""
### Final Param
"""

# ╔═╡ 5002c70a-9a59-497d-afed-9dd1342c43a0
begin
	
	wo1fparam = load_tsv(finalfitpath("withoutAmp", "R01"))
	wo2fparam = load_tsv(finalfitpath("withoutAmp", "R02"))

	w1fparam = load_tsv(finalfitpath("withAmp", "R01"))
	w2fparam = load_tsv(finalfitpath("withAmp", "R02"))

	allplasmids = intersect(names(wo1fparam), names(wo2fparam), names(w1fparam), names(w2fparam))
end

# ╔═╡ 9b0c5375-f051-47fb-ac9c-232cb0a36391
copy = sum(-log.(1 .- Vector(w1fparam[5, :])))

# ╔═╡ dbeb996f-9590-4c6d-9cc6-4bf78262e80a
histogram(Vector(wo1fparam[6,:]))

# ╔═╡ 4bb8ac4e-b80a-4f09-bddc-6da3e4112b87
histogram(Vector(wo2fparam[6,:]))

# ╔═╡ 397ce2d2-2d43-4d84-a512-4ba6b714c6c0
begin
	cmap1 = countmap(wo1fparam[6,:])
	cmap2 = countmap(wo2fparam[6,:])
	
	thresh1 = sort(collect(keys(cmap1)))
	thresh2 = sort(collect(keys(cmap2)))

	allthresh = convert(Array{Integer}, sort(union(thresh1, thresh2)))
	
	nthresh1 = size(wo1fparam, 2)
	nthresh2 = size(wo2fparam, 2)
	
	freq1 = [get(cmap1, t, 0)/nthresh1 for t in allthresh]
	freq2 = [get(cmap2, t, 0)/nthresh2 for t in allthresh]

	groupedbar(allthresh, [freq1 freq2], lw=0, legend=false, tickfontsize=16, xticks=1:maximum(allthresh), size=(350, 300))
end

# ╔═╡ 551a80ed-1d48-4ada-91f1-9ea92a6bc6c6
begin
	woBClist = intersect(names(wo1fparam), names(wo2fparam))
	scatter(Vector(wo1fparam[6,woBClist]), Vector(wo2fparam[6,woBClist]))
end

# ╔═╡ ff7135f3-80e2-4728-ac16-ff70ef3ce037
begin
	cmap1a = countmap(w1fparam[6,:])
	cmap2a = countmap(w2fparam[6,:])
	
	thresh1a = sort(collect(keys(cmap1a)))
	thresh2a = sort(collect(keys(cmap2a)))

	allthresha = convert(Array{Integer}, sort(union(thresh1a, thresh2a)))
	
	nthresh1a = size(w1fparam, 2)
	nthresh2a = size(w2fparam, 2)
	
	freq1a = [get(cmap1a, t, 0)/nthresh1a for t in allthresha]
	freq2a = [get(cmap2a, t, 0)/nthresh2a for t in allthresha]

	groupedbar(allthresha, [freq1a freq2a], lw=0, legend=false, tickfontsize=16, size=(900, 300), barposition=:dodge)
end

# ╔═╡ cfa31038-2057-4708-8e99-7b3544849351
begin
	wBClist = intersect(names(w1fparam), names(w2fparam))
	
	scatter(Vector(w1fparam[6,wBClist]), Vector(w2fparam[6,wBClist]), alpha = 0.5, markerstrokewidth=0, color=2, label="with Amplification")
	scatter!(Vector(wo1fparam[6,woBClist]), Vector(wo2fparam[6,woBClist]), alpha = 0.6, markerstrokewidth=0, color=1, label="w/o Amplification")
	
	plot!(tickfontsize=16, legendfontsize=12, size = (500, 250))
	# savefig("../Plots/Threshold_corr.png")
end

# ╔═╡ 3f8d55e0-0a61-41f1-a2be-9da0f29bea84
md"""
### ROC curves, Error Rates, etc.
"""

# ╔═╡ 7a11e58e-01ca-40b3-8c9a-c4a06b3fde19
begin
	ROCrep = "R01"
	wROCparamdf = load_tsv(finalfitpath("withAmp", ROCrep))
	wROCcountdf = load_tsv(countpath("withAmp", ROCrep))
	wROCBClist = names(wROCparamdf)

	woROCparamdf = load_tsv(finalfitpath("withoutAmp", ROCrep))
	woROCcountdf = load_tsv(countpath("withoutAmp", ROCrep))
	woROCBClist = names(woROCparamdf)
end

# ╔═╡ 817e41a9-2ea0-4895-bb4c-26efc676ed38
begin
	function ROCcurve(ROCBC, ROCcountdf, ROCparamdf)
		maxcount = maximum(ROCcountdf[:,ROCBC])
		ROCparam = ROCparamdf[1:5, ROCBC]
		
		ROCgex, ROCgcont, ROCgseq = genfuncs(ROCparam...)
		ROCgexco(z) = ROCgex(z)*ROCgcont(z)
	
		pexco = taylor_expand(ROCgexco, order = maxcount).coeffs
		pco = taylor_expand(ROCgcont, order = maxcount).coeffs
		pexco[pexco .<=0] .= 0
		pco[pco .<=0] .= 0
	
		TPR = [1-sum(pexco[1:t]) for t in 0:maxcount]
	
		FPR = [1-sum(pco[1:t]) for t in 0:maxcount]

		return TPR, FPR
	end
	
end

# ╔═╡ 3ff602b5-a54d-4f9d-87c9-e036caa5c6a1
begin

	TPR1, FPR1 = ROCcurve(woROCBClist[1], woROCcountdf, woROCparamdf)
	fig = plot(FPR1, TPR1, color = 1, alpha = 0.3, label = "w/o Amplification")

	for ROCBC in woROCBClist[2:end]
		TPR, FPR = ROCcurve(ROCBC, woROCcountdf, woROCparamdf)
		plot!(fig, FPR, TPR, color = 1, alpha = 0.3, label = false)
	end

	TPR2, FPR2 = ROCcurve(wROCBClist[1], wROCcountdf, wROCparamdf)
	plot!(fig, FPR2, TPR2, color = 2, alpha = 0.2, label="with Amplification")

	for ROCBC in wROCBClist[2:end]
		TPR, FPR = ROCcurve(ROCBC, wROCcountdf, wROCparamdf)
		plot!(fig, FPR, TPR, color = 2, alpha = 0.2, label = false)
	end
	plot!(fig, tickfontsize = 16, legendfontsize = 12, legendtraceorder="reversed", size = (300, 200), xticks = [0, 0.5, 1], yticks = [0, 0.5, 1])
	# savefig("../Plots/ROCcurve_$(ROCrep).png")
end

# ╔═╡ 8424f5bf-2aee-4050-a487-c9e3d0ecbb35
begin
	woFPfrac = []
	for ROCBC in woROCBClist
		ROCparam = woROCparamdf[1:5, ROCBC]
		fROC = ROCparam[end]
			
		ROCgex, ROCgcont, ROCgseq = genfuncs(ROCparam...)
		FPprob = (1-fROC)*(1-ROCgcont(0))
		TPprob = fROC*(1-ROCgex(0)*ROCgcont(0))
		push!(woFPfrac, FPprob/(FPprob+TPprob))
	end

	wFPfrac = []
	for ROCBC in wROCBClist
		ROCparam = wROCparamdf[1:5, ROCBC]
		fROC = ROCparam[end]
			
		ROCgex, ROCgcont, ROCgseq = genfuncs(ROCparam...)
		FPprob = (1-fROC)*(1-ROCgcont(0))
		TPprob = fROC*(1-ROCgex(0)*ROCgcont(0))
		push!(wFPfrac, FPprob/(FPprob+TPprob))
	end
	histogram(woFPfrac, nbins = 20, alpha = 0.9, lw=0, label = "w/o Amplification")
	histogram!(wFPfrac, nbins=20, alpha = 0.7, lw=0, label = "with Amplification")
	plot!(tickfontsize = 16, legendfontsize = 12, size =(300, 200), xticks = [0, 0.4, 0.8], yticks = [0, 10, 20])

	# savefig("../Plots/FPfrac_$(ROCrep).png")
end

# ╔═╡ Cell order:
# ╠═9f56ed72-a60a-11ed-0a9c-078cd1e09182
# ╠═04c7ebc3-d911-4f6d-ac1d-c21262506972
# ╠═7eb6af35-b356-4b9d-a203-22913e239eab
# ╟─adffddb0-d063-490a-a17d-4587fa350749
# ╠═ab46ed63-338f-4297-94d1-0bda0b9a9704
# ╟─33361e24-b008-4996-b36e-90ffbd1c9874
# ╠═1c1e08d4-dc8d-46bb-a275-b8146b1f75c0
# ╠═87d16f2c-b2a9-44e6-a4ec-88eda65a8739
# ╠═dc7489d6-35f5-46f5-a978-e16abbf8072d
# ╠═95ee2eb2-a3f1-4e0a-99d9-3913d6334e71
# ╠═29302f0f-4ccb-4905-a8d4-fc79934809df
# ╟─b13cb5c8-e27b-47ef-b523-a95c73b4a23f
# ╠═68822eff-a59b-486b-9ac3-5c9ae7af734d
# ╠═aed79458-b324-4662-a1c3-7368768a4047
# ╠═972a908d-5fa8-406e-8832-4f30af6ffdb6
# ╠═e2f55db3-c8d5-4f69-b540-59e7ca551abd
# ╠═653b2c74-2dca-4988-865c-02d99b1ee49f
# ╟─293a0eb4-d300-4e67-8388-dd489e2dd427
# ╠═d7d7ebda-ffe7-4a58-9362-7df922c5c2fd
# ╠═f36c3230-92cd-4c9e-b01c-590bbe5776bd
# ╠═10cdb8f2-95c2-429f-b7e4-a39a99f6db65
# ╠═bdd3a1f1-8e69-4559-8f35-4c8182340425
# ╠═277256d1-b2e4-4d52-8efd-b74b477a1fb4
# ╠═8aa03f0c-a67a-42c1-a913-8a7a8f591b2f
# ╟─7f274c2e-dc27-4cb4-92dc-4a649ca3e3e3
# ╠═ff34c6db-35f8-4e36-9f6c-f75227b39c5c
# ╠═ece369cf-5820-42c5-b4ea-5c7fa9f6cd01
# ╠═7332a6d2-ed2c-4619-85fc-c1aa26f228e6
# ╠═5ef3aa92-7638-464c-8714-1fe78f6aaa3c
# ╠═01ea329b-4781-401b-aa8d-a7ee347a5b81
# ╠═d2858b3d-e9b1-49d5-8310-136fb8cba12f
# ╠═a55d2c29-284e-418c-a683-c793d3db8444
# ╠═43c4454f-b901-4a53-a9fc-6fe0e32520db
# ╠═b8fd4822-e4d0-4cf8-9efa-589f9205ae28
# ╠═f030de7f-2af4-4dff-bca5-caff7ac3a0fc
# ╠═a29d9e64-6325-442e-9d9a-2e26730f818c
# ╠═3ca5df27-ae5e-423d-9228-b67389132215
# ╟─c3828bf3-8401-477e-899a-b6ac375552ff
# ╠═ab30e5f7-a0ea-4e5e-8124-7a1052eef7e5
# ╠═523a26af-9b92-45a3-9f2f-8b6f50b5a790
# ╠═8a6b0fe0-9be0-4b3d-80c7-ac4bae5ac822
# ╠═c6757638-8730-457e-9dd5-0c12c78ba4cf
# ╠═c484975f-8bc9-413b-b28a-35947491cf60
# ╠═d09e3965-d0fc-499e-b303-804ec35e7023
# ╠═071c0c0a-bc39-46bf-86bc-a88338e5e060
# ╠═ac85e63f-df7c-4af1-9331-348c571bef7a
# ╠═a52dde76-615a-4a19-832d-9e3ce635bc75
# ╠═65d6643b-fea1-4816-8a5c-da485a3e9d5d
# ╠═16738e01-a06f-4330-ab9f-80aea36fa226
# ╠═ecbe3dc0-af16-452e-8f6d-0a11b0ca22f4
# ╠═53532863-e43f-4fec-8f03-05a31c0ebe89
# ╠═6ddefc81-f4aa-4102-bcca-345371e1491e
# ╠═cc5792d3-54f1-4f54-bd63-c290cfe00174
# ╠═af6006d3-614b-4e2b-b9ce-0dbfb0b3ec48
# ╠═c1bc847d-b29b-4865-9a49-281adaa19062
# ╟─6377a680-8d8e-4675-858f-427bb75db022
# ╠═70400124-5b17-4ad6-bddb-4a2715d0d7d2
# ╟─1044643c-884f-4017-a7c8-55a129374cb9
# ╠═fdcea879-6ae0-4e74-b931-934e551f2fed
# ╠═725a2c87-b9c3-418c-8fc1-48c448713400
# ╠═dc1c2c00-1a17-461e-8dfb-3ee6a68c9fb8
# ╠═e310b53e-de2a-4295-9b69-6b7657e29951
# ╠═dce8b0ee-6bf9-41da-8096-30cf02139cbc
# ╠═4f3b95e9-ca2a-41c8-a0e8-4a149642e4b8
# ╠═41929758-8a1a-4e5d-8078-36e715fbb2a9
# ╠═eb493270-6eec-4883-9237-c0b8862a332b
# ╠═1fced61e-e061-42e5-b153-0d64f7d3d4b8
# ╠═3db6e14f-aca3-44ed-bd3c-12319468b1ed
# ╟─3f9816ad-7661-45e4-88fe-6af36353720b
# ╠═aa008f08-94dd-4bfb-9648-e6b247a9b7dc
# ╠═1e595a25-dd80-4878-8cad-5666b9f5a817
# ╠═ebdf1dec-89e3-41ef-bf85-69e112c420fe
# ╠═efea7b8a-ab20-4aef-8156-e3f664a3dab7
# ╠═928e4c83-67a7-4c1c-8c82-0552f3348d82
# ╠═78084e5d-1b00-4740-ada8-322cb53a3bdb
# ╟─d8f00d51-b809-4b37-a05c-8450998c58a9
# ╠═5002c70a-9a59-497d-afed-9dd1342c43a0
# ╠═9b0c5375-f051-47fb-ac9c-232cb0a36391
# ╠═dbeb996f-9590-4c6d-9cc6-4bf78262e80a
# ╠═4bb8ac4e-b80a-4f09-bddc-6da3e4112b87
# ╠═397ce2d2-2d43-4d84-a512-4ba6b714c6c0
# ╠═551a80ed-1d48-4ada-91f1-9ea92a6bc6c6
# ╠═ff7135f3-80e2-4728-ac16-ff70ef3ce037
# ╠═cfa31038-2057-4708-8e99-7b3544849351
# ╠═3f8d55e0-0a61-41f1-a2be-9da0f29bea84
# ╠═7a11e58e-01ca-40b3-8c9a-c4a06b3fde19
# ╠═817e41a9-2ea0-4895-bb4c-26efc676ed38
# ╠═3ff602b5-a54d-4f9d-87c9-e036caa5c6a1
# ╠═8424f5bf-2aee-4050-a487-c9e3d0ecbb35
