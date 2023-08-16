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
	using CSV, DataFrames
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

# ╔═╡ 02d9a29a-b282-48c2-a883-c06295f576b9
begin
	replist = ["R01", "R02"]
	condlist = ["withoutAmp", "withAmp"]
end

# ╔═╡ 7dfdca3e-ad92-4c15-8709-fa32b3496197
begin
	datapath = "../Data"
	withoutAmpPath = joinpath(datapath, "withoutAmp_R01R02")
	withAmpPath = joinpath(datapath, "withAmp_R01R02")
	
	wo1celldf = DataFrame(CSV.File(joinpath(withoutAmpPath, "withoutAmp_R01_CellCounts.tsv"), delim="\t"))
	wo1noBCcells = parse(Int, readline(joinpath(withoutAmpPath, "withoutAmp_R01_noBCcells.txt")))
	# wo1finalcelldf = DataFrame(CSV.File(joinpath(withoutAmpPath, "withoutAmp_R01_Labeled.tsv"), delim="\t"))

	wo2celldf = DataFrame(CSV.File(joinpath(withoutAmpPath, "withoutAmp_R02_CellCounts.tsv"), delim="\t"))
	wo2noBCcells = parse(Int, readline(joinpath(withoutAmpPath, "withoutAmp_R02_noBCcells.txt")))
	# wo2finalcelldf = DataFrame(CSV.File(joinpath(withoutAmpPath, "withoutAmp_R02_Labeled.tsv"), delim="\t"))
	
	post1celldf = DataFrame(CSV.File(joinpath("./Data/PostAmp_R01R02/PostAmp_R01_CellCounts.tsv"), delim="\t"))
	post1noBCcells = parse(Int, readline("./Data/PostAmp_R01R02/PostAmp_R01_noBCcells.txt"))
	post1finalcelldf = DataFrame(CSV.File("./Data/PostAmp_R01R02/PostAmp_R01_Labeled.tsv", delim="\t"))
	
	post2celldf = DataFrame(CSV.File(joinpath("./Data/PostAmp_R01R02/PostAmp_R02_CellCounts.tsv"), delim="\t"))
	post2noBCcells = parse(Int, readline("./Data/PostAmp_R01R02/PostAmp_R02_noBCcells.txt"))
	post2finalcelldf = DataFrame(CSV.File("./Data/PostAmp_R01R02/PostAmp_R02_Labeled.tsv", delim="\t"))

	R01ddPCR = parse(Float64, readline("./Data/PreAmp_R01R02/R01_ddPCR.txt"))
	R02ddPCR = parse(Float64, readline("./Data/PreAmp_R01R02/R02_ddPCR.txt"))
end

# ╔═╡ 50580bae-24ff-40c8-98f9-5600589572fc
R01ddPCR, R02ddPCR

# ╔═╡ 86042e1e-3ed1-4fd4-8296-8f25c61c3f16
md"""
### Summary Figures
"""

# ╔═╡ 96f7230d-1f62-45f8-91e8-0b4f7dd28348
begin
	wo1noisefitdf = DataFrame(CSV.File(joinpath(withoutAmpPath, "GeneFitData", "contaminantfits","withoutAmp_R01_Parameters.tsv"), delim="\t"))
	w1noisefitdf = DataFrame(CSV.File(joinpath(withAmpPath, "GeneFitData", "contaminantfits","withAmp_R01_Parameters.tsv"), delim="\t"))

	wo2noisefitdf = DataFrame(CSV.File(joinpath(withoutAmppath, "GeneFitData", "contaminantfits","withoutAmp_R02_Parameters.tsv"), delim="\t"))
	w2noisefitdf = DataFrame(CSV.File(joinpath(withAmpPath, "GeneFitData", "contaminantfits","withAmp_R02_Parameters.tsv"), delim="\t"))
end

# ╔═╡ 4778354d-4208-4b65-ba8f-67f5f86233b9
pre1noisefitdf

# ╔═╡ 5961d6a2-4b13-4865-b6be-0e8ece20d5f2
post1noisefitdf

# ╔═╡ e0a3512d-2b66-4f5c-9bda-8c7a0712ecbf
begin
	pre1names = names(pre1noisefitdf)
	pre2names = names(pre2noisefitdf)
	
	pre1νvalues = Vector(pre1noisefitdf[1,:])
	pre1γvalues = Vector(pre1noisefitdf[2,:])

	pre2νvalues = Vector(pre2noisefitdf[1,:])
	pre2γvalues = Vector(pre2noisefitdf[2,:])
	
	pre1prod = median(pre1νvalues.*pre1γvalues)
	pre2prod = median(pre2νvalues.*pre2γvalues)


	post1names = names(post1noisefitdf)
	post2names = names(post2noisefitdf)
	
	post1νvalues = Vector(post1noisefitdf[1,:])
	post1γvalues = Vector(post1noisefitdf[2,:])

	post2νvalues = Vector(post2noisefitdf[1,:])
	post2γvalues = Vector(post2noisefitdf[2,:])

	post1prod = median(post1νvalues.*post1γvalues)
	post2prod = median(post2νvalues.*post2γvalues)
end

# ╔═╡ 0c15417f-d723-4672-8039-c33564b12f61
colorpalette = theme_palette(:auto)

# ╔═╡ b501569c-59e7-46d4-a059-02708c44eb99
bluegrad = cgrad([:white, colorpalette[1], :black])

# ╔═╡ 07dd1ce6-7791-4b37-9695-24ecb8dd6fcb
myblues = cgrad([get(bluegrad, 0.1), colorpalette[1], get(bluegrad, 0.8)])

# ╔═╡ 35134891-e3dd-459c-8703-97235be886ea
orangegrad = cgrad([:white, colorpalette[2], :black])

# ╔═╡ b09d1534-311a-4fd2-83d2-d1a9628b0c34
myoranges = cgrad([get(orangegrad, 0.1), colorpalette[2], get(orangegrad, 0.8)])

# ╔═╡ ff72b3c2-c4e1-458a-9f57-34fbcd694327
begin
	scatter(pre1νvalues, pre1γvalues, markerstrokecolor=:auto, label="Replicate 1")
	linebounds = 10^-4:10^-5:0.0047
	plot!(linebounds, pre1prod./linebounds, color=1, lw=2, label=false)

	scatter!(pre2νvalues, pre2γvalues, color=2, markerstrokecolor=:auto, label="Replicate 2", tickfontsize=16, xticks=[0, 0.002, 0.004], yticks=[0, 100, 200], xlim=[0, 0.005], ylim = [0, 200], size=(600, 350))
	plot!(linebounds, pre2prod./linebounds, color=2, lw=2, label=false, legendfontsize=12)
	# savefig("./Plots/PreScatter.png")
end

# ╔═╡ 487249ca-ddcf-4dbd-8c79-dde5639cfb29
begin
	
	pre1numcounts = sum.(eachcol(pre1celldf[!,pre1names]))
	# pre1numcounts = sum.(c .!= 0 for c in eachcol(pre1celldf[!,pre1names]))
	pre1perm = sortperm(pre1numcounts)
	pre2numcounts = sum.(eachcol(pre2celldf[!,pre2names]))
	# pre2numcounts = sum.(c .!= 0 for c in eachcol(pre2celldf[!,pre2names]))
	pre2perm = sortperm(pre2numcounts)
	
	scatter(pre1νvalues[pre1perm], pre1γvalues[pre1perm], zcolor=log.(pre1numcounts[pre1perm]), colorbar=false, c=myblues, markerstrokewidth=0, xlim=(0, 0.005))
	# plot!(linebounds, pre1prod./linebounds, color=1, lw=2)


	scatter!(pre2νvalues[pre2perm], pre2γvalues[pre2perm], zcolor=log.(pre2numcounts[pre2perm]), colorbar=false, c=myoranges, markerstrokewidth=0, alpha=0.7)
	# plot!(linebounds, pre2prod./linebounds, color=2, lw=2)
	plot!(size=(800, 500), ylims=(0, 200), xlim=(0, 0.0047), legend=false, tickfontsize=16)

end

# ╔═╡ 7c18dc4c-fab0-48ba-8693-1e0158cb029e
# begin
# 	minalpha = 0.3

# 	pre1alpha = pre1numcounts[pre1perm]
# 	pre1alpha = pre1alpha .- pre1alpha[1]
# 	pre1alpha = pre1alpha ./ pre1alpha[end]
# 	pre1alpha = minalpha .+ pre1alpha.*(1-minalpha)
	
	
# 	scatter(pre1νvalues[pre1perm], pre1γvalues[pre1perm], alpha=pre1alpha, c=1, markerstrokewidth=0, xlim=(0, 0.0047), ylim=(0, 200))


# 	pre2alpha = pre2numcounts[pre2perm]
# 	pre2alpha = pre2alpha .- pre2alpha[1]
# 	pre2alpha = pre2alpha ./ pre2alpha[end]
# 	pre2alpha = minalpha .+ pre2alpha.*(1-minalpha)
	
# 	scatter!(pre2νvalues[pre2perm], pre2γvalues[pre2perm], alpha=pre2alpha, markerstrokewidth=0)

# end

# ╔═╡ 88e75c2b-7b0f-4625-979c-6b48bd034d06
# begin
# 	scatter(1 ./pre1numcounts, pre1γvalues, xlim=(0, 0.00011), ylim=(0, 100), markerstrokewidth=0)
# 	scatter!(1 ./pre2numcounts, pre2γvalues, xlim=(0, 0.00011), ylim=(0, 100), markerstrokewidth=0)
# end

# ╔═╡ 78ca2110-df58-4032-a5f5-21d6a9fcae59
median(pre1γvalues), median(pre2γvalues)

# ╔═╡ c9b28149-b848-4cae-8de7-70533628c7dd
# begin
# 	scatter(1 ./pre1numcounts, pre1νvalues, markerstrokewidth=0)
# 	scatter!(1 ./pre2numcounts, pre2νvalues, markerstrokewidth=0)
# 	plot!(xlim=(0, 0.00011), ylim=(0, 0.006))
# end

# ╔═╡ b98ae866-3d06-4126-b58d-537eaab6b504
begin
	scatter(post1νvalues, post1γvalues, markerstrokecolor=:auto)
	
	scatter!(post2νvalues, post2γvalues, markerstrokecolor=:auto)
end

# ╔═╡ d35b5774-cb3d-4214-9a70-48dbad96ca84
median(post1γvalues), median(post2γvalues)

# ╔═╡ 054d8a56-a86d-4883-be0e-901b01893570
begin
	post1numcounts = sum.(eachcol(post1celldf[!,post1names]))
	# post1numcounts = sum.(c .!= 0 for c in eachcol(post1celldf[!,post1names]))

	post1perm = sortperm(post1numcounts)
	
	post2numcounts = sum.(eachcol(post2celldf[!,post2names]))
	# post2numcounts = sum.(c .!= 0 for c in eachcol(post2celldf[!,post2names]))

	post2perm = sortperm(post2numcounts)
	
	scatter(post1νvalues[post1perm], post1γvalues[post1perm], zcolor=log.(post1numcounts[post1perm]), colorbar=false, c=myblues, markerstrokewidth=0)
end

# ╔═╡ 7c28ffc9-3720-4808-ae17-cabad2d89c68
# begin
# 	scatter(1 ./post1numcounts, post1γvalues, ylim=(0, 50), markerstrokewidth=0)
# 	scatter!(1 ./post2numcounts, post2γvalues, ylim=(0, 50), markerstrokewidth=0)
# end

# ╔═╡ c1e1a053-6d49-4ae8-a1d1-c6fa5f712c3c
# begin
# 	scatter(1 ./post1numcounts, post1νvalues, markerstrokewidth=0)
# 	scatter!(1 ./post2numcounts, post2νvalues, markerstrokewidth=0)
# 	plot!(xlim=(0, 0.00014), ylim=(0, 0.02))
# end

# ╔═╡ 60f09855-a237-4259-8c6d-98eaf946c40d
# begin
# 	post1alpha = post1numcounts[post1perm]
# 	post1alpha = post1alpha .- post1alpha[1]
# 	post1alpha = post1alpha ./ post1alpha[end]
# 	post1alpha = minalpha .+ post1alpha.*(1-minalpha)

# 	post2alpha = post2numcounts[post2perm]
# 	post2alpha = post2alpha .- post2alpha[1]
# 	post2alpha = post2alpha ./ post2alpha[end]
# 	post2alpha = minalpha .+ post2alpha.*(1-minalpha)
# end

# ╔═╡ 09fbc6b9-94e4-4109-9173-571bf648f35a
begin
	scatter(pre1νvalues[pre1perm], pre1γvalues[pre1perm], label=nothing, zcolor=log.(pre1numcounts[pre1perm]), color=myblues, markerstrokewidth=0)
	pre1linebounds = [1.1*10^-4, 4*10^-3]
	plot!(pre1linebounds, pre1prod./pre1linebounds, color = 1, lw=2, xlim=(8*10^-5, 2*10^-2), ylim=(6, 2.5*10^2), label=:none)
	# scatter!([NaN], [NaN], markerstrokewidth=0, color=colorpalette[1])  # Add dummy series


	scatter!(pre2νvalues[pre2perm], pre2γvalues[pre2perm], label=nothing, zcolor=log.(pre2numcounts[pre2perm]), color=myoranges, markerstrokewidth=0)
	pre2linebounds = [10^-4, 3*10^-3]
	plot!(pre2linebounds, pre2prod./pre2linebounds, color = 2, lw=2, label=:none)
	# scatter!([NaN], [NaN], markerstrokewidth=0, color=colorpalette[2])  # Add dummy series
	
	scatter!(post1νvalues, post1γvalues,xscale=:log10, yscale=:log10, label=:none, zcolor=log.(post1numcounts[post1perm]), color=myblues, markerstrokewidth=0)
	post1linebounds = [6*10^(-3), 2*10^-2]
	plot!(post1linebounds, post1prod./post1linebounds, color = 1, lw=2, label=:none)

	scatter!(post2νvalues, post2γvalues, markerstrokecolor=:auto, label=:none, tickfontsize=15, zcolor=log.(post2numcounts[post2perm]), color=myoranges, markerstrokewidth=0, colorbar=:false)
	post2linebounds = [6*10^(-3), 1.6*10^-2]
	plot!(post2linebounds, post2prod./post2linebounds, color = 2, lw=2, label=:none, size=(600, 350), legendfontsize=12)
	
	# savefig("./Plots/LogScatter.png")
end

# ╔═╡ 1cc1ee57-3d3b-48aa-941e-8bf70b3e3132
# begin
# 	pre1post1names = [p for p in pre1names if p∈ post1names]
	
# 	pre1ρμvalues = Vector(pre1noisefitdf[3,pre1post1names])
	
# 	post1ρμvalues = Vector(post1noisefitdf[3,pre1post1names])

# 	scatter(pre1ρμvalues, post1ρμvalues)
# 	plot!([0, 250], [0, 250])
# end

# ╔═╡ 037c7154-c014-476a-9b36-e2605884be04
md"""
### Parameters after Regularisation
"""

# ╔═╡ 7345e9fd-438e-4117-89a8-b7a88cb70dc4
begin
	pre1fitdf = DataFrame(CSV.File(joinpath("./Data/PreAmp_R01R02/", "GeneFitData", "allfits","PreAmp_R01_Parameters.tsv"), delim="\t"))
	post1fitdf = DataFrame(CSV.File(joinpath("./Data/PostAmp_R01R02/", "GeneFitData", "allfits","PostAmp_R01_Parameters.tsv"), delim="\t"))

	pre2fitdf = DataFrame(CSV.File(joinpath("./Data/PreAmp_R01R02/", "GeneFitData", "allfits","PreAmp_R02_Parameters.tsv"), delim="\t"))
	post2fitdf = DataFrame(CSV.File(joinpath("./Data/PostAmp_R01R02/", "GeneFitData", "allfits","PostAmp_R02_Parameters.tsv"), delim="\t"))
end

# ╔═╡ e101d9ac-1b7d-4cb3-896a-edc2b9ecbcde
begin
	pre1allνvalues = Vector(pre1fitdf[1, :])
	pre1allγvalues = Vector(pre1fitdf[2, :])
	pre1allθvalues = Vector(pre1fitdf[6, :])

	pre2allνvalues = Vector(pre2fitdf[1, :])
	pre2allγvalues = Vector(pre2fitdf[2, :])
	pre2allθvalues = Vector(pre2fitdf[6, :])

	post1allνvalues = Vector(post1fitdf[1, :])
	post1allγvalues = Vector(post1fitdf[2, :])
	post1allθvalues = Vector(post1fitdf[6, :])

	post2allνvalues = Vector(post2fitdf[1, :])
	post2allγvalues = Vector(post2fitdf[2, :])
	post2allθvalues = Vector(post2fitdf[6, :])

end

# ╔═╡ 9441c065-2d92-43ad-964a-51f2838971c2
histogram(post2allθvalues)

# ╔═╡ 53c56589-b43e-4e23-9bfb-6b1097715ad0
scatter(Vector(pre1fitdf[5,:]).*Vector(pre1fitdf[3,:]), Vector(pre1fitdf[6,:]))

# ╔═╡ 482c62fa-12f5-48fd-bdb7-9ce62bceaabc
md"""
### Labeling Histograms
"""

# ╔═╡ 1614673f-5bc8-4907-a89c-95d826dbedf7
begin
	celldf = pre1celldf
	noBCcells = pre1noBCcells
	ddPCR = R01ddPCR
	finalcelldf = pre1finalcelldf
	noisefitdf = pre1noisefitdf
	
	numuniqueplasmids = vcat([sum(val != 0 for val in row) for row in eachrow(celldf[:, 2:end])], zeros(Int, noBCcells))
end

# ╔═╡ 777d7c02-a033-4060-a11d-3b70278edafa
ddPCR

# ╔═╡ a336c29d-502c-4727-977b-578c9798ce53
mplasmid = mean(numuniqueplasmids)

# ╔═╡ 97919f1b-5428-44dc-928d-b9ee2f141a9e
begin
	hist = fit(Histogram, numuniqueplasmids, nbins=80)
	bar(hist.edges[1][1:end-1], hist.weights, legend=false, linecolor=:match, xlim=(-0.7, hist.edges[1][end]+2), ylim=(0, 1.05*maximum(hist.weights)), size=(600, 200), tickfontsize=16, xticks=[0, 20, 40, 60], yticks=[0, 800, 1600])
	vline!([mplasmid], line=(2, 4))
	vline!([ddPCR], line=(:black, 4))

	# Fill the region between a and b
	fill_between_x = [ddPCR, mplasmid, mplasmid, ddPCR]
	fill_between_y = [0, 0, 2*maximum(hist.weights), 2*maximum(hist.weights)]
	plot!(fill_between_x, fill_between_y, fill=(0, 0.4, 1))
	# savefig("./Plots/CopyHist_PreBin_$(condition)_$(replicate).png")
end

# ╔═╡ 8e6e6608-f389-4c9d-9730-e13fb91a2319
begin
	finalnumplasmids = vcat([sum(val != 0 for val in row) for row in eachrow(finalcelldf[:, 2:end])], zeros(Int, noBCcells))
	finalmplasmid = mean(finalnumplasmids)
end

# ╔═╡ 1454ac4e-5c4f-4f36-bbfb-d95574695e0d
begin
	histfinal = fit(Histogram, finalnumplasmids, nbins=60)
	bar(histfinal.edges[1][1:end-1], histfinal.weights, legend=false, linecolor=:match, xlim=(-0.7, 22), ylim=(0, 2200), size=(300, 200), tickfontsize=14, xticks=[0, 10, 20], yticks=[0, 1000, 2000])
	vline!([finalmplasmid], line=(2, 4))
	vline!([ddPCR], line=(:black, 4))

	# Fill the region between a and b
	fill_between_xfinal = [ddPCR, finalmplasmid, finalmplasmid, ddPCR]
	fill_between_yfinal = [0, 0, 2*maximum(histfinal.weights), 2*maximum(histfinal.weights)]
	plot!(fill_between_xfinal, fill_between_yfinal, fill=(0, 0.4, :blue))
	# savefig("./Plots/CopyHist_PostBin_$(condition)_$(replicate).png")
end

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

# ╔═╡ 7b028b95-79ee-4133-a986-64be7579c1df
begin
	TFlist = names(celldf)[2:end]
end

# ╔═╡ 127cd7d8-241d-40ff-a9f6-5e2e2f307d56
begin
	TFkept = []

	# Filter for genes with at least 200 counts greater than 10
	for TF in TFlist
		if sum(celldf[:, TF] .>= 10) > 200
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
	testTF = TFkept[15]
	testcountdict = formcountdict(testTF, celldf, 0)
	testmax = maximum(keys(testcountdict))
	testnumcounts = sum(values(testcountdict))
	testcounts, testcountfreq = formcountfreq(testcountdict)
end

# ╔═╡ 725a2c87-b9c3-418c-8fc1-48c448713400
plot(testcounts, testcountfreq, ylim=(0, 0.003), label=false, xlabel="Counts", ylabel="Frequency")

# ╔═╡ dc1c2c00-1a17-461e-8dfb-3ee6a68c9fb8
testTF

# ╔═╡ e310b53e-de2a-4295-9b69-6b7657e29951
plotTF = "Plasmid_16"

# ╔═╡ dce8b0ee-6bf9-41da-8096-30cf02139cbc
begin
	plotcountdict = formcountdict(plotTF, celldf, 0)
	plotmax = maximum(keys(plotcountdict))
	plotnumcounts = sum(values(plotcountdict))
	plotcounts, plotcountfreq = formcountfreq(plotcountdict)
end

# ╔═╡ 4f3b95e9-ca2a-41c8-a0e8-4a149642e4b8
begin
	plot(plotcounts, plotcountfreq, ylim=(0, 1), xlim=(0, plotmax+10), label=false, color = 1, lw=10, tickfontsize=16)
	plot!(size=(600, 300))
	# savefig("./Plots/AllCounts_R01_P16.png")
end

# ╔═╡ 41929758-8a1a-4e5d-8078-36e715fbb2a9
begin
	plot(plotcounts[1:6], plotcountfreq[1:6], fillrange=0, fillalpha=0.3, fillcolor=4, label=false)
	plot!(plotcounts, plotcountfreq, label=false, color = 1, lw=4, tickfontsize=16)  
	
	plot!(ylim=(0, 0.008), xlim=(0, 110),xticks = [0, 100], yticks=[0.0, 0.004, 0.008], size=(600, 300))
	# savefig("./Plots/CountsZoomed_R01_P16.png")
end

# ╔═╡ eb493270-6eec-4883-9237-c0b8862a332b
sum(plotcountfreq[2:6])

# ╔═╡ 1fced61e-e061-42e5-b153-0d64f7d3d4b8
begin
	ν, γ, ρμ, α, f, threshold = noisefitdf[:, plotTF]
	
	gex(z) = GNB(z, ρμ, α)
	gexν(z) = GNB(z, ρμ*ν, α)
	gcont(z) = Gpois(gexν(z), f*γ)
	gseq(z) = gcont(z)*(f*gex(z) + (1-f))

	pnoexlist = taylor_expand(z->gcont(z)*gex(z), order = plotcounts[end]).coeffs
    pseqlist = taylor_expand(gseq, order = plotcounts[end]).coeffs
    pnolistlong = taylor_expand(gcont, order = plotcounts[end]).coeffs
end

# ╔═╡ 68bff31e-b982-4183-8c5a-8585aee5f47d
ν

# ╔═╡ 3db6e14f-aca3-44ed-bd3c-12319468b1ed
begin
	plot(plotcounts, plotcountfreq, color = 1, lw=3, label="Counts")

	plot!(plotcounts, (1-f).*pnolistlong, lw=3, label="Contaminant")
	plot!(plotcounts, f.*pnoexlist, lw=3, label="Exp+Cont", c=5)
	
	# plot!(plotcounts, pseqlist.+4e-5, lw=3, label="Total")
	
	plot!(ylim=(0, 0.008), xlim=(0, 110),xticks = [0, 100], yticks=[0.0, 0.004, 0.008], size=(600, 300), tickfontsize=16, legendfontsize=12)
	savefig("./Plots/Fit_Pre_R01_P16.png")
end

# ╔═╡ 31d069e8-3838-4686-9f05-8ae423f70c47
# begin
# 	testlower = [0., 0., 0., 0., 0.]
# 	testupper = [1., 1e4, 1e4, 10, 1.]
# 	testthresh = 10
	
# 	testres = fitgene_MLE(testcountdict, testlower, testupper, MoMinitial(testcountdict, testthresh))
# end

# ╔═╡ 426cd89c-25de-4c59-98cd-a0a976adf186
# begin
# 	νi, γi, μi, ri, fi = MoMinitial(testcountdict, testthresh)
# 	gnoi(z) = Gnoise(z, νi, fi*γi, μi, ri)
# 	gexi(z) = GNB(z, μi, ri)
#     gnoexi(z) = gexi(z)*gnoi(z)
#     gseqi(z) = fi*gnoexi(z) + (1-fi)*gnoi(z)
	
# 	ν, γ, μ, r, f = testres.minimizer

# 	gno(z) = Gnoise(z, ν, f*γ, μ, r)
# 	gex(z) = GNB(z, μ, r)
#     gnoex(z) = gex(z)*gno(z)
#     gseq(z) = f*gnoex(z) + (1-f)*gno(z)

	
	
# 	pnoexlist=taylor_expand(z->f*gnoex(z), order=testmax).coeffs
# 	peak = max(maximum(pnoexlist), 0.0015/4)

# 	plot(testcounts, testcountfreq, ylim=(0, 0.003), xlim=(0, 0.75*testmax), label="Counts", xlabel="$(testTF) Count\nν=$(@sprintf("%.3g", ν)) , γ=$(@sprintf("%.1f", γ)) , μ=$(@sprintf("%.1f", μ)) , r=$(@sprintf("%.3g", r)), f=$(@sprintf("%.3g", f))", ylabel="Fraction of Droplets")
	
# 	plot!(0:testmax, taylor_expand(z->(1-f)*gno(z), order=testmax).coeffs, label="Noise", lw=2)

# 	plot!(0:testmax, pnoexlist, label="Expression", lw=2)
	
# 	plot!(0:testmax, taylor_expand(z->gseq(z), order=testmax).coeffs .+ peak*2e-2, label="Mixture", lw=2)
	
# 	# plot!(0:testmax, taylor_expand(z->Gseq(z, MoMinitial(testcountdict, testthresh)...), order=testmax).coeffs, label="Initial", lw=2)
# end

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
# ╠═7dfdca3e-ad92-4c15-8709-fa32b3496197
# ╠═50580bae-24ff-40c8-98f9-5600589572fc
# ╠═86042e1e-3ed1-4fd4-8296-8f25c61c3f16
# ╠═96f7230d-1f62-45f8-91e8-0b4f7dd28348
# ╠═4778354d-4208-4b65-ba8f-67f5f86233b9
# ╠═5961d6a2-4b13-4865-b6be-0e8ece20d5f2
# ╠═e0a3512d-2b66-4f5c-9bda-8c7a0712ecbf
# ╠═0c15417f-d723-4672-8039-c33564b12f61
# ╠═b501569c-59e7-46d4-a059-02708c44eb99
# ╠═07dd1ce6-7791-4b37-9695-24ecb8dd6fcb
# ╠═35134891-e3dd-459c-8703-97235be886ea
# ╠═b09d1534-311a-4fd2-83d2-d1a9628b0c34
# ╠═ff72b3c2-c4e1-458a-9f57-34fbcd694327
# ╠═487249ca-ddcf-4dbd-8c79-dde5639cfb29
# ╟─7c18dc4c-fab0-48ba-8693-1e0158cb029e
# ╟─88e75c2b-7b0f-4625-979c-6b48bd034d06
# ╠═78ca2110-df58-4032-a5f5-21d6a9fcae59
# ╟─c9b28149-b848-4cae-8de7-70533628c7dd
# ╠═b98ae866-3d06-4126-b58d-537eaab6b504
# ╠═d35b5774-cb3d-4214-9a70-48dbad96ca84
# ╠═054d8a56-a86d-4883-be0e-901b01893570
# ╟─7c28ffc9-3720-4808-ae17-cabad2d89c68
# ╟─c1e1a053-6d49-4ae8-a1d1-c6fa5f712c3c
# ╟─60f09855-a237-4259-8c6d-98eaf946c40d
# ╠═09fbc6b9-94e4-4109-9173-571bf648f35a
# ╠═1cc1ee57-3d3b-48aa-941e-8bf70b3e3132
# ╠═037c7154-c014-476a-9b36-e2605884be04
# ╠═7345e9fd-438e-4117-89a8-b7a88cb70dc4
# ╠═e101d9ac-1b7d-4cb3-896a-edc2b9ecbcde
# ╠═9441c065-2d92-43ad-964a-51f2838971c2
# ╠═53c56589-b43e-4e23-9bfb-6b1097715ad0
# ╠═482c62fa-12f5-48fd-bdb7-9ce62bceaabc
# ╠═1614673f-5bc8-4907-a89c-95d826dbedf7
# ╠═777d7c02-a033-4060-a11d-3b70278edafa
# ╠═a336c29d-502c-4727-977b-578c9798ce53
# ╠═97919f1b-5428-44dc-928d-b9ee2f141a9e
# ╠═8e6e6608-f389-4c9d-9730-e13fb91a2319
# ╠═1454ac4e-5c4f-4f36-bbfb-d95574695e0d
# ╟─6377a680-8d8e-4675-858f-427bb75db022
# ╠═70400124-5b17-4ad6-bddb-4a2715d0d7d2
# ╠═7b028b95-79ee-4133-a986-64be7579c1df
# ╠═127cd7d8-241d-40ff-a9f6-5e2e2f307d56
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
# ╠═68bff31e-b982-4183-8c5a-8585aee5f47d
# ╠═3db6e14f-aca3-44ed-bd3c-12319468b1ed
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
