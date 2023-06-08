### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 98ee3970-a64d-11ed-2c51-7bb5dda97592
begin
	using Pkg
	Pkg.activate("Project.toml")
end

# ╔═╡ 215b45f2-fd7d-4e0e-87ed-921ae04d059c
begin
	using StatsBase
	using CSV, DataFrames
	using Plots, Plots.PlotMeasures, StatsPlots
end

# ╔═╡ 0aea22d8-2ecd-4f7e-9484-218177f2a30c
celldf = DataFrame(CSV.File("GB_Cell_TF_counts.csv", delim="\t"))

# ╔═╡ 1481ece7-bfef-4f83-a83b-0a364a8113b7
notcelldf = DataFrame(CSV.File("GB_NotCell_TF_counts.csv", delim="\t"))

# ╔═╡ 111c3fb2-7dbd-4c99-a449-aaf1bfe88010
zerobackground = parse(Int, readline("GB_Zerobackground.txt"))

# ╔═╡ b78a7fbd-3352-4472-997e-459c2cbee083
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

# ╔═╡ e987b55c-f5b9-4af2-a188-529dbd8ee197
begin
	TFlist = names(celldf)[2:end]
	
	TFidx = 9

	cellcountdict = formcountdict(TFlist[TFidx], celldf, 0)
	notcellcountdict = formcountdict(TFlist[TFidx], notcelldf, zerobackground)

	cellcounts, cellfreq = formcountfreq(cellcountdict)
	notcellcounts, notcellfreq = formcountfreq(notcellcountdict)

end

# ╔═╡ 0c434835-ffa5-4d0e-a165-b75b93381f16
begin
	plot(cellcounts, cellfreq, ylim = (0, 1), xlim=(0, cellcounts[end]), linewidth=10, label=false, tickfontsize=18)# xlabel="Count", ylabel="Fraction of Droplets")
	# plot(counts, countfreq, ylim = (0, 1), xlim=(0, counts[end]), linewidth=2, label=false, xlabel="Count", ylabel="Fraction of Droplets")
	

	plot!(cellcounts, cellfreq, ylim = (0, 0.003), xlim = (0, 0.6*cellcounts[end]), linewidth = 3, inset=(1, bbox(0.2, 0.23, 0.75, 0.75, :bottom, :left)), subplot = 2,label="Droplets Containing Cells", thickness_scaling = 1, xticks = [0, 100, 200], yticks=[0.000, 0.002,], tickfontsize = 18)
	
	plot!(notcellcounts, notcellfreq, subplot = 2, linewidth = 3, label="Empty Droplets", thickness_scaling = 1, legendfontsize=18, size=(900, 750))

	# savefig("Plots/Rawcounts_Inset_Plasmid_9.png")
end

# ╔═╡ d1e47ab7-8ae7-4146-8bf1-f1b789d7cc21
begin
	# Splitting up the previous plot
	plot(cellcounts, cellfreq, ylim = (0, 1), xlim=(0, cellcounts[end]), linewidth=10, label=false, tickfontsize=12, size=(500, 300))

	# savefig("Plots/Rawcounts_fullaxes_Plasmid_9.png")
end

# ╔═╡ 94c8b828-8f06-4100-b6e3-49e5a8ca908b
begin
	plot(cellcounts, cellfreq, ylim = (0, 0.003), xlim = (0, 0.6*cellcounts[end]), linewidth = 3, label="Droplets Containing Cells", thickness_scaling = 1, xticks = [0, 100, 200], yticks=[0.000, 0.003,], tickfontsize = 12)
	
	plot!(notcellcounts, notcellfreq, linewidth = 3, label="Empty Droplets", thickness_scaling = 1, legendfontsize=12, size=(500, 300))

	# savefig("Plots/Rawcounts_zoomedaxes_Plasmid_9.png")
end

# ╔═╡ 2a857da5-43f6-4b57-b029-10caec0127a0
begin
	plot(cellcounts, cellfreq, ylim = (0, 0.0015), xlim = (0, 0.7*cellcounts[end]), linewidth=2, label=false) #xlabel="Count", ylabel="Fraction of Droplets")
	
	plot!(notcellcounts, notcellfreq, linewidth = 2, label="Droplets Containing Cells", thickness_scaling = 1, legend=false, size=(500, 300), tickfontsize=12)

	# savefig("Plots/Rawcounts_Plasmid_$(TFidx).png")
end

# ╔═╡ a44c1a44-0829-4f2f-97d4-3081b43b3299
begin
	plot(cellcounts, cellfreq, ylim = (0, 0.005), xlim = (0, 20), linewidth=2, label=false, size=(500, 300), tickfontsize=12) #xlabel="Count", ylabel="Fraction of Droplets")
	
	# plot!(notcellcounts, notcellfreq, linewidth = 2, label="Droplets Containing Cells", thickness_scaling = 1, legend=false, size=(500, 300), tickfontsize=12)

	# savefig("Plots/Rawcounts_Noise_Plasmid_$(TFidx).png")
end

# ╔═╡ ff50ae47-27aa-46e4-819a-66fdd0d28c40
begin
	notcellmean = Vector{Float64}()
	cellmean = Vector{Float64}()
	
	for TF in TFlist
	
		notcellcountdict = formcountdict(TF, notcelldf, zerobackground)
	
		push!(notcellmean, sum(c*nc for (c, nc) in notcellcountdict) / sum(values(notcellcountdict)))
		
		cellcountdict = formcountdict(TF, celldf, 0)
	
		push!(cellmean, sum(c*nc for (c, nc) in cellcountdict) / sum(values(cellcountdict)))
		
	end
	
	scatter(cellmean, notcellmean, legend=false, xlim = (0, 12), ylim = (0, 0.11), xticks=[0, 4, 8, 12], yticks=[0., 0.05, 0.1], tickfontsize=12, markerstrokewidth=0)# ylabel = "Mean Count for Empty Droplets", xlabel = "Mean Count for Droplets Containing Cells")
	plot!([0, 12], [0, (cellmean \ notcellmean)*12], lw=2, size=(500, 300))
	# savefig("./Plots/BGCellCorrelation.png")
end

# ╔═╡ Cell order:
# ╠═98ee3970-a64d-11ed-2c51-7bb5dda97592
# ╠═215b45f2-fd7d-4e0e-87ed-921ae04d059c
# ╠═0aea22d8-2ecd-4f7e-9484-218177f2a30c
# ╠═1481ece7-bfef-4f83-a83b-0a364a8113b7
# ╠═111c3fb2-7dbd-4c99-a449-aaf1bfe88010
# ╠═b78a7fbd-3352-4472-997e-459c2cbee083
# ╠═e987b55c-f5b9-4af2-a188-529dbd8ee197
# ╠═0c434835-ffa5-4d0e-a165-b75b93381f16
# ╠═d1e47ab7-8ae7-4146-8bf1-f1b789d7cc21
# ╠═94c8b828-8f06-4100-b6e3-49e5a8ca908b
# ╠═2a857da5-43f6-4b57-b029-10caec0127a0
# ╠═a44c1a44-0829-4f2f-97d4-3081b43b3299
# ╠═ff50ae47-27aa-46e4-819a-66fdd0d28c40
