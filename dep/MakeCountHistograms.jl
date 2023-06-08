### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ a5881b94-d209-11ed-12d2-d11029fabab3
begin
	using Pkg
	Pkg.activate("Project.toml")
end

# ╔═╡ 7943e527-77cb-4366-b172-d0b9f785068f
begin
	using StatsBase
	using CSV, DataFrames
end

# ╔═╡ db15663b-17df-4c10-bd37-43ff14acdc54
begin
	replist = ["R11", "R12", "R13"]
	markerlist = ["NG", "FN", "CD45", "NKp46"]
end

# ╔═╡ f8275dc2-2c0e-49f4-abcb-85100a13897c
begin
	rep = replist[3]
	marker = markerlist[4]
end

# ╔═╡ 950c9228-dc30-4c0f-a59a-6851f68b4a0c
barcodedf = DataFrame(CSV.File("../RandomDesigns/ScreenData/NKv/bb20230124104525_NKv_R11R12R13_counts/bb20230124104525_NKv_$(rep)_$(marker)_eTF_count.txt", delim="\t"))

# ╔═╡ 19a3c813-ff5c-4ee7-8329-9ce3658392eb
outputpath = "../RandomDesigns/ScreenData/NKv/bb20230124104525_NKv_counthistograms/"

# ╔═╡ 517dcff9-3459-4658-bc64-91f4318d8543
gexdf = DataFrame(CSV.File("../RandomDesigns/ScreenData/NKv/bb20230124104525_NKv_R11R12R13_GEX_dropletlabels/bb20230124103838_NKv_$(rep)_$(marker)_GEX_prefiltering_table_cellsandnocells.txt", delim="\t"))

# ╔═╡ 2e4acf7f-1f2e-486c-801b-78886d4b4a19
# Pull cells from gex file
cellidentifiers = gexdf[gexdf[:, end] .== true, "orig.ident"]

# ╔═╡ e2ce217b-0d5c-4a3f-a015-47df08694d7b
cellmask = [c ∈ cellidentifiers for c in barcodedf[:,:cell_barcode]]

# ╔═╡ 57e37940-18fe-43ff-ade2-47985cf9ff40
zerobcdroplets = 0

# ╔═╡ 3f576a80-6fe0-4f9d-a4c9-1e21faca7eb8
begin
	cellBCdf = DataFrame()
	notcellBCdf = DataFrame()
	
	append!(cellBCdf, barcodedf[cellmask, :])
	append!(notcellBCdf, barcodedf[.!cellmask, :])
end

# ╔═╡ 93431543-50d2-485e-9edf-20d5f2581220
begin
	function formcountdict(BC, BCdf, zerobackground)

		genecounts = BCdf[:, BC]
		countdict = countmap(genecounts)
		countdict[0] += zerobackground
		
		return countdict
	end
end

# ╔═╡ b2123a35-0241-44a7-9e91-a1fe4838a6ce
begin
	# Make histograms from countmaps and find the maximum count values needed
	BClist = names(cellBCdf)[2:end]

	cellBCdict = Dict()
	emptydropletBCdict = Dict()
	
	for BC in BClist
		cellhist = formcountdict(BC, cellBCdf, 0)
		droplethist = formcountdict(BC, notcellBCdf, zerobcdroplets)

		cellBCdict[BC] = cellhist
		emptydropletBCdict[BC] = droplethist
	end

	cellmaxcount = maximum([maximum(keys(BCdict)) for BCdict in values(cellBCdict)])
	dropletmaxcount = maximum([maximum(keys(BCdict)) for BCdict in values(emptydropletBCdict)])
end

# ╔═╡ 2902076f-454a-4361-9734-f68519340f7b
begin
	# Place histograms into a dataframe as columns
	cellBChistdf = DataFrame()
	emptydropletBChistdf = DataFrame()
	
	for BC in BClist
		cellBChistdf[!, BC] = [haskey(cellBCdict[BC], c) ? cellBCdict[BC][c] : 0 for c in 0:cellmaxcount]
		emptydropletBChistdf[!, BC] = [haskey(emptydropletBCdict[BC], c) ? emptydropletBCdict[BC][c] : 0 for c in 0:dropletmaxcount]
	end
end

# ╔═╡ bcb8b979-d34f-440c-a2d1-986a36f8ddcb
begin
	mkpath(outputpath)
	CSV.write(outputpath*"bb20230124104525_NKv_$(rep)_$(marker)_cell_count_histograms.tsv",cellBChistdf, delim='\t')
	CSV.write(outputpath*"bb20230124104525_NKv_$(rep)_$(marker)_emptydroplet_count_histograms.tsv",emptydropletBChistdf, delim='\t')
end

# ╔═╡ Cell order:
# ╠═a5881b94-d209-11ed-12d2-d11029fabab3
# ╠═7943e527-77cb-4366-b172-d0b9f785068f
# ╠═db15663b-17df-4c10-bd37-43ff14acdc54
# ╠═f8275dc2-2c0e-49f4-abcb-85100a13897c
# ╠═950c9228-dc30-4c0f-a59a-6851f68b4a0c
# ╠═19a3c813-ff5c-4ee7-8329-9ce3658392eb
# ╠═517dcff9-3459-4658-bc64-91f4318d8543
# ╠═2e4acf7f-1f2e-486c-801b-78886d4b4a19
# ╠═e2ce217b-0d5c-4a3f-a015-47df08694d7b
# ╠═57e37940-18fe-43ff-ade2-47985cf9ff40
# ╠═3f576a80-6fe0-4f9d-a4c9-1e21faca7eb8
# ╠═93431543-50d2-485e-9edf-20d5f2581220
# ╠═b2123a35-0241-44a7-9e91-a1fe4838a6ce
# ╠═2902076f-454a-4361-9734-f68519340f7b
# ╠═bcb8b979-d34f-440c-a2d1-986a36f8ddcb
