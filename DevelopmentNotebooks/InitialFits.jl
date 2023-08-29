### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 3612abe4-3052-11ee-1dc3-f967a4c493e1
begin
	using Pkg
	Pkg.activate("../Project.toml")
end

# ╔═╡ 774cb3e9-8396-4644-bc0b-10f91933434d
begin
	using StatsBase
	using TaylorSeries
	using CSV, DataFrames, JSON
	using Plots, Plots.PlotMeasures, StatsPlots
	using Printf
	using Colors
end

# ╔═╡ c198432e-fd6c-4af2-811e-383817bbbecf
include("../FitUtilities.jl")

# ╔═╡ 72d5a56b-1251-413b-be5e-70e9f39d5e61
begin
	cond = "withAmp"
	datapath = joinpath("../Data", "$(cond)_R01R02")
	rep = "R01"

	contfitdf = DataFrame(CSV.File(joinpath(datapath, "GeneFitData", "contaminantfits","$(cond)_$(rep)_Parameters_Full.tsv"), delim="\t"))
	selected = JSON.parse(readline(joinpath(datapath, "GeneFitData", "contaminantfits", "$(cond)_$(rep)_SelectedFits_Full.json")))

	# wo2noisefitdf = DataFrame(CSV.File(joinpath(withoutAmpPath, "GeneFitData", "contaminantfits","withoutAmp_R02_Parameters_Full.tsv"), delim="\t"))
	# w2noisefitdf = DataFrame(CSV.File(joinpath(withAmpPath, "GeneFitData", "contaminantfits","withAmp_R02_Parameters_Full.tsv"), delim="\t"))

end

# ╔═╡ e1281b06-7198-4e27-802a-0fceb5cb49ac
begin
	paramdf = contfitdf[:, selected]
	νvec = Vector(paramdf[1,:])
	γvec = Vector(paramdf[2,:])
	ρμvec = Vector(paramdf[3,:])
	αvec = Vector(paramdf[4,:])
	fvec = Vector(paramdf[5,:])
end

# ╔═╡ 9a12dc93-31b7-456d-92b4-714f649c360b
histogram(αvec, nbins= 15)

# ╔═╡ ba11d865-17fe-4c89-b971-d0eb83286c7e
histogram(fvec)

# ╔═╡ 2a1a1fb1-1b80-40aa-b728-d3ba04bc5e8a
mean(fvec)

# ╔═╡ 66e23f36-cdab-4f7a-996e-45c92dda88e4
begin
	countdf = DataFrame(CSV.File(joinpath(datapath, "$(cond)_$(rep)_CellCounts.tsv"), delim="\t"))
	noBCcells = parse(Int, readline(joinpath(datapath, "$(cond)_$(rep)_noBCcells.txt")))
end

# ╔═╡ f38a89e8-f7d3-4a20-9d68-100be7d1bb92
begin
	meandict = Dict(plasmid=>sum(countdf[:, plasmid])/(size(countdf, 1)+noBCcells) for plasmid in selected)
end

# ╔═╡ 9e874ed0-17af-42fe-ab38-3b3187511a30
scatter([meandict[p] for p in selected], fvec.*ρμvec)

# ╔═╡ Cell order:
# ╠═3612abe4-3052-11ee-1dc3-f967a4c493e1
# ╠═774cb3e9-8396-4644-bc0b-10f91933434d
# ╠═c198432e-fd6c-4af2-811e-383817bbbecf
# ╠═72d5a56b-1251-413b-be5e-70e9f39d5e61
# ╠═e1281b06-7198-4e27-802a-0fceb5cb49ac
# ╠═9a12dc93-31b7-456d-92b4-714f649c360b
# ╠═ba11d865-17fe-4c89-b971-d0eb83286c7e
# ╠═2a1a1fb1-1b80-40aa-b728-d3ba04bc5e8a
# ╠═66e23f36-cdab-4f7a-996e-45c92dda88e4
# ╠═f38a89e8-f7d3-4a20-9d68-100be7d1bb92
# ╠═9e874ed0-17af-42fe-ab38-3b3187511a30
