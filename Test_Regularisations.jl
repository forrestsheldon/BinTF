### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 7878da9c-fe22-11ed-35ff-afa72bda41b4
begin
	using Pkg
	Pkg.activate("Project.toml")
end

# ╔═╡ 3be08f88-58e2-47e1-80d4-70bf3e31f824
begin
	using DataFrames, CSV
	using TaylorSeries, Optim
	using Plots
	
	include("fitutils.jl")
end

# ╔═╡ e699d95c-ad8f-4c06-840c-a70a8861b47a
begin
	#######################################################
	# Set paths and parameters for screen
	#######################################################
	replist = readlines("replicates.txt")
	condlist = readlines("condition.txt")
	datapathprefix = readline("pathinfo.txt")
	
	rep = "R01"
	cond = "PreAmp"
	
	datapath = joinpath(datapathprefix, "$(cond)_"*join(replist), "GeneFitData")
	outputpath = joinpath(datapath, "allfits")
	
	cellhistpath = joinpath(datapath, "counthistograms", "$(cond)_$(rep)_cell_count_histograms.tsv")
end

# ╔═╡ 3894a63c-a2c6-4fbd-b6ce-ed959e9266c9
begin
	#######################################################
	# Parameters for run
	#######################################################
	
	# to fit the global noise parameters keep only the genes with at least
	# excludedcounts counts >= excludedthresh
	excludedthresh = 6
	excludedcounts = 30
	
	# For the initial conditions of the fit, use this threshold to separate noise
	# below and expression above
	highthresh = 10
	lowthresh = 2
end

# ╔═╡ f6f1cd10-aa02-4b6f-895c-32ec22fa0803
begin
	#######################################################
	# Data Processing Loops
	#######################################################
	runident = "$(cond)_$(rep)"
	# Load the Noise Constant set by the last fit
	ν_const, γ_const = [parse(Float64, s) for s in readlines(joinpath(datapath, "noisefits", runident*"_NoiseParameters.txt"))]
	νγ_const = ν_const*γ_const
	
	########
	# Alternative Noise parameters
	regparam = [parse(Float64, s) for s in readlines(joinpath(datapath, "noisefits", runident*"_NoiseParameters2.txt"))]
	
	# Load Data
	println("Loading data for $(cond) $(rep)")
	cellhistogramdf = DataFrame(CSV.File(cellhistpath))
	
	BClist = names(cellhistogramdf)
end

# ╔═╡ 51a1c8ad-2c05-4e18-9303-3f0010de15a1
begin
	# Move dataframe vectors to dictionaries
	cellcountdicts = Dict{String, Dict{Int, Int}}()
	
	# Form count dictionaries
	for BC in BClist
	    cellcountdicts[BC] = Dict{Int, Int}()
	    
	    cellcounts = cellhistogramdf[!, BC]
	    
	    for (idx, numcells) in enumerate(cellcounts)
	        count = idx-1
	        if numcells != 0
	            cellcountdicts[BC][count] = numcells
	        end
	    end
	end
end

# ╔═╡ bfd165fc-bed5-4b9d-92d8-4e652ce390bf
begin
	BC = BClist[5]
	celldict = cellcountdicts[BC]
	
	if sum(nc for (c, nc) in celldict if c >= excludedthresh; init=0) > excludedcounts
        momthresh = highthresh
    else
        momthresh = lowthresh
    end

    νi, γi, μi, ri, fi = MoMinitial(celldict, momthresh)
    
    
    totalcounts = sum(values(celldict))
    try
        global fmax = min(1., (totalcounts - celldict[0] + celldict[1]^2/celldict[2])/totalcounts)
    catch KeyError
        try
            global fmax = min(1., (totalcounts - celldict[0] + celldict[1])/totalcounts)
        catch KeyError
            global fmax = 1
        end
    end

	lower = [0., 0., 0., 0., 0.]
	upper = [1., 1e3, 1e3, 10, fmax]
	initial = [ν_const, γ_const, μi, ri, fi]
end

# ╔═╡ bcb24809-8e51-444d-920a-4c165fd33b44
begin
	maxcount = maximum(keys(celldict))
	sumcounts = sum(nc for (c, nc) in celldict)
	νγ_med, νγ_MAD, hprbl_med, hprbl_MAD = regparam

	function objective(param)
        ν, γ, μ, r, f = param
        ν = max(1e-20, ν)
        γ = max(1e-20, γ)
        gν(z) = Gpois(z, ν)
        gex(z) = GNB(z, μ, r)
        gno(z) = Gpois(gex(gν(z)), f*γ)
        gnoex(z) = gex(z)*gno(z)
        gseq(z) = f*gnoex(z) + (1-f)*gno(z)

        plist = taylor_expand(gseq, order = maxcount).coeffs
        hprbl = asinh((ν-γ)/(2*sqrt(ν*γ)))
        
        sum(-nc*log(plist[c+1]) for (c, nc) in celldict) / sumcounts+ (ν*γ - νγ_med)^2/(2*(1.48*νγ_MAD)^2)+ (hprbl - hprbl_med)^2/(2*(1.48*hprbl_MAD)^2)

    end

	objective(initial)
end

# ╔═╡ a8e79472-807c-4959-a9f1-593c7926c6d0
Optim.optimize(objective, lower, upper, initial)

# ╔═╡ 336442e1-a717-4e8b-9da6-e9aec532dc47
begin
	fitresults = fitgene_MLE_reg(celldict, regparam, lower, upper, initial)
	ν, γ, μ, r, f = fitresults.minimizer
end

# ╔═╡ 8b3b50f5-57b3-4a16-9bc3-a55eeefc3257
begin
	# mixture fits
	cellcountvec, cellcountfreq = formcountfreq(celldict)
	gnof(z) = Gnoise(z, ν, f*γ, μ, r)
	gexf(z) = GNB(z, μ, r)
	gnoexf(z) = gexf(z)*gnof(z)
	gseqf(z) = f*gnoexf(z) + (1-f)*gnof(z)
	pnoexflist = taylor_expand(gnoexf, order = cellcountvec[end]).coeffs
	pseqflist = taylor_expand(gseqf, order = cellcountvec[end]).coeffs
	pnoflistlong = taylor_expand(gnof, order = cellcountvec[end]).coeffs
	
	threshold = cellcountvec[findfirst(f.*pnoexflist[1:length(cellcountvec)] .> (1-f).*pnoflistlong[1:length(cellcountvec)])]
	
end

# ╔═╡ 6f2eacee-50ad-44ec-bae2-e5d361b73b3f
begin
	plot(cellcountvec, cellcountfreq, label = "Counts",ylim=(0, 0.001))
	plot!(cellcountvec, pseqflist, label = "Pseq")
	plot!(cellcountvec, f.*pnoexflist, label = "Pex")
	plot!(cellcountvec, (1-f).*pnoflistlong, label="Pno")
end

# ╔═╡ 58961776-8846-409d-9f2e-add3bda46e42
begin
		function objective_parts(param)
	        ν, γ, μ, r, f = param
	        ν = max(1e-20, ν)
	        γ = max(1e-20, γ)
	        gν(z) = Gpois(z, ν)
	        gex(z) = GNB(z, μ, r)
	        gno(z) = Gpois(gex(gν(z)), f*γ)
	        gnoex(z) = gex(z)*gno(z)
	        gseq(z) = f*gnoex(z) + (1-f)*gno(z)
	
	        plist = taylor_expand(gseq, order = maxcount).coeffs
	        hprbl = asinh((ν-γ)/(2*sqrt(ν*γ)))
	        
	        sum(-nc*log(plist[c+1]) for (c, nc) in celldict) / sumcounts, (ν*γ - νγ_med)^2/(2*(1.48*νγ_MAD)^2), (hprbl - hprbl_med)^2/(2*(1.48*hprbl_MAD)^2)
	
	    end
	objective_parts(fitresults.minimizer)
end

# ╔═╡ Cell order:
# ╠═7878da9c-fe22-11ed-35ff-afa72bda41b4
# ╠═3be08f88-58e2-47e1-80d4-70bf3e31f824
# ╠═e699d95c-ad8f-4c06-840c-a70a8861b47a
# ╠═3894a63c-a2c6-4fbd-b6ce-ed959e9266c9
# ╠═f6f1cd10-aa02-4b6f-895c-32ec22fa0803
# ╠═51a1c8ad-2c05-4e18-9303-3f0010de15a1
# ╠═bfd165fc-bed5-4b9d-92d8-4e652ce390bf
# ╠═bcb24809-8e51-444d-920a-4c165fd33b44
# ╠═a8e79472-807c-4959-a9f1-593c7926c6d0
# ╠═336442e1-a717-4e8b-9da6-e9aec532dc47
# ╠═8b3b50f5-57b3-4a16-9bc3-a55eeefc3257
# ╠═6f2eacee-50ad-44ec-bae2-e5d361b73b3f
# ╠═58961776-8846-409d-9f2e-add3bda46e42
