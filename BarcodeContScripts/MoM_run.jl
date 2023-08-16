using Pkg
Pkg.activate("Project.toml")

using DataFrames, CSV
using TaylorSeries

include("../FitUtilities.jl")

#######################################################
# Set paths and parameters for screen
#######################################################
replist = readlines("replicates.txt")
condlist = readlines("condition.txt")
datapathprefix = readline("pathinfo.txt")

rep = ARGS[1]
cond = ARGS[2]
contdist = ARGS[3]
thresh = parse(Int64, ARGS[4])
BC = ARGS[5]

datapath = joinpath(datapathprefix, "$(cond)_"*join(replist))
outputpath = joinpath(datapath, "GeneFitData", "allfits")
countpath = joinpath(datapath, "$(cond)_$(rep)_CellCounts.tsv")

contparam = [parse(Float64, s) for s in readlines(joinpath(datapath, "GeneFitData", "contaminantfits", "$(cond)_$(rep)_FixedParameters_$(contdist).txt"))]

#######################################################
# Load Counts
#######################################################

println("Loading data for $(cond) $(rep) $(BC)")
allcountsdf = DataFrame(CSV.File(countpath, delim='\t'))
noBCcells = parse(Int, readline(joinpath(datapath, "$(cond)_$(rep)_noBCcells.txt")))

BCcountdict = countmap(allcountsdf[!, BC])
BCcountdict[0] += noBCcells

#######################################################
# Run Fits
#######################################################
if contdist == "Full"
    ν, γ, α = contparam

    ρμ, f = MoMFullfμ(BCcountdict, ν, γ, α)

    fitparam = [ν, γ, ρμ, α, f]

elseif contdist == "Poisson"
    λ, α = contparam

    ρμ, f = MoMPoissonfμ(BCcountdict, λ, α)

    fitparam = [λ, ρμ, α, f]
end

gex, gcont, gseq = genfuncs(fitparam...)
f = fitparam[end]

cellcountvec, cellcountfreq = formcountfreq(BCcountdict)
pnoexlist = taylor_expand(z->gcont(z)*gex(z), order = cellcountvec[end]).coeffs
pseqlist = taylor_expand(gseq, order = cellcountvec[end]).coeffs
pnolistlong = taylor_expand(gcont, order = cellcountvec[end]).coeffs

threshold = cellcountvec[findfirst(f.*pnoexlist[1:length(cellcountvec)] .> (1-f).*pnolistlong[1:length(cellcountvec)])]	

#######################################################
# Save Results
#######################################################
# Check if parameter file exists and load if so
outputparamfilepath = joinpath(outputpath, "$(cond)_$(rep)_Parameters_$(contdist).tsv")

if isfile(outputparamfilepath)
    parameterdf = DataFrame(CSV.File(outputparamfilepath, delim='\t'))
else
    parameterdf=DataFrame()
end

parameterdf[!, BC] = [fitparam..., threshold]
CSV.write(outputparamfilepath, parameterdf, delim='\t')
   
println("\tGenerating Plots.")
# for plotting, save each barcode fit in different dataframes
mixturefitdf = DataFrame()
mixturefitdf[!, BC] = cellcountfreq
mixturefitdf[!, "$(BC)_contaminant"] = (1-f).*pnolistlong
mixturefitdf[!, "$(BC)_expression"] = f.*pnoexlist
mixturefitdf[!, "$(BC)_mixture"] = pseqlist

CSV.write(joinpath(outputpath, "$(cond)_$(rep)_MixtureFit_$(BC)_$(contdist).tsv"), mixturefitdf, delim='\t')
println("\tOutputs saved.")



