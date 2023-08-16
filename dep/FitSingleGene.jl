using Pkg
Pkg.activate("Project.toml")

using DataFrames, CSV
using TaylorSeries, Optim

include("../FitUtilities.jl")

#######################################################
# Set paths and parameters for fit
#######################################################
replist = readlines("replicates.txt")
condlist = readlines("condition.txt")
datapathprefix = readline("pathinfo.txt")

rep = ARGS[1]
cond = ARGS[2]
contdist = ARGS[3]
initialthresh = parse(Int64, ARGS[4])
BC = ARGS[5]


datapath = joinpath(datapathprefix, "$(cond)_"*join(replist))
outputpath = joinpath(datapath, "GeneFitData", "allfits")
mkpath(outputpath)

countpath = joinpath(datapath, "$(cond)_$(rep)_CellCounts.tsv")
parampath = joinpath(datapath, "GeneFitData", "contaminantfits", "$(cond)_$(rep)_ContaminantParameters_$(contdist).txt")

#######################################################
# Load Counts
#######################################################

println("Loading data for $(cond) $(rep) $(BC)")
allcountsdf = DataFrame(CSV.File(countpath, delim='\t'))
noBCcells = parse(Int, readline(joinpath(datapath, "$(cond)_$(rep)_noBCcells.txt")))

BCcountdict = countmap(allcountsdf[!, BC])
BCcountdict[0] += noBCcells

contparam = [parse(Float64, s) for s in readlines(parampath)]

#######################################################
# Run Fit
#######################################################
if contdist == "Full"
    γ_const, γ_MAD, νγ_const, νγ_MAD = contparam

    νi, γi, ρμi, αi, fi = MoMFull(BCcountdict, initialthresh)
    
    # [ν, γ, μ, α, f]
    lower = [0., 0., 0., 0., 0.]
    upper = [1., 1e3, 1e3, 1e2, 1]
    initial = [νγ_const/γ_const, γ_const, min(max(ρμi, 20), 100.), min(max(αi, 0.1), 5.), min(max(fi, 0.05), 0.9)]

    fitresults = fitgene_MLE_Full_reg(BCcountdict, contparam, lower, upper, initial)

    print("\tFit complete. ")

elseif contdist == "Poisson"
    νγ_const, νγ_MAD = contparam

    λi, ρμi, αi, fi = MoMPoisson(BCcountdict, initialthresh)
    # [λ, ρμ, α, f]
    lower = [0., 0., 0., 0.]
    upper = [1e2, 1e3, 1e2, 1.]
    initial = [νγ_const, min(max(ρμi, 20), 100.), min(max(αi, 0.1), 5.), min(max(fi, 0.05), 0.9)]

    fitresults = fitgene_MLE_Poisson_reg(BCcountdict, contparam, lower, upper, initial)
    
    print("Fit complete. ")

end
    
#######################################################
# Save Results
#######################################################
gex, gcont, gseq = genfuncs(fitresults.minimizer...)
f = fitresults.minimizer[end]

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

parameterdf[!, BC] = [fitresults.minimizer..., threshold]
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


