using Pkg
Pkg.activate("Project.toml")

using DataFrames, CSV
using TaylorSeries, Optim

include("../FitUtilities.jl")

#######################################################
# Set parameters for screen
#######################################################
replist = readlines("replicates.txt")
condlist = readlines("condition.txt")
datapathprefix = readline("pathinfo.txt")

rep = ARGS[1]
cond = ARGS[2]
contdist = ARGS[3]

#######################################################
# Parameters for run
#######################################################
# to fit the global noise parameters keep only the genes with at least
# excludedcounts counts >= excludedthresh
excludedcounts = parse(Int, ARGS[4])
excludedthresh = parse(Int, ARGS[5])

# For the initial conditions of the fit, use this threshold to separate noise
# below and expression above
initialthresh = excludedthresh

#######################################################
# Parameters for run
#######################################################
# function loadBCcounts(filepath, BC)
#     df = DataFrame(CSV.File(filepath, delim='\t', select=[Symbol(BC)]))
#     df[!, Symbol(BC)]
# end

function load_tsv(filepath)
    DataFrame(CSV.File(filepath, delim='\t'))
end

#######################################################
# Data Processing Loops
#######################################################

datapath = joinpath(datapathprefix, "$(cond)_"*join(replist))
outputpath = joinpath(datapath, "GeneFitData", "contaminantfits")
    
countpath = joinpath(datapath, "$(cond)_$(rep)_CellCounts.tsv")

noBCcells = parse(Int, readline(joinpath(datapath, "$(cond)_$(rep)_noBCcells.txt")))

runident = "$(cond)_$(rep)"

#####################
# Initial Count filter
# Load Data
println("Filtering Plasmids for $(runident)")

allcountsdf = load_tsv(countpath)
BClist = names(allcountsdf)[2:end]
println("\tFound $(length(BClist)) barcode types")

# Filter for sufficient counts to fit noise

overthresh = mapcols(col -> sum(c >= excludedthresh for c in col), allcountsdf[!, BClist])
keptBC = BClist[Vector(overthresh[1, :]) .> excludedcounts]


println("Found $(length(keptBC)) barcodes with sufficient counts for fitting the contaminant model")


parameterdf = DataFrame()
mkpath(outputpath)

# Loop over Barcodes that were kept
for BC in keptBC
    print("Beginning Fit for $(BC):")
    BCcountdict = countmap(allcountsdf[!, BC])
    BCcountdict[0] += noBCcells

    
    # fitting occurs here
    try
        if contdist == "Full"
            νi, γi, ρμi, αi, fi = MoMFull(BCcountdict, initialthresh)

            # totalcellcounts = sum(values(BCcountdict))
            # fmax = (totalcellcounts - BCcountdict[0] + BCcountdict[1])/totalcellcounts
            # [ν, γ, ρμ, α, f]
            lower = [0., 0., 0., 0., 0.]
            upper = [1., 1e3, 1e3, 10, 1.]
            initial = [νi, γi, ρμi, αi, fi]

        
            fitresults = fitgene_MLE_Full(BCcountdict, lower, upper, initial)
            
            print("\tFit complete.")

            # Generate outputs from fit
            # ν, γ, ρμ, α, f = fitresults.minimizer
            
            # dist for threshold and plots
            # gex(z) = GNB(z, ρμ, α)
            # gexν(z) = GNB(z, ρμ*ν, α)
            # gcont(z) = Gpois(gexν(z), f*γ)
            # gseq(z) = gcont(z)*(f*gex(z) + (1-f))

        elseif contdist == "Poisson"
            λi, ρμi, αi, fi = MoMPoisson(BCcountdict, initialthresh)
            # [λ, ρμ, α, f]
            lower = [0., 0., 0., 0.]
            upper = [1e2, 1e3, 1e2, 1.]
            initial = [λi, ρμi, αi, fi]

            fitresults = fitgene_MLE_Poisson(BCcountdict, lower, upper, initial)
            
            print("\tFit complete.")

        end

        gex, gcont, gseq = genfuncs(fitresults.minimizer...)
        f = fitresults.minimizer[end]

        cellcountvec, cellcountfreq = formcountfreq(BCcountdict)
        pnoexlist = taylor_expand(z->gcont(z)*gex(z), order = cellcountvec[end]).coeffs
        pseqlist = taylor_expand(gseq, order = cellcountvec[end]).coeffs
        pnolistlong = taylor_expand(gcont, order = cellcountvec[end]).coeffs
        
        threshold = cellcountvec[findfirst(f.*pnoexlist[1:length(cellcountvec)] .> (1-f).*pnolistlong[1:length(cellcountvec)])]	
        #save the parameters in a dataframe with genes across columns and parameters in rows

        parameterdf[!, BC] = [fitresults.minimizer..., threshold]
            
        print("\tGenerating Plots.")
        # for plotting, save each barcode fit in different dataframes
        mixturefitdf = DataFrame()
        
        mixturefitdf[!, BC] = cellcountfreq
        mixturefitdf[!, "$(BC)_contaminant"] = (1-f).*pnolistlong
        mixturefitdf[!, "$(BC)_expression"] = f.*pnoexlist
        mixturefitdf[!, "$(BC)_mixture"] = pseqlist

        CSV.write(joinpath(outputpath, runident*"_MixtureFit_$(BC)_$(contdist).tsv"), mixturefitdf, delim='\t')
        println("\tOutputs saved.")
    catch DomainError
        println("\tError encountered in Fit. Skipping...")
    end
end

CSV.write(joinpath(outputpath, runident*"_Parameters_$(contdist).tsv"), parameterdf, delim='\t')
