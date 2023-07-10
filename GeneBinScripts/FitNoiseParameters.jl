using Pkg
Pkg.activate("Project.toml")

using DataFrames, CSV
using TaylorSeries, Optim

include("../fitutils.jl")

#######################################################
# Set parameters for screen
#######################################################
replist = readlines("replicates.txt")
condlist = readlines("condition.txt")
datapathprefix = readline("pathinfo.txt")

rep = ARGS[1]
cond = ARGS[2]

#######################################################
# Parameters for run
#######################################################
# to fit the global noise parameters keep only the genes with at least
# excludedcounts counts >= excludedthresh
excludedcounts = parse(Int, ARGS[3])
excludedthresh = parse(Int, ARGS[4])


# For the initial conditions of the fit, use this threshold to separate noise
# below and expression above
initialthresh = 10



#######################################################
# Data Processing Loops
#######################################################

# Data is histogram files
datapath = joinpath(datapathprefix, "$(cond)_"*join(replist), "GeneFitData")
outputpath = joinpath(datapath, "noisefits")
    
cellhistpath = joinpath(datapath, "counthistograms", "$(cond)_$(rep)_cell_count_histograms.tsv")

runident = "$(cond)_$(rep)"

# Load Data
println("Loading data for $(runident)")
cellhistogramdf = DataFrame(CSV.File(cellhistpath))
# Filter for sufficient counts to fit noise
BClist = names(cellhistogramdf)
keptBC = []
for BC in BClist
    cellcounts = cellhistogramdf[!, BC]
    if sum(cellcounts[excludedthresh+1:end]) > excludedcounts
        push!(keptBC, BC)
    end
end
println("Found $(length(keptBC)) barcodes with sufficient counts for fitting the noise model")

cellcountdicts = Dict{String, Dict{Int, Int}}()

# Form count dictionaries
for BC in keptBC
    cellcountdicts[BC] = Dict{Int, Int}()
    
    cellcounts = cellhistogramdf[!, BC]

    for (idx, numcells) in enumerate(cellcounts)
        count = idx-1
        if numcells != 0
            cellcountdicts[BC][count] = numcells
        end
    end
end


parameterdf = DataFrame()

mkpath(outputpath)

# Loop over Barcodes that were kept
for BC in keptBC
    print("Beginning Fit for $(BC):")
    celldict = cellcountdicts[BC]
    
    νi, γi, ρμi, αi, fi = MoMinitial(celldict, initialthresh)

    totalcellcounts = sum(values(celldict))
    # totalemptycounts = sum(values(emptydict))
    
    fmax = (totalcellcounts - celldict[0] + celldict[1])/totalcellcounts
    

    # [ν, γ, ρμ, α, f]
    lower = [0., 0., 0., 0., 0.]
    upper = [1., 1e3, 1e3, 10, fmax]
    initial = [νi, γi,ρμi, αi, fi]

    # fitting occurs here
    try
        fitresults = fitgene_MLE(celldict, lower, upper, initial)
        
        println("\tFit complete. Generating outputs")
        # Generate outputs from fit
        ν, γ, ρμ, α, f = fitresults.minimizer
        
        # mixture fitting
        cellcountvec, cellcountfreq = formcountfreq(celldict)
        gex(z) = GNB(z, ρμ, α)
        gexν(z) = GNB(z, ρμ*ν, α)
        gcont(z) = Gpois(gexν(z), f*γ)
        gseq(z) = gcont(z)*(f*gex(z) + (1-f))

        pnoexlist = taylor_expand(z->gcont(z)*gex(z), order = cellcountvec[end]).coeffs
        pseqlist = taylor_expand(gseq, order = cellcountvec[end]).coeffs
        pnolistlong = taylor_expand(gcont, order = cellcountvec[end]).coeffs
        
        threshold = cellcountvec[findfirst(f.*pnoexlist[1:length(cellcountvec)] .> (1-f).*pnolistlong[1:length(cellcountvec)])]	
        #save the parameters in a dataframe with genes across columns and parameters in rows

        parameterdf[!, BC] = [fitresults.minimizer..., threshold]
        
        # for plotting, save each barcode fit in different dataframes
        mixturefitdf = DataFrame()
        
        mixturefitdf[!, BC] = cellcountfreq
        mixturefitdf[!, BC*"_contaminant"] = (1-f).*pnolistlong
        mixturefitdf[!, BC*"_expression"] = f.*pnoexlist
        mixturefitdf[!, BC*"_mixture"] = pseqlist

        CSV.write(joinpath(outputpath, runident*"_MixtureFit_$(BC).tsv"), mixturefitdf, delim='\t')
        println("\tOutputs saved.")
    catch DomainError
        println("\tError encountered in Fit. Skipping...")
    end
end

CSV.write(joinpath(outputpath, runident*"_Parameters.tsv"), parameterdf, delim='\t')
