using Pkg
Pkg.activate("Project.toml")

using DataFrames, CSV
using TaylorSeries, Optim

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

runident = "$(cond)_$(rep)"


datapath = joinpath(datapathprefix, "$(cond)_"*join(replist))
outputpath = joinpath(datapath, "GeneFitData", "allfits")

countpath = joinpath(datapath, "$(cond)_$(rep)_CellCounts.tsv")
parampath = joinpath(datapath, "GeneFitData", "contaminantfits", "$(cond)_$(rep)_FixedParameters_$(contdist).txt")

noBCcells = parse(Int, readline(joinpath(datapath, "$(cond)_$(rep)_noBCcells.txt")))

#######################################################
# Parameters for run
#######################################################

function loadBCcounts(filepath, BC)
    df = DataFrame(CSV.File(filepath, delim='\t', select=[Symbol(BC)]))
    df[!, Symbol(BC)]
end

#######################################################
# Data Processing Loops
#######################################################
# Load the Noise Constant set by the last fit
contparam = [parse(Float64, s) for s in readlines(parampath)]

# Load Data
println("Loading data for $(cond) $(rep)")
allcountsdf = load_tsv(countpath)
BClist = names(allcountsdf)[2:end]


println("\tFound $(length(BClist)) barcode types. Attempting fits.")


parameterdf = DataFrame()
    
mkpath(outputpath)
erroredBC = []
for BC in BClist

    print("Beginning Fit for $(BC): ")

    BCcountdict = countmap(allcountsdf[!, BC])
    BCcountdict[0] += noBCcells

    try
        if contdist == "Full"
    
           ν, γ, α = contparam

            ρμ, f = MoMFullfμ(BCcountdict, ν, γ, α)
            
            print("\tFit complete. ")
            gex, gcont, gseq = genfuncs(ν, γ, ρμ, α, f)


        elseif contdist == "Poisson"
            λ, α = contparam

            ρμ, f = MoMPoissonfμ(BCcountdict, λ, α)


            print("Fit complete. ")
            gex, gcont, gseq = genfuncs(λ, ρμ, α, f)

        end

        cellcountvec, cellcountfreq = formcountfreq(BCcountdict)
        pnoexlist = taylor_expand(z->gcont(z)*gex(z), order = cellcountvec[end]).coeffs
        pseqlist = taylor_expand(gseq, order = cellcountvec[end]).coeffs
        pnolistlong = taylor_expand(gcont, order = cellcountvec[end]).coeffs
        
        threshold = cellcountvec[findfirst(f.*pnoexlist[1:length(cellcountvec)] .> (1-f).*pnolistlong[1:length(cellcountvec)])]	
        #save the parameters in a dataframe with genes across columns and parameters in rows
        if contdist == "Full"
            parameterdf[!, BC] = [ν, γ, ρμ, α, f, threshold]
        elseif contdist == "Poisson"
            parameterdf[!, BC] = [λ, ρμ, α, f, threshold]
        end
            
        print("\tGenerating Plots.")
        # for plotting, save each barcode fit in different dataframes
        mixturefitdf = DataFrame()
        
        mixturefitdf[!, BC] = cellcountfreq
        mixturefitdf[!, "$(BC)_contaminant"] = (1-f).*pnolistlong
        mixturefitdf[!, "$(BC)_expression"] = f.*pnoexlist
        mixturefitdf[!, "$(BC)_mixture"] = pseqlist

        CSV.write(joinpath(outputpath, runident*"_MixtureFit_$(BC)_$(contdist).tsv"), mixturefitdf, delim='\t')
        println("\tOutputs saved.")
    catch
        println("\tError encountered in fit. Skipping")
        push!(erroredBC, BC)
    end
end

CSV.write(joinpath(outputpath, runident*"_Parameters_$(contdist).tsv"), parameterdf, delim='\t')

if length(erroredBC) > 0
    println("Fitting Failed for $(erroredBC).")
    println("Check counts to see if fit is possible.")
end
