using Pkg
Pkg.activate("Project.toml")

using DataFrames, CSV
using TaylorSeries, Optim

include("../fitutils.jl")

#######################################################
# Set paths and parameters for screen
#######################################################
replist = readlines("replicates.txt")
condlist = readlines("condition.txt")
datapathprefix = readline("pathinfo.txt")

rep = ARGS[1]
cond = ARGS[2]

datapath = joinpath(datapathprefix, "$(cond)_"*join(replist), "GeneFitData")
outputpath = joinpath(datapath, "allfits")

cellhistpath = joinpath(datapath, "counthistograms", "$(cond)_$(rep)_cell_count_histograms.tsv")


#######################################################
# Parameters for run
#######################################################

# to check which threshold to use in MoM
excludedthresh = 6
excludedcounts = 30

# For the initial conditions of the fit, use this threshold to separate noise
# below and expression above
highthresh = 10
lowthresh = 2

#######################################################
# Data Processing Loops
#######################################################
runident = "$(cond)_$(rep)"
# Load the Noise Constant set by the last fit
contparam = [parse(Float64, s) for s in readlines(joinpath(datapath, "contaminantfits", runident*"_ContaminantParameters.txt"))]
γ_const, γ_MAD, νγ_const, νγ_MAD = contparam

# Load Data
println("Loading data for $(cond) $(rep)")
cellhistogramdf = DataFrame(CSV.File(cellhistpath))

BClist = names(cellhistogramdf)


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

#Refit genes with fixed noise and limits on the value of f given by a geometric increase in the number of counts

parameterdf = DataFrame()
    
mkpath(outputpath)
erroredBC = []
for BC in BClist

    print("Beginning Fit for $(BC): ")
    celldict = cellcountdicts[BC]

    if sum(nc for (c, nc) in celldict if c >= excludedthresh; init=0) > excludedcounts
        momthresh = highthresh
    else
        momthresh = lowthresh
    end

    νi, γi, ρμi, αi, fi = MoMinitial(celldict, momthresh)
    
    
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


     # [ν, γ, μ, α, f]
     lower = [0., 0., 0., 0., 0.]
     upper = [1., 1e3, 1e3, 10, fmax]
     initial = [νγ_const/γ_const, γ_const, min(ρμi, 100.), min(αi, 5.), fi]
    
    # fitting occurs here
    try
        fitresults = fitgene_MLE_reg(celldict, contparam, lower, upper, initial)

        print("Fit complete. Generating outputs. ")
        # Generate outputs from fit
 
        ν, γ, ρμ, α, f = fitresults.minimizer

        gex(z) = GNB(z, ρμ, α)
        gexν(z) = GNB(z, ρμ*ν, α)
        gcont(z) = Gpois(gexν(z), f*γ)
        gseq(z) = gcont(z)*(f*gex(z) + (1-f))
                
        # mixture fits
        cellcountvec, cellcountfreq = formcountfreq(celldict)
        pnoexlist = taylor_expand(z->gcont(z)*gex(z), order = cellcountvec[end]).coeffs
        pseqlist = taylor_expand(gseq, order = cellcountvec[end]).coeffs
        pnolistlong = taylor_expand(gcont, order = cellcountvec[end]).coeffs
        
        threshold = cellcountvec[findfirst(f.*pnoexlist[1:length(cellcountvec)] .> (1-f).*pnolistlong[1:length(cellcountvec)])]

        #save the parameters in a dataframe with genes across columns and parameters in rows
        parameterdf[!, BC] = [ν, γ, ρμ, α, f, threshold]
        
        # for plotting, save each barcode fit in different dataframes
        mixturefitdf = DataFrame()        
        mixturefitdf[!, BC] = cellcountfreq
        mixturefitdf[!, BC*"_noise"] = (1-f).*pnolistlong
        mixturefitdf[!, BC*"_expression"] = f.*pnoexlist
        mixturefitdf[!, BC*"_mixture"] = pseqlist
        
        CSV.write(joinpath(outputpath, runident*"_MixtureFit_$(BC).tsv"), mixturefitdf, delim='\t')
        println(" Outputs saved.")
    catch
        println("\tError encountered in fit. Skipping")
        push!(erroredBC, BC)
    end
end

CSV.write(joinpath(outputpath, runident*"_Parameters.tsv"), parameterdf, delim='\t')

if length(erroredBC) > 0
    println("Fitting Failed for $(erroredBC).")
    println("Check counts to see if fit is possible.")
end
