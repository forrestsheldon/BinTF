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

# to fit the global noise parameters keep only the genes with at least
# excludedcounts counts >= excludedthresh
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
ν_const, γ_const = [parse(Float64, s) for s in readlines(joinpath(datapath, "noisefits", runident*"_NoiseParameters.txt"))]
νγ_const = ν_const*γ_const

########
# Alternative Noise parameters
regparam = [parse(Float64, s) for s in readlines(joinpath(datapath, "noisefits", runident*"_NoiseParameters2.txt"))]

# Load Data
println("Loading data for $(cond) $(rep)")
cellhistogramdf = DataFrame(CSV.File(cellhistpath))

BClist = names(cellhistogramdf)
# keptBC = []
# for BC in BClist
#     cellcounts = cellhistogramdf[!, BC]
#     if sum(cellcounts[excludedthresh+1:end]) > excludedcounts
#         push!(keptBC, BC)
#     else
#         println("$(BC) excluded for insufficient counts")
#     end
# end
# println("Found $(length(keptBC)) barcodes with $(excludedcounts) counts over $(excludedthresh-1) for fitting")


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

    # # [ν, μ, r, f]
    # lower = [1e-4, 0., 0., 0.]
    # upper = [1., 1e3, 100, fmax]
    # initial = [νγ_const/γi, μi, ri, fi]

     # [ν, γ, μ, r, f]
     lower = [0., 0., 0., 0., 0.]
     upper = [1., 1e3, 1e3, 10, fmax]
     initial = [ν_const, γ_const, min(μi, 100.), min(ri, 5.), fi]
    
    # fitting occurs here
    try
        # fitresults = fitgene_MLE_νγ(celldict, νγ_const, lower, upper, initial)
        fitresults = fitgene_MLE_reg(celldict, regparam, lower, upper, initial)

        print("Fit complete. Generating outputs. ")
        # Generate outputs from fit
        # ν, μ, r, f = fitresults.minimizer
        # γ = νγ_const/ν
        ν, γ, μ, r, f = fitresults.minimizer
        gnof(z) = Gnoise(z, ν, f*γ, μ, r)
                
        # mixture fits
        cellcountvec, cellcountfreq = formcountfreq(celldict)
        gexf(z) = GNB(z, μ, r)
        gnoexf(z) = gexf(z)*gnof(z)
        gseqf(z) = f*gnoexf(z) + (1-f)*gnof(z)
        pnoexflist = taylor_expand(gnoexf, order = cellcountvec[end]).coeffs
        pseqflist = taylor_expand(gseqf, order = cellcountvec[end]).coeffs
        pnoflistlong = taylor_expand(gnof, order = cellcountvec[end]).coeffs
        
        threshold = cellcountvec[findfirst(f.*pnoexflist[1:length(cellcountvec)] .> (1-f).*pnoflistlong[1:length(cellcountvec)])]

        #save the parameters in a dataframe with genes across columns and parameters in rows
        parameterdf[!, BC] = [ν, γ, μ, r, f, threshold]
        
        # for plotting, save each barcode fit in different dataframes
        mixturefitdf = DataFrame()        
        mixturefitdf[!, BC] = cellcountfreq
        mixturefitdf[!, BC*"_noise"] = (1-f).*pnoflistlong
        mixturefitdf[!, BC*"_expression"] = f.*pnoexflist
        mixturefitdf[!, BC*"_mixture"] = pseqflist
        
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
