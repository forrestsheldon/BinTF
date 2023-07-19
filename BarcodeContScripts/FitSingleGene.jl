using Pkg
Pkg.activate("Project.toml")

using DataFrames, CSV
using TaylorSeries, Optim

include("../fitutils.jl")

#######################################################
# Set paths and parameters for fit
#######################################################
replist = readlines("replicates.txt")
condlist = readlines("condition.txt")
datapathprefix = readline("pathinfo.txt")

rep = ARGS[1]
cond = ARGS[2]
initialthresh = parse(Int64, ARGS[3])
BC = ARGS[4]


datapath = joinpath(datapathprefix, "$(cond)_"*join(replist), "GeneFitData")
outputpath = joinpath(datapath, "allfits")
mkpath(outputpath)

runident = "$(cond)_$(rep)"
ν_const, γ_const = [parse(Float64, s) for s in readlines(joinpath(datapath, "noisefits", runident*"_NoiseParameters.txt"))]
νγ_const = ν_const*γ_const

########
# Alternative Noise parameters
regparam = [parse(Float64, s) for s in readlines(joinpath(datapath, "noisefits", runident*"_NoiseParameters2.txt"))]

#######################################################
# Load Counts
#######################################################
cellhistpath = joinpath(datapath, "counthistograms", "$(cond)_$(rep)_cell_count_histograms.tsv")

println("Loading data for $(cond) $(rep) $(BC)")
cellcounts = DataFrame(CSV.File(cellhistpath, select= [BC]))[!, BC]

celldict = Dict{Int, Int}()

for (idx, numcells) in enumerate(cellcounts)
    count = idx - 1
    if numcells != 0
        celldict[count] = numcells
    end
end

#######################################################
# Run Fit
#######################################################

νi, γi, μi, ri, fi = MoMinitial(celldict, initialthresh)

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

# [ν, γ, μ, r, f]
lower = [0., 0., 0., 0., 0.]
upper = [1., 1e3, 1e3, 10, fmax]
initial = [ν_const, γ_const, min(μi, 100.), min(ri, 5.), fi]

# fitting occurs here
# fitresults = fitgene_MLE_νγ(celldict, νγ_const, lower, upper, initial)

#######################################################
# Generate Outputs
#######################################################
fitresults = fitgene_MLE_reg(celldict, regparam, lower, upper, initial)
ν, γ, μ, r, f = fitresults.minimizer

gex(z) = GNB(z, μ, r)
gno(z) = Gnoise(z, ν, f*γ, μ, r)
gnoex(z) = gno(z)*gex(z)
gseq(z) = f*gnoex(z) + (1-f)*gno(z)


cellcountvec, cellcountfreq = formcountfreq(celldict)
pnoexlist = taylor_expand(gnoex, order = cellcountvec[end]).coeffs
pseqlist = taylor_expand(gseq, order = cellcountvec[end]).coeffs
pnolistlong = taylor_expand(gno, order = cellcountvec[end]).coeffs

threshold = cellcountvec[findfirst(f.*pnoexlist[1:length(cellcountvec)] .> (1-f).*pnolistlong[1:length(cellcountvec)])]

#######################################################
# Save Results
#######################################################
# Check if parameter file exists and load if so
outputparamfilepath = joinpath(outputpath, "$(cond)_$(rep)_Parameters.tsv")

parameterdf = DataFrame(CSV.File(outputparamfilepath, delim='\t'))
#save the parameters in a dataframe with genes across columns and parameters in rows
parameterdf[!, BC] = [ν, γ, μ, r, f, threshold]

CSV.write(outputparamfilepath, parameterdf, delim='\t')

# for plotting, save each barcode fit in different dataframes
mixturefitdf = DataFrame()
mixturefitdf[!, BC] = cellcountfreq
mixturefitdf[!, BC*"_noise"] = (1-f).*pnolistlong
mixturefitdf[!, BC*"_expression"] = f.*pnoexlist
mixturefitdf[!, BC*"_mixture"] = pseqlist

CSV.write(joinpath(outputpath, "$(cond)_$(rep)_MixtureFit_$(BC).tsv"), mixturefitdf, delim='\t')
println(" Outputs saved.")




