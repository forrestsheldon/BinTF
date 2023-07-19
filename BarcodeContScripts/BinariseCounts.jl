using Pkg
Pkg.activate("Project.toml")

using CSV, DataFrames
using TaylorSeries

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
fitpath = joinpath(datapath, "finalfits")
countpath = joinpath(datapathprefix, "$(cond)_"*join(replist), "$(cond)_$(rep)_CellCounts.tsv")
outputpath = joinpath(datapathprefix, "$(cond)_"*join(replist))

# Load parameters set by previous steps
parameter_df = DataFrame(CSV.File(joinpath(fitpath, "$(cond)_$(rep)_Parameters.tsv"), delim='\t'))
maxcountdict = Dict{String, Int}()

# Load Data
# Filter cells
cellBCdf = DataFrame(CSV.File(countpath))

# Create empty DataFrames
numcells = size(cellBCdf, 1)
binarydf = DataFrame("cell_barcode" => cellBCdf[!, "cell_barcode"], [BC => zeros(Int8, numcells) for BC in names(parameter_df)]...)
falseposdf = DataFrame("cell_barcode" => cellBCdf[!, "cell_barcode"], [BC => zeros(numcells) for BC in names(parameter_df)]...)
falsenegdf = DataFrame("cell_barcode" => cellBCdf[!, "cell_barcode"], [BC => zeros(numcells) for BC in names(parameter_df)]...)

# Fill them in
for BC in names(parameter_df)
    ν, γ, ρμ, α, f, thresh = parameter_df[!, BC]
    thresh = Integer(thresh)

    gex(z) = GNB(z, ρμ, α)
    gexν(z) = GNB(z, ρμ*ν, α)
    gcont(z) = Gpois(gexν(z), f*γ)
    gseq(z) = gcont(z)*(f*gex(z) + (1-f))
    gcontex(z) = gcont(z)*gex(z)

    maxcount = maximum(cellBCdf[!, BC])

    maxcountdict[BC] = maxcount

    pnolist = taylor_expand(gcont, order=maxcount).coeffs
    pnoexlist = taylor_expand(gcontex, order=maxcount).coeffs

    if sum(pnolist) > 1 pnolist = pnolist ./ sum(pnolist) end
    if sum(pnoexlist) > 1 pnoexlist = pnoexlist ./ sum(pnoexlist) end

    # See notebook for calculated error rates
    for (idx, count) in enumerate(cellBCdf[!, BC])
        if count >= thresh
            binarydf[idx, BC] = 1
            falseposdf[idx, BC] = pnolist[count+1]
            falsenegdf[idx, BC] = 1-pnoexlist[count+1]
        else
            binarydf[idx, BC] = 0
            falseposdf[idx, BC] = 1-pnolist[count+1]
            falsenegdf[idx, BC] = pnoexlist[count+1]
        end
    end

    CSV.write(joinpath(outputpath, "$(cond)_$(rep)_Labeled.tsv"), binarydf, delim='\t')
    CSV.write(joinpath(outputpath, "$(cond)_$(rep)_FalseNeg.tsv"), falsenegdf, delim='\t')
    CSV.write(joinpath(outputpath, "$(cond)_$(rep)_FalsePos.tsv"), falseposdf, delim='\t')

end


#Also want average error rates
avgerrdf = DataFrame()
for BC in names(parameter_df)
    ν, γ, ρμ, α, f, thresh = parameter_df[!, BC]
    thresh = Integer(thresh)

    gex(z) = GNB(z, ρμ, α)
    gexν(z) = GNB(z, ρμ*ν, α)
    gcont(z) = Gpois(gexν(z), f*γ)
    gseq(z) = gcont(z)*(f*gex(z) + (1-f))
    gcontex(z) = gcont(z)*gex(z)

    pnolist = taylor_expand(gcont, order=maxcountdict[BC]).coeffs
    pnoexlist = taylor_expand(gcontex, order=maxcountdict[BC]).coeffs
    
    if sum(pnolist) > 1 pnolist = pnolist ./ sum(pnolist) end
    if sum(pnoexlist) > 1 pnoexlist = pnoexlist ./ sum(pnoexlist) end
    
    #sum 0:thresh-1
    ϵT = sum(pnoexlist[1:thresh])
    # sum(thresh:Infty) = 1-sum(0:thresh-1)
    ϵnT = 1-sum(pnolist[1:thresh])

    avgerrdf[!, BC] = [ϵT, ϵnT]
end

CSV.write(joinpath(outputpath, "$(cond)_$(rep)_AvgErr.tsv"), avgerrdf, delim='\t')
