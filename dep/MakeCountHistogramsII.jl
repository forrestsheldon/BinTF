using Pkg
Pkg.activate("Project.toml")

using StatsBase
using CSV, DataFrames

#######################################################
# Set paths and parameters for screen
#######################################################
replist = readlines("replicates.txt")
condlist = readlines("condition.txt")
datapathprefix = readline("pathinfo.txt")



#######################################################
# Function Definitions
#######################################################
function formcountdict(BC, BCdf, zerobackground)

    genecounts = BCdf[:, BC]
    countdict = countmap(genecounts)
    countdict[0] += zerobackground
    
    return countdict
end

#######################################################
# Data Processing Loops
#######################################################
for cond in condlist
    # Data is original Count files
    datapath = joinpath(datapathprefix, "$(cond)_"*join(replist))
    outputpath = joinpath(datapath,"GeneFitData", "counthistograms")

    for rep in replist
        println("Beginning Binning for $(rep) $(cond)")

        # zerobcdroplets = parse(Int, readline(joinpath(datapath, "$(cond)_$(rep)_ZeroCounts.txt")))

        cellBCdf = DataFrame(CSV.File(joinpath(datapath, "$(cond)_$(rep)_CellCounts.tsv"), delim='\t'))
        # notcellBCdf = DataFrame(CSV.File(joinpath(datapath, "$(cond)_$(rep)_EmptyCounts.tsv"), delim='\t'))

        # Make histograms from countmaps and find the maximum count values needed
        BClist = names(cellBCdf)[2:end]

        cellBCdict = Dict()
        # emptydropletBCdict = Dict()
        
        for BC in BClist
            cellhist = formcountdict(BC, cellBCdf, 0)
            # droplethist = formcountdict(BC, notcellBCdf, zerobcdroplets)

            cellBCdict[BC] = cellhist
            # emptydropletBCdict[BC] = droplethist
        end
        print("Binning Complete, ")

        cellmaxcount = maximum([maximum(keys(BCdict)) for BCdict in values(cellBCdict)])
        # dropletmaxcount = maximum([maximum(keys(BCdict)) for BCdict in values(emptydropletBCdict)])

        # Place histograms into a dataframe as columns
        cellBChistdf = DataFrame()
        # emptydropletBChistdf = DataFrame()
        
        for BC in BClist
            cellBChistdf[!, BC] = [haskey(cellBCdict[BC], c) ? cellBCdict[BC][c] : 0 for c in 0:cellmaxcount]
            # emptydropletBChistdf[!, BC] = [haskey(emptydropletBCdict[BC], c) ? emptydropletBCdict[BC][c] : 0 for c in 0:dropletmaxcount]
        end

        mkpath(outputpath)
        CSV.write(joinpath(outputpath, "$(cond)_$(rep)_cell_count_histograms.tsv"), cellBChistdf, delim='\t')
        # CSV.write(joinpath(outputpath, "$(cond)_$(rep)_emptydroplet_count_histograms.tsv"), emptydropletBChistdf, delim='\t')
        println("Results Saved")
    end
end

