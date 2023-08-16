import streamlit as st
import subprocess
import os

def app(datadir, rep, cond, contdist):
    st.title("Data Requirements")
    st.text("""For these scripts to run, they need access to condition.txt which gives
the different conditions that the screen may have been run/sequenced under, e.g.
\twithoutAmp
\twithAmp
replicates.txt which is a list of replicates e.g.
\tR01
\tR02
and pathinfo.txt which contains the path to the data e.g.
\t/path/Data/

This will assume that the count data is in tab separated files like,
\t/path/Data/withoutAmp_R01R02/withoutAmp_R01_CellCounts.tsv
\t/path/Data/withoutAmp_R01R02/withoutAmp_R01_noBCcells.txt
with corresponding withAmp and R02 files.

The CellCounts files should contain a first column of cell-barcodes
and subsequent columns for each plasmid's barcode counts.
            
The noBCcells file should contain an integer with the number of cells captured
by GEX that did not contain any barcodes (and thus are usually excluded from the
CellCounts file).
            
Additionally, some plots will use an estimate of the number of unique plasmids in each cell.
For the data here, these are from ddPCR genomic reads that count the number of insertions
of any plasmid. These are floating point numbers in
\t/path/Data/withoutAmp_R01R02/R01_ddPCR.txt
with corresponding entries for withAMP and R02

            
The analysis runs for each replicate and condition separately, but there are opportunities to
check replicates against each other.
            
The Noise model can be selected as the Full model from the paper, or the Poisson
limiting distribution depending on your dataset.
""")