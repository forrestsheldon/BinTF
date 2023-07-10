import streamlit as st
import subprocess
import os

def app(datadir, rep, cond):
    st.title("Make Histograms")
    st.text("This will bin counts to allow plotting and fitting.")
    st.text("""For these scripts to run, they need access to condition.txt which gives
the conditions for the screen, e.g.
\tPreAmp
\tPostAmp
replicates.txt which is a list of replicates e.g.
\tR01
\tR02
and pathinfo.txt which contains the path to the data e.g.
\t/path/Data/
which will assume that the count data is in tab separated files like,
\t/path/Data/PreAmp_R01R02/PreAmp_R01_CellCounts.tsv
\t/path/Data/PreAmp_R01R02/PreAmp_R01_ZeroCounts.txt
with corresponding R02 files.

Note that by default this will run over all replicates and conditions given.""")


    if st.button("Histogram Counts"):
        script_path = "./GeneBinScripts/MakeCountHistograms.jl"  # Replace with the actual path to your Julia script

        # Run the subprocess with real-time output
        with subprocess.Popen(["julia", script_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
            for line in proc.stdout:
                st.write(line.strip())

        # Wait for the subprocess to finish and get the return code
        return_code = proc.wait()

        if return_code == 0:
            st.success("Binning Complete!")
        else:
            st.error("The script encountered an error.")
            st.write("Error message:")
            st.write(proc.stderr)