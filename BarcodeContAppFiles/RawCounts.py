import os
import streamlit as st
import pandas as pd
import plotly.graph_objs as go
import subprocess

def loadlist(fname):
    with open(fname) as f:
        list = [line.rstrip() for line in f]
    return list

def get_plasmid_names(file_path):
    df = pd.read_csv(file_path, sep='\t', nrows=0)  # Load only the first (header) row
    return df.columns.tolist()[1:]

def load_plasmid_counts(filepath, selected_plasmids):
    return pd.read_csv(filepath, sep='\t', usecols=selected_plasmids)

def plot_data(rawcount_df, noBCcells, normalize):

    fig = go.Figure()
    maxcells = 0
    maxcount = 0

    for plasmid_col in rawcount_df:
        counthist = rawcount_df[plasmid_col].value_counts().sort_index()
        counthist.loc[0] += noBCcells
        
        if normalize:
            counthist = counthist / counthist.sum()

        fig.add_trace(go.Scatter(x=counthist.index, y=counthist.values, mode="lines", name=plasmid_col))
        
        maxcells = max(maxcells, max(counthist.values))
        maxcount = max(maxcount, max(counthist.index))

    return fig, maxcells, maxcount

def app(datadir, rep, cond, contdist):
    #####################
    # Set Path Variables
    #####################
    # Data is histogram files
    # countdir = os.path.join(datadir, "GeneFitData", "counthistograms")


    #####################
    # Functions
    #####################
    def get_CellCount_path(rep, cond):
        return os.path.join(datadir, f"{cond}_{rep}_CellCounts.tsv")
    
    def get_noBCcells(rep, cond):
        with open(os.path.join(datadir, f"{cond}_{rep}_noBCcells.txt"), 'r') as file:
            noBCcells = int(file.read().strip())
        return noBCcells
        


    
    #####################
    # Pick which run to Plot
    #####################
    st.title("Raw Count Explorer")
    st.text("""Check counts for anomalous genes and make note of special cases.
    Zoom axes using controls at the bottom of the sidebar.
    Normalising counts with the checkbox, most genes will be viewable with y-limit near 0.01.""")

    file_path_cell = get_CellCount_path(rep, cond)
    noBCcells = get_noBCcells(rep, cond)

    plasmid_names = get_plasmid_names(file_path_cell)
    plasmid_names.sort()

    #####################
    # Pick which Gene to Plot
    #####################
    st.sidebar.subheader("Select plasmids to plot:")
    selected_plasmids = [col for col in plasmid_names if st.sidebar.checkbox(col)]

    st.sidebar.subheader("Set axis limits:")
    normalize = st.sidebar.checkbox("Normalize by sum")


    if selected_plasmids:

        rawcount_df = load_plasmid_counts(file_path_cell, selected_plasmids)

        # Prepare data for the plot
        fig, maxcells, maxcount = plot_data(rawcount_df, noBCcells, normalize)

        if normalize:
            step_size = 0.0001
            y_min, y_max = 0., 1.
            fmtstring = "%.3f"
        else:
            step_size = 1
            y_min, y_max = 0, maxcells
            fmtstring="%d"

        x_min, x_max = st.sidebar.number_input("X-axis min value:", value=0, step=1), st.sidebar.number_input("X-axis max value:", value=maxcount, step=1)
        y_min, y_max = st.sidebar.number_input("Y-axis min value:", value=y_min, step=step_size), st.sidebar.number_input("Y-axis max value:", value=y_max, step=step_size, format=fmtstring)

        # Update axis limits
        fig.update_xaxes(range=[x_min, x_max])
        fig.update_yaxes(range=[y_min, y_max])

        # Display the line plot
        st.write(fig)

    else:
        st.warning("Please select at least one column to plot.")


    #####################
    # Run Fits
    #####################

    st.title("Run Contaminant Fits")

    st.text("""This will run the fitting algorithm on high count
genes to estimate the contaminant parameters.""")
    numcounts = st.number_input("Select genes that appear in at least this many cells:", value=200, step = 1)
    initthresh = st.number_input("With a count at least this high:", value = 10, step = 1)

    if st.button("Fit Contaminants: This will take a while"):
        script_path = "./BarcodeContScripts/FitContaminantParameters.jl"  # Replace with the actual path to your Julia script

        with subprocess.Popen(["julia", script_path, rep, cond, contdist, str(numcounts), str(initthresh)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
            for line in proc.stdout:
                st.write(line.strip())

        # Wait for the subprocess to finish and get the return code
        return_code = proc.wait()

        if return_code == 0:
            st.success("Fitting Complete!")
        else:
            st.error("The script encountered an error.")
            st.write("Error message:")
            st.write(proc.stderr)



