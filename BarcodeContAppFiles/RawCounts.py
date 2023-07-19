import os
import streamlit as st
import pandas as pd
import plotly.express as px
import subprocess

def loadlist(fname):
    with open(fname) as f:
        list = [line.rstrip() for line in f]
    return list

def app(datadir, rep, cond):
    #####################
    # Set Path Variables
    #####################
    # Data is histogram files
    histdir = os.path.join(datadir, "GeneFitData", "counthistograms")


    #####################
    # Functions
    #####################
    def get_file_path(rep, cond):
        return os.path.join(histdir, f"{cond}_{rep}_cell_count_histograms.tsv")

    def load_tsv_file(file_path):
        return pd.read_csv(file_path, sep='\t')

    def plot_data(df_cell, selected_columns, normalize):
        df_plot = df_cell[selected_columns].reset_index()
        
        if normalize:
            df_plot[selected_columns] = df_plot[selected_columns].div(df_plot[selected_columns].sum(axis=0), axis=1)

        fig = px.line(df_plot, x="index", y=selected_columns, title="Line Plot of Selected Columns")

        return fig

    
    #####################
    # Pick which run to Plot
    #####################
    st.title("Raw Count Explorer")
    st.text("""Check counts for anomalous genes and make note of special cases.
    Zoom axes using controls at the bottom of the sidebar.
    Normalising counts with the checkbox, most genes will be viewable with y-limit near 0.01.""")
 

    file_path_cell = get_file_path(rep, cond)

    df_cell = load_tsv_file(file_path_cell)

    columns = df_cell.columns.tolist()

    #####################
    # Pick which Gene to Plot
    #####################
    st.sidebar.subheader("Select columns to plot:")
    selected_columns = [col for col in columns if st.sidebar.checkbox(col)]

    if selected_columns:
        st.sidebar.subheader("Set axis limits:")
        normalize = st.sidebar.checkbox("Normalize by sum")

        # Prepare data for the plot
        fig = plot_data(df_cell, selected_columns, normalize)

        if normalize:
            step_size = 0.0001
            y_min, y_max = 0., 1.
            fmtstring = "%.3f"
        else:
            step_size = 1
            y_min, y_max = 0, df_cell.max().max()
            fmtstring="%d"

        x_min, x_max = st.sidebar.number_input("X-axis min value:", value=df_cell.index.min(), step=1), st.sidebar.number_input("X-axis max value:", value=df_cell.index.max(), step=1)
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
    numcounts = st.number_input("Select plasmids with at least this many counts:", value=200, step = 1)
    initthresh = st.number_input("With at least this magnitude:", value = 10, step = 1)

    if st.button("Fit Contaminants: This will take a while"):
        script_path = "./BarcodeContScripts/FitContaminantParameters.jl"  # Replace with the actual path to your Julia script

        with subprocess.Popen(["julia", script_path, rep, cond, str(numcounts), str(initthresh)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
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



