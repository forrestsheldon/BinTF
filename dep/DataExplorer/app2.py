import os
import streamlit as st
import pandas as pd
import plotly.express as px

DATA_DIR = "."  # Use this to refer to the current directory
replist = ["R11", "R12"]
markerlist = ["NG", "B"]

fileprefix = "bb20230124104525_NKv"
filesuffixes = ["cell_count_histograms.tsv", "emptydroplet_count_histograms.tsv"]

st.title("Screen Count Viewer and Plotter")

def get_file_path(rep, marker, file_suffix):
    return os.path.join(DATA_DIR, f"{fileprefix}_{rep}_{marker}_{file_suffix}")

def load_tsv_file(file_path):
    return pd.read_csv(file_path, sep='\t')

def plot_data(df_cell, df_droplet, selected_columns, normalize):
    df_plot = pd.concat([df_cell[selected_columns].reset_index(), df_droplet[selected_columns].reset_index()], keys=['cell', 'droplet'], names=['type']).reset_index()
    
    if normalize:
        df_plot[selected_columns] = df_plot[selected_columns].div(df_plot[selected_columns].sum(axis=0), axis=1)

    fig = px.line(df_plot, x="index", y=selected_columns, color='type', title="Line Plot of Selected Columns", hover_data=["type"])

    return fig

selected_rep = st.selectbox("Select a replicate:", replist)
selected_marker = st.selectbox("Select a marker:", markerlist)

file_path_cell = get_file_path(selected_rep, selected_marker, filesuffixes[0])
file_path_droplet = get_file_path(selected_rep, selected_marker, filesuffixes[1])

df_cell = load_tsv_file(file_path_cell)
df_droplet = load_tsv_file(file_path_droplet)

columns = df_cell.columns.tolist()

st.sidebar.subheader("Select columns to plot:")
selected_columns = [col for col in columns if st.sidebar.checkbox(col)]

if selected_columns:
    st.sidebar.subheader("Set axis limits:")
    normalize = st.sidebar.checkbox("Normalize by sum")

    # Prepare data for the plot
    fig = plot_data(df_cell, df_droplet, selected_columns, normalize)

    if normalize:
        step_size = 0.0001
        y_min, y_max = 0., 1.
    else:
        step_size = 1
        y_min, y_max = 0, max(df_cell.max().max(), df_droplet.max().max())

    x_min, x_max = st.sidebar.number_input("X-axis min value:", value=df_cell.index.min(), step=1), st.sidebar.number_input("X-axis max value:", value=df_cell.index.max(), step=1)
    y_min, y_max = st.sidebar.number_input("Y-axis min value:", value=y_min, step=step_size), st.sidebar.number_input("Y-axis max value:", value=y_max, step=step_size)

    # Update axis limits
    fig.update_xaxes(range=[x_min, x_max])
    fig.update_yaxes(range=[y_min, y_max])

    # Display the line plot
    st.write(fig)

else:
    st.warning("Please select at least one column to plot.")
