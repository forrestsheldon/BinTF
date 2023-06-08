import os
import streamlit as st
import pandas as pd
import plotly.express as px

DATA_DIR = "."  # Use this to refer to the current directory

def get_tsv_files(data_dir):
    return [f for f in os.listdir(data_dir) if f.endswith('.tsv')]

def load_tsv_file(file_path):
    return pd.read_csv(file_path, sep='\t')

st.title("TSV File Viewer and Plotter")

# Get the list of TSV files and create a select box to choose a file
tsv_files = get_tsv_files(DATA_DIR)
selected_file = st.selectbox("Select a TSV file:", tsv_files)

if selected_file:
    file_path = os.path.join(DATA_DIR, selected_file)
    df = load_tsv_file(file_path)

    # Checkboxes for column selection in the sidebar
    columns = df.columns.tolist()
    st.sidebar.subheader("Select columns to plot:")
    selected_columns = [col for col in columns if st.sidebar.checkbox(col)]

    if selected_columns:
        # Create a line plot using the selected columns
        fig = px.line(df.reset_index(), x="index", y=selected_columns, title="Line Plot of Selected Columns")

        # Input fields for setting axis limits in the sidebar
        st.sidebar.subheader("Set axis limits:")
        x_min, x_max = st.sidebar.number_input("X-axis min value:", value=df.index.min(), step=1), st.sidebar.number_input("X-axis max value:", value=df.index.max(), step=1)
        y_min, y_max = st.sidebar.number_input("Y-axis min value:", value=df[selected_columns].min().min(), step=1), st.sidebar.number_input("Y-axis max value:", value=df[selected_columns].max().max(), step=1)

        # Update axis limits
        fig.update_xaxes(range=[x_min, x_max])
        fig.update_yaxes(range=[y_min, y_max])

        # Display the line plot
        st.write(fig)
    else:
        st.warning("Please select at least one column to plot.")
else:
    st.warning("No TSV files found in the specified directory.")
