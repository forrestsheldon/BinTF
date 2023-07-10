import os
import streamlit as st
import glob
import re
import pandas as pd
import plotly.express as px
import json

def app(datadir, rep, cond):
    #####################
    # Set Path Variables
    #####################
    fitdir = os.path.join(datadir, "GeneFitData", "noisefits")

    #####################
    # Functions
    #####################

    def get_mixture_path(rep, cond, gene):
        return os.path.join(fitdir, f"{cond}_{rep}_MixtureFit_{gene}.tsv")
    
    def get_param_path(rep, cond):
        return os.path.join(fitdir, f"{cond}_{rep}_Parameters.tsv")

    def load_tsv_file(file_path):
        return pd.read_csv(file_path, sep='\t')

    #####################
    # Pick which run to Plot
    #####################
    st.title("Noise Fit Explorer")
    st.write("Examine Noise Fits for high count genes and remove any that seem to be problematic.")

    paramdf = load_tsv_file(get_param_path(rep, cond))
    st.dataframe(paramdf)

    #####################
    # Select Genes for Fit
    #####################
    # # Pull gene list by first finding all the mixture fit files
    # mixture_file_pattern = f"{cond}_{rep}_MixtureFit_*.tsv"
    # mixture_file_list = glob.glob(os.path.join(fitdir, mixture_file_pattern))

    # # Regular expression pattern to extract the wildcard part
    # mixture_pattern = fr'{cond}_{rep}_MixtureFit_(.+)\.tsv'

    # gene_list = [re.search(mixture_pattern, os.path.basename(file)).group(1) for file in mixture_file_list]
    # gene_list.sort()

    gene_list = paramdf.columns.tolist()
    gene_list.sort()

    # Create a list of checkboxes in the sidebar
    st.sidebar.title("Select Good Fits:")
    selected_genes = {}

    for gene in gene_list:
        # Use the wildcard part as the label for the checkbox
        selected_genes[gene] = st.sidebar.checkbox(gene, value=True)

    # Filter the selected files
    selected_genes_list = [gene for gene, selected in selected_genes.items() if selected]

    # Add a button to write the selected genes to a file
    if st.sidebar.button("Save Selected Genes"):
        # Write the selected genes to a file as a list of strings
        with open(os.path.join(fitdir, f"{cond}_{rep}_SelectedFits.json"), "w") as f:
            json.dump(selected_genes_list, f)

        st.success("Selected genes saved to selected_genes.json")

    st.sidebar.title("Axis Limits")
    x_max = st.sidebar.number_input("X-axis max", value=100)
    y_max = st.sidebar.number_input("Y-axis max", value=0.01, step=0.001, format="%.3f")


    #####################
    # Plotting
    #####################

    plotting_gene = st.selectbox("Select a gene to plot:", gene_list)

    mix_path = get_mixture_path(rep, cond, plotting_gene)
    mix_df = load_tsv_file(mix_path)

    # # Shift indices by 1
    # shifted_indices = mix_df.index - 1

    # Melt DataFrame to make it suitable for plotting using Plotly Express
    melted_mix_df = mix_df.reset_index().melt(id_vars='index', var_name='Column', value_name='Value')

    # Create line plot using Plotly Express
    mix_fig = px.line(melted_mix_df, x='index', y='Value', color='Column')
    # Set x-axis and y-axis limits
    mix_fig.update_xaxes(range=[0, x_max], title_text='Count')
    mix_fig.update_yaxes(range=[0, y_max], title_text='Fraction of Cells')

    st.title("Mixture Fit")
    st.write(mix_fig)




