import os
import shutil
import streamlit as st
import glob
import re
import pandas as pd
import plotly.express as px
import json
import subprocess

def loadlist(fname):
    with open(fname) as f:
        replist = [line.rstrip() for line in f]
    return replist

def app(datadir, rep, cond):

    fitdir = os.path.join(datadir, "GeneFitData", "allfits")
    histdir = os.path.join(datadir, "GeneFitData", "counthistograms")
    outputdir = os.path.join(datadir, "GeneFitData", "finalfits")

    os.makedirs(outputdir, exist_ok=True)

    st.title("Fit Tuning")
    st.text("""Despite our best efforts, some genes will still display problematic fits.
This can result in a threshold of 0 (everything is included) or too large
of a threshold (almost everything is noise). This page allows you to adjust
the fraction of cells that contain the gene, f in the sidebar and tune the
fit of these problematic genes to your liking.

No genes should display a threshold of 0.
Reducing f will often help to increase the threshold. Increasing f will often
decrease the threshold for genes where noise is dominant.

For problematic fits, you can re-run the fit using a different initial condition or use a
method of moments run to obtain an approximate fit. Everything below the threshold is
assumed to be noise and everything equal to or above is assumed to be expression.""")
    
    cellhistpath = os.path.join(histdir, f"{cond}_{rep}_cell_count_histograms.tsv")
    df_cellhists = pd.read_csv(cellhistpath, sep='\t')
    gene_list = df_cellhists.columns.tolist()
    gene_list.sort()

    removedpath = os.path.join(outputdir, f"{cond}_{rep}_RemovedPlasmids.txt")
    try:
        removed_genes = loadlist(removedpath)
    except:
        removed_genes = []

    gene_list = [g for g in gene_list if g not in removed_genes]


    # # Pull gene list by first finding all the mixture fit files
    mix_file_pattern = f"{cond}_{rep}_MixtureFit_*.tsv"
    # mix_file_list = glob.glob(os.path.join(fitdir, mix_file_pattern))
    # # Regular expression pattern to extract the wildcard part
    mix_pattern = fr'{cond}_{rep}_MixtureFit_(.+)\.tsv'
    # gene_list = [re.search(mix_pattern, os.path.basename(file)).group(1) for file in mix_file_list]
    # gene_list.sort()

    # Pull finalised gene list by first finding all the mixture fit files
    final_file_list = glob.glob(os.path.join(outputdir, mix_file_pattern))
    # Regular expression pattern to extract the wildcard part
    finalised_gene_list = [re.search(mix_pattern, os.path.basename(file)).group(1) for file in final_file_list]
    unfinalised_genes = [g for g in gene_list if g not in finalised_gene_list]
    unfinalised_genes.sort()


    if unfinalised_genes:
        g_index = gene_list.index(unfinalised_genes[0])
    else:
        g_index = 0

    selected_gene = st.selectbox("Select a gene to plot:", gene_list, index=g_index)


    def load_tsv_file(file_path):
        return pd.read_csv(file_path, sep='\t')
            
    #Check if Fit has been Tuned previously by seeing if the Fit already exists in final
    mixture_path_final = os.path.join(outputdir, f"{cond}_{rep}_MixtureFit_{selected_gene}.tsv")
    mixture_path_init = os.path.join(fitdir, f"{cond}_{rep}_MixtureFit_{selected_gene}.tsv")
    param_path_final = os.path.join(outputdir, f"{cond}_{rep}_Parameters.tsv")
    param_path_init = os.path.join(fitdir, f"{cond}_{rep}_Parameters.tsv")

    if os.path.exists(mixture_path_final):
        st.text("Tuned Parameters Exist: Loading")

        param_df = load_tsv_file(param_path_final)
        mixture_df = load_tsv_file(mixture_path_final)
        yval = 0.01
    elif os.path.exists(mixture_path_init):
        st.text("No Tuned Parameter File Found: Loading Fit")
        param_df = load_tsv_file(param_path_init)
        mixture_df = load_tsv_file(mixture_path_init)
        yval = 0.01
    else:
        st.text("No Fit Found: Plotting Raw Counts")
        param_df = pd.DataFrame({selected_gene: [0., 0., 0., 0., 0., 0]})
        mixture_df = df_cellhists[[selected_gene]]
        yval = 10.
    
        

    ################
    # Sidebar Setup
    st.sidebar.title("Axis Limits")
    x_max = st.sidebar.number_input("X-axis max", value=100)
    y_max = st.sidebar.number_input("Y-axis max", value=yval, step=0.001, format="%.3f")

    
    geneparam = param_df[selected_gene].values

    st.sidebar.text(f"Volume Ratio ν:\t\t {round(geneparam[0], 4)}")
    st.sidebar.text(f"Bursting Rate γ:\t {round(geneparam[1], 1)}")
    st.sidebar.text(f"Mean Expression μ:\t {round(geneparam[2], 1)}")
    st.sidebar.text(f"Aggregation r:\t\t {round(geneparam[3], 4)}")
    st.sidebar.text(f"Threshold t:\t\t {int(geneparam[5])}")

    new_f = st.sidebar.number_input(f"Fraction f:\t", value=geneparam[4], step = 0.001, format="%.3f")
    new_thresh = st.sidebar.number_input(f"MoM thresh:\t", value=5, step = 1, format="%d")

    ##############
    # Plotting

    # Melt DataFrame to make it suitable for plotting using Plotly Express
    melted_mix_df = mixture_df.reset_index().melt(id_vars='index', var_name='Column', value_name='Value')

    # Create line plot using Plotly Express
    mix_fig = px.line(melted_mix_df, x='index', y='Value', color='Column')
    # Set x-axis and y-axis limits
    mix_fig.update_xaxes(range=[0, x_max])
    mix_fig.update_yaxes(range=[0, y_max])

    st.title("Mixture Fit")
    st.write(mix_fig)


    #################
    # Tune Fit

    if st.sidebar.button("Accept Fit"):
        # just copy parameters and fits
        if os.path.exists(param_path_final):
            param_df = load_tsv_file(param_path_final)
            param_df[selected_gene] = geneparam
        else:
            param_df = pd.DataFrame({selected_gene: geneparam})
        param_df.to_csv(param_path_final, sep='\t', index=False)

        if not os.path.exists(mixture_path_final):
            shutil.copy(mixture_path_init, mixture_path_final)
        st.experimental_rerun() 
            
    if st.sidebar.button("Rerun Fit"):
        # Rerun fit with new initial condition
        script_path = "./GeneBinScripts/FitSingleGene.jl"  # Replace with the actual path to your Julia script

        # Run the subprocess with real-time output
        with subprocess.Popen(["julia", script_path, rep, cond, str(new_thresh), selected_gene], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
            for line in proc.stdout:
                st.write(line.strip())

        # Wait for the subprocess to finish and get the return code
        return_code = proc.wait()

        if return_code == 0:
            st.success("Fit Completed")
            st.experimental_rerun() 
        else:
            st.error("The script encountered an error.")
            st.write("Error message:")
            st.write(proc.stderr)
           
    if st.sidebar.button("Rerun with fixed f"):
        # Run Tuning Script
        script_path = "./GeneBinScripts/FitSingleGene_fixedf.jl"  # Replace with the actual path to your Julia script

        # Run the subprocess with real-time output
        with subprocess.Popen(["julia", script_path, rep, cond, selected_gene, str(new_f), str(geneparam[0]), str(geneparam[1]), str(geneparam[2]), str(geneparam[3])], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
            for line in proc.stdout:
                st.write(line.strip())

        # Wait for the subprocess to finish and get the return code
        return_code = proc.wait()

        if return_code == 0:
            st.success("Final Fit Saved")
            st.experimental_rerun() 
        else:
            st.error("The script encountered an error.")
            st.write("Error message:")
            st.write(proc.stderr)
        
    if st.sidebar.button("Use MoM Fit"):
        # Run MoM Script
        script_path = "./GeneBinScripts/MoM_run.jl"  # Replace with the actual path to your Julia script

        # Run the subprocess with real-time output
        with subprocess.Popen(["julia", script_path, rep, cond, str(new_thresh), selected_gene], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
            for line in proc.stdout:
                st.write(line.strip())

        # Wait for the subprocess to finish and get the return code
        return_code = proc.wait()

        if return_code == 0:
            st.success("Final Fit Saved")
            st.experimental_rerun() 
        else:
            st.error("The script encountered an error.")
            st.write("Error message:")
            st.write(proc.stderr)
    
    if st.sidebar.button("Reset Fit"):
        try:
            os.remove(mixture_path_final)
        except OSError:
            pass
        st.experimental_rerun()

    if st.sidebar.button("Remove Gene"):
        removed_genes.append(selected_gene)
        f = open(removedpath, 'w')
        f.writelines(removed_genes)
        f.close()
        st.experimental_rerun()

    ##########################
    # Display Unfinalised Genes
    st.sidebar.markdown("Unchecked Genes:")
    for gene in unfinalised_genes:
        st.sidebar.markdown(f"- {gene}")


    if st.sidebar.button("Binarise Counts"):
        script_path = "./GeneBinScripts/BinariseCounts.jl"  # Replace with the actual path to your Julia script

        # Run the subprocess with real-time output
        with subprocess.Popen(["julia", script_path, rep, cond], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
            for line in proc.stdout:
                st.write(line.strip())

        # Wait for the subprocess to finish and get the return code
        return_code = proc.wait()

        if return_code == 0:
            st.success("Binarisation Complete!")


        else:
            st.error("The script encountered an error.")
            st.write("Error message:")
            st.write(proc.stderr)