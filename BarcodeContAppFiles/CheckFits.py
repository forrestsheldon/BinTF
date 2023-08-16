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
        list = [line.rstrip() for line in f]
    return list

def get_plasmid_names(file_path):
    df = pd.read_csv(file_path, sep='\t', nrows=0)  # Load only the first (header) row
    return df.columns.tolist()[1:]

def load_plasmid_counts(filepath, selected_plasmids):
    return pd.read_csv(filepath, sep='\t', usecols=selected_plasmids)

def app(datadir, rep, cond, contdist):

    fitdir = os.path.join(datadir, "GeneFitData", "allfits")
    outputdir = os.path.join(datadir, "GeneFitData", "finalfits")

    os.makedirs(outputdir, exist_ok=True)

    def get_noBCcells(rep, cond):
        with open(os.path.join(datadir, f"{cond}_{rep}_noBCcells.txt"), 'r') as file:
            noBCcells = int(file.read().strip())
        return noBCcells

    st.title("Fit Tuning")
    st.text("""Despite our best efforts, some genes will still display problematic fits.""")
    
    file_path_cell = os.path.join(datadir, f"{cond}_{rep}_CellCounts.tsv")
    plasmid_names = get_plasmid_names(file_path_cell)
    plasmid_names.sort()
    noBCcells = get_noBCcells(rep, cond)


    removedpath = os.path.join(outputdir, f"{cond}_{rep}_RemovedPlasmids_{contdist}.txt")
    try:
        removed_plasmids = loadlist(removedpath)
    except:
        removed_plasmids = []

    plasmid_list = [g for g in plasmid_names if g not in removed_plasmids]


    # # Pull gene list by first finding all the mixture fit files
    mix_file_pattern = f"{cond}_{rep}_MixtureFit_*_{contdist}.tsv"
    mix_pattern = fr'{cond}_{rep}_MixtureFit_(.+)_{contdist}.tsv'

    # Pull finalised gene list by first finding all the mixture fit files
    final_file_list = glob.glob(os.path.join(outputdir, mix_file_pattern))
    # Regular expression pattern to extract the wildcard part
    finalised_gene_list = [re.search(mix_pattern, os.path.basename(file)).group(1) for file in final_file_list]

    unfinalised_genes = [g for g in plasmid_list if g not in finalised_gene_list]
    unfinalised_genes.sort()


    if unfinalised_genes:
        g_index = plasmid_list.index(unfinalised_genes[0])
    else:
        g_index = 0

    selected_gene = st.selectbox("Select a gene to plot:", plasmid_list, index=g_index)


    def load_tsv_file(file_path):
        return pd.read_csv(file_path, sep='\t')
            
    #Check if Fit has been Tuned previously by seeing if the Fit already exists in final
    mixture_path_final = os.path.join(outputdir, f"{cond}_{rep}_MixtureFit_{selected_gene}_{contdist}.tsv")
    mixture_path_init = os.path.join(fitdir, f"{cond}_{rep}_MixtureFit_{selected_gene}_{contdist}.tsv")

    param_path_final = os.path.join(outputdir, f"{cond}_{rep}_Parameters_{contdist}.tsv")
    param_path_init = os.path.join(fitdir, f"{cond}_{rep}_Parameters_{contdist}.tsv")

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
        rawcount_df = load_plasmid_counts(file_path_cell, [selected_gene])
        param_df = pd.DataFrame({selected_gene: [0., 0., 0., 0., 0., 0]})

        mixture_df = rawcount_df[selected_gene].value_counts().sort_index()
        mixture_df.loc[0] += noBCcells
        yval = 10.
    
        

    ################
    # Sidebar Setup
    st.sidebar.title("Axis Limits")
    x_max = st.sidebar.number_input("X-axis max", value=100)
    y_max = st.sidebar.number_input("Y-axis max", value=yval, step=0.001, format="%.3f")

    
    geneparam = param_df[selected_gene].values

    if contdist == "Full":
        st.sidebar.text(f"Volume Ratio ν:\t\t {round(geneparam[0], 4)}")
        st.sidebar.text(f"Bursting Rate γ:\t {round(geneparam[1], 1)}")
        st.sidebar.text(f"Mean Expression ρμ:\t {round(geneparam[2], 1)}")
        st.sidebar.text(f"Dispersion α:\t\t {round(geneparam[3], 4)}")
        st.sidebar.text(f"Mixing fraction f:\t{round(geneparam[4], 3)}")
        st.sidebar.text(f"Threshold t:\t\t {int(geneparam[5])}")
    elif contdist == "Poisson":
        st.sidebar.text(f"Poisson Param νγ:\t {round(geneparam[0], 4)}")
        st.sidebar.text(f"Mean Expression ρμ:\t {round(geneparam[1], 1)}")
        st.sidebar.text(f"Dispersion α:\t\t {round(geneparam[2], 4)}")
        st.sidebar.text(f"Mixing fraction f:\t{round(geneparam[3], 3)}")
        st.sidebar.text(f"Threshold t:\t\t {int(geneparam[4])}")

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

    if st.sidebar.button("MoM Fit"):
        # Run MoM Script
        script_path = "./BarcodeContScripts/MoM_run.jl"  # Replace with the actual path to your Julia script

        # Run the subprocess with real-time output
        with subprocess.Popen(["julia", script_path, rep, cond, contdist, str(new_thresh), selected_gene], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
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
            
    # if st.sidebar.button("Rerun Fit w/ new threshold"):
    #     # Rerun fit with new initial condition
    #     script_path = "./BarcodeContScripts/FitSingleGene.jl"  # Replace with the actual path to your Julia script

    #     # Run the subprocess with real-time output
    #     with subprocess.Popen(["julia", script_path, rep, cond, contdist, str(new_thresh), selected_gene], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
    #         for line in proc.stdout:
    #             st.write(line.strip())

    #     # Wait for the subprocess to finish and get the return code
    #     return_code = proc.wait()

    #     if return_code == 0:
    #         st.success("Fit Completed")
    #         st.experimental_rerun() 
    #     else:
    #         st.error("The script encountered an error.")
    #         st.write("Error message:")
    #         st.write(proc.stderr)
           
    # if st.sidebar.button("Rerun with fixed f"):
    #     # Run Tuning Script
    #     script_path = "./BarcodeContScripts/FitSingleGene_fixedf.jl"  # Replace with the actual path to your Julia script

    #     # Run the subprocess with real-time output
    #     with subprocess.Popen(["julia", script_path, rep, cond, selected_gene, str(new_f), str(geneparam[0]), str(geneparam[1]), str(geneparam[2]), str(geneparam[3])], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
    #         for line in proc.stdout:
    #             st.write(line.strip())

    #     # Wait for the subprocess to finish and get the return code
    #     return_code = proc.wait()

    #     if return_code == 0:
    #         st.success("Final Fit Saved")
    #         st.experimental_rerun() 
    #     else:
    #         st.error("The script encountered an error.")
    #         st.write("Error message:")
    #         st.write(proc.stderr)
    
    if st.sidebar.button("Reset Fit"):
        try:
            os.remove(mixture_path_final)
        except OSError:
            pass
        st.experimental_rerun()

    if st.sidebar.button("Remove Gene"):
        removed_plasmids.append(selected_gene)
        f = open(removedpath, 'w')
        f.writelines(p+'\n' for p in removed_plasmids)
        f.close()
        st.experimental_rerun()

    ##########################
    # Display Unfinalised Genes
    st.sidebar.markdown("Unchecked Genes:")
    for gene in unfinalised_genes:
        st.sidebar.markdown(f"- {gene}")


    if st.sidebar.button("Label Counts"):
        script_path = "./BarcodeContScripts/LabelCounts.jl"  # Replace with the actual path to your Julia script

        # Run the subprocess with real-time output
        with subprocess.Popen(["julia", script_path, rep, cond, contdist], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
            for line in proc.stdout:
                st.write(line.strip())

        # Wait for the subprocess to finish and get the return code
        return_code = proc.wait()

        if return_code == 0:
            st.success("Labeling Complete!")

        else:
            st.error("The script encountered an error.")
            st.write("Error message:")
            st.write(proc.stderr)