import os
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import json
import subprocess
from statistics import median
from numpy import arcsinh, cosh, sinh, sqrt

def app(datadir, rep, cond):
    #####################
    # Set Path Variables
    #####################
    fitdir = os.path.join(datadir, "GeneFitData", "contaminantfits")

    #####################
    # Functions
    #####################
    parampath = os.path.join(fitdir, f"{cond}_{rep}_Parameters.tsv")
    genepath = os.path.join(fitdir, f"{cond}_{rep}_SelectedFits.json")
    
    def load_tsv_file(file_path):
        return pd.read_csv(file_path, sep='\t')

    st.title("Calculate Contaminant Parameters")
    st.text("""This will select contaminant parameters based on the fit genes. It uses the median and median
    absolute deviation to calculate the parameter ranges.""")

    #####################
    # Pick Run and Plot
    #####################

    with open(genepath, "r") as f:
        selected_genes = json.load(f)

    param_df = load_tsv_file(parampath)

    nuvalues = param_df.loc[0, selected_genes]
    gamvalues = param_df.loc[1, selected_genes]
    nugamvalues = [n*g for n, g in zip(nuvalues, gamvalues)]
    
    gammedian = median(gamvalues)
    gamMAD = median(abs(g - gammedian) for g in gamvalues)

    nugammedian = median(nugamvalues)
    nugamMAD = median(abs(ng - nugammedian) for ng in nugamvalues)

    nuest = nugammedian/gammedian
    nuMAD = sqrt(1/gammedian**2*nugamMAD**2 + nugammedian**2/gammedian**4*gamMAD**2)

    scat_fig = px.scatter(x = nuvalues, y = gamvalues, title="ν vs. γ")
    scat_fig.update_layout(xaxis_title="ν", yaxis_title="γ")
    scat_fig.update_traces(hovertemplate='<b>%{text}</b><br><br>' +
        'ν: %{x}<br>' +
        'γ: %{y}<br>',
        text=selected_genes)

    # c = st.number_input("Enter the value for constant νγ", value=nugammedian, step=0.001, format="%.3f")


    # hprblparam = [arcsinh((n-g)/(2*sqrt(n*g))) for n,g in zip(nuvalues, gamvalues)]
    # hprblmedian = median(hprblparam)
    # hprblMAD = median(abs(t - hprblmedian) for t in hprblparam)

    # Create x values for the curve
    numin = min(nuvalues)
    numax = max(nuvalues)
    nustep = (numax-numin)/100
    x_curve = [numin+n*nustep for n in range(100)]
    y_curve = [nugammedian / x for x in x_curve]

    # Add the curve to the scatter plot
    scat_fig.add_trace(go.Scatter(x=x_curve, y=y_curve, mode='lines', name='νγ=c'))

    # gamma = (cosh(hprblmedian) - sinh(hprblmedian))*sqrt(nugammedian)
    # nu = nugammedian/gamma

    scat_fig.add_trace(go.Scatter(x=[nuest], y=[gammedian], mode='markers', marker=dict(color='red',), name='Estimate'))

    log_scale = st.checkbox("Switch axes to log scale")
    if log_scale:
        scat_fig.update_xaxes(type="log")
        scat_fig.update_yaxes(type="log")

    st.write(scat_fig)


    st.text(f"Volume Ratio ν:\t\t {round(nuest, 4)} +- {round(1.48*nuMAD, 4)}")
    st.text(f"Burst Rate γ:\t\t {round(gammedian, 4)} +- {round(1.48*gamMAD, 4)}")

    #####################
    # Save Contaminant Parameters
    #####################
    st.text("Save the parameters. These will be used to refit the remaining genes. The fits will be regularised by the values of γ and νγ  from the contaminant fits.")
    
    if st.button("Save Contaminant parameters"):
        with open(os.path.join(fitdir, f"{cond}_{rep}_ContaminantParameters.txt"), "w") as f:
            f.write(f"{gammedian}\n")
            f.write(f"{gamMAD}\n")
            f.write(f"{nugammedian}\n")
            f.write(f"{nugamMAD}")
        st.success("The parameters have been saved.")

    st.title("Run Full Fit")

    if st.button("Full Gene Fit: will take a while"):
        script_path = "./GeneBinScripts/FitAllGenes.jl"
        
        with subprocess.Popen(["julia", script_path, rep, cond], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) as proc:
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




