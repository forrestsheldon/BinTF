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
    fitdir = os.path.join(datadir, "GeneFitData", "noisefits")

    #####################
    # Functions
    #####################
    parampath = os.path.join(fitdir, f"{cond}_{rep}_Parameters.tsv")
    genepath = os.path.join(fitdir, f"{cond}_{rep}_SelectedFits.json")
    
    def load_tsv_file(file_path):
        return pd.read_csv(file_path, sep='\t')

    st.title("Set Noise Parameter")
    st.text("""This will select noise parameters based on the fit genes. It uses the median and median
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
    
    nugammedian = median(nugamvalues)
    nugamMAD = median(abs(ng - nugammedian) for ng in nugamvalues)

    scat_fig = px.scatter(x = nuvalues, y = gamvalues, title="ν vs. γ")
    scat_fig.update_layout(xaxis_title="ν", yaxis_title="γ")
    scat_fig.update_traces(hovertemplate='<b>%{text}</b><br><br>' +
        'ν: %{x}<br>' +
        'γ: %{y}<br>',
        text=selected_genes)

    # c = st.number_input("Enter the value for constant νγ", value=nugammedian, step=0.001, format="%.3f")


    hprblparam = [arcsinh((n-g)/(2*sqrt(n*g))) for n,g in zip(nuvalues, gamvalues)]
    hprblmedian = median(hprblparam)
    hprblMAD = median(abs(t - hprblmedian) for t in hprblparam)

    # Create x values for the curve
    numin = min(nuvalues)
    numax = max(nuvalues)
    nustep = (numax-numin)/100
    x_curve = [numin+n*nustep for n in range(100)]
    y_curve = [nugammedian / x for x in x_curve]

    # Add the curve to the scatter plot
    scat_fig.add_trace(go.Scatter(x=x_curve, y=y_curve, mode='lines', name='νγ=c'))

    gamma = (cosh(hprblmedian) - sinh(hprblmedian))*sqrt(nugammedian)
    nu = nugammedian/gamma

    scat_fig.add_trace(go.Scatter(x=[nu], y=[gamma], mode='markers', marker=dict(color='red',), name='Estimate'))

    log_scale = st.checkbox("Switch axes to log scale")
    if log_scale:
        scat_fig.update_xaxes(type="log")
        scat_fig.update_yaxes(type="log")

    st.write(scat_fig)

    st.title("Values Along Histogram")

    st.text("""The values from the above histogram are binned along the arc length of the curve.
The median is marked by the line.""")

    # Create the histogram
    histfig = px.histogram(hprblparam, nbins=20)

    histfig.add_shape(
        type='line',
        x0=hprblmedian, y0=0,
        x1=hprblmedian, y1=1,
        yref='paper',
        line=dict(color='red')
        )

    # Change x-axis label
    histfig.update_xaxes(title_text='Hyperbola Parameter')

    # Remove the legend
    histfig.update_layout(showlegend=False)

    # Display the plot in Streamlit
    st.plotly_chart(histfig)

    st.text(f"Volume Ratio ν:\t\t {round(nu, 4)}")
    st.text(f"Burst Rate γ:\t\t {round(gamma, 4)}")

    #####################
    # Save Noise Parameters
    #####################
    st.text("Save the noise parameters. This will be used to refit the remaining genes.")
    
    if st.button("Save Noise parameters"):
        with open(os.path.join(fitdir, f"{cond}_{rep}_NoiseParameters.txt"), "w") as f:
            f.write(f"{nu}\n")
            f.write(f"{gamma}")
        with open(os.path.join(fitdir, f"{cond}_{rep}_NoiseParameters2.txt"), "w") as f:
            f.write(f"{nugammedian}\n")
            f.write(f"{nugamMAD}\n")
            f.write(f"{hprblmedian}\n")
            f.write(f"{hprblMAD}")
        st.success("The noise parameters have been saved.")

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




