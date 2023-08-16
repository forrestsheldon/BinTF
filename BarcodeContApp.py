# main.py
import os
import streamlit as st
import BarcodeContAppFiles.Intro as Intro
import BarcodeContAppFiles.RawCounts as RawCounts
import BarcodeContAppFiles.InitialLabels as InitialLabels
import BarcodeContAppFiles.ContaminantFitExplorer as ContaminantFitExplorer
import BarcodeContAppFiles.SelectContaminantParameters as SelectContaminantParameters
import BarcodeContAppFiles.CheckFits as CheckFits
import BarcodeContAppFiles.FinalLabels as FinalLabels
import BarcodeContAppFiles.CompareRepsFinal as CompareRepsFinal

def loadlist(fname):
    try:
        with open(fname) as f:
            mylist = [line.rstrip() for line in f]
    except FileNotFoundError:
        mylist = [None]
    return mylist

replist = loadlist("replicates.txt")
condlist = loadlist("condition.txt")
datapathprefix = loadlist("pathinfo.txt")[0]

def main():
    rep = st.sidebar.selectbox("Select a replicate:", replist)
    cond = st.sidebar.selectbox("Select a Condition:", condlist)
    contdist = st.sidebar.selectbox("Select a Contamination Distribution:", ["Full", "Poisson"])

    datadir = os.path.join(datapathprefix, f"{cond}_{''.join(replist)}")


    st.sidebar.title("Analysis Steps")
    app_pages = {
        "Data Requirements": Intro,
        "Initial Labels": InitialLabels,
        "Examine Raw Counts": RawCounts,
        "Select Contamination Fit Genes": ContaminantFitExplorer,
        "Fit Contamination Parameters": SelectContaminantParameters,
        "Check Fits": CheckFits,
        "Final Labels": FinalLabels,
        "Compare Replicates Final": CompareRepsFinal
        }
    selected_app = st.sidebar.selectbox("Select a page:", options=list(app_pages.keys()))

    # Clean up the main area before displaying the selected app
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("---")

    app_pages[selected_app].app(datadir, rep, cond, contdist)

if __name__ == "__main__":
    main()