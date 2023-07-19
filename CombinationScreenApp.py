# main.py
import os
import streamlit as st

def loadlist(fname):
    with open(fname) as f:
        replist = [line.rstrip() for line in f]
    return replist

replist = loadlist("replicates.txt")
markerlist = loadlist("markers.txt")
screenident, GEXident, datapath = loadlist("pathinfo.txt")

DATA_DIR_prefix = os.path.join(datapath, f"{screenident}_{''.join(replist)}")

def main():
    st.sidebar.title("Analysis Steps")
    app_pages = {
        "Basic Combinations": BasicCombinations,

        }
    selected_app = st.sidebar.selectbox("Select a page:", options=list(app_pages.keys()))

    # Clean up the main area before displaying the selected app
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown("---")

    app_pages[selected_app].app(datapath, screenident, replist)

if __name__ == "__main__":
    main()