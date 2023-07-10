import os
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import numpy as np


def load_tsv(datapath):
    # Load the TSV file
    df = pd.read_csv(datapath, sep='\t')
    return df

def loadlist(fname):
    with open(fname) as f:
        replist = [line.rstrip() for line in f]
    return replist

def app(datadir, rep, cond):
    st.title('Compare Replicates')
    replist = loadlist("replicates.txt")

    fitdir = os.path.join(datadir, "GeneFitData", "noisefits")

    plasmid_names = []
    transfected_fraction = []
    for rep in replist:
        param_df = load_tsv(os.path.join(fitdir, f"{cond}_{rep}_Parameters.tsv"))
        cell_df = load_tsv(os.path.join(datadir, f"{cond}_{rep}_CellCounts.tsv"))
        repnames = []
        fractions = []
        for plasmid in param_df.columns:
            threshold = param_df[plasmid].iloc[5]
            f = (cell_df[plasmid] >= threshold).mean()
            repnames.append(plasmid)
            fractions.append(f)
        plasmid_names.append(repnames)
        transfected_fraction.append(fractions)

    # create a figure with two bar charts
    fig = go.Figure()

    # add bar chart for nonzero_fraction of the first replicate
    fig.add_trace(go.Bar(
        x=plasmid_names[0],
        y=transfected_fraction[0],
        name=f'Fraction of cells labeled ({replist[0]})',
        marker_color='blue'
    ))

    # add bar chart for nonzero_fraction of the second replicate
    fig.add_trace(go.Bar(
        x=plasmid_names[1],
        y=transfected_fraction[1],
        name=f'Fraction of cells labeled ({replist[1]})',
        marker_color='orange'
    ))

    # add a layout
    fig.update_layout(
        title='Comparison of Labeled Fractions between Replicates',
        xaxis_title='Plasmids',
        yaxis_title='Fraction',
        barmode='group',  # bars are placed beside each other
        width=1200,  # width in pixels
        height=600,  # height in pixels
        xaxis_tickangle=-45  # rotate labels by 45 degrees
    )

    st.plotly_chart(fig)
    
    rep1fractions = dict(zip(plasmid_names[0], transfected_fraction[0]))
    rep2fractions = dict(zip(plasmid_names[1], transfected_fraction[1]))
    plasmidsall = sorted(set(rep1fractions.keys()).union(set(rep2fractions.keys())))
    
    nonzerofractionsall = [[], []]
    for plasmid in plasmidsall:
        nonzerofractionsall[0].append(rep1fractions.get(plasmid, 0))
        nonzerofractionsall[1].append(rep2fractions.get(plasmid, 0))

    # create a figure with a scatter plot
    fig2 = go.Figure()

    # add scatter plot for fractions of nonzero entries
    fig2.add_trace(go.Scatter(
        x=nonzerofractionsall[0],
        y=nonzerofractionsall[1],
        mode='markers',
        name='Fraction of cells labeled',
        marker_color='blue',
        hovertemplate='<b>%{text}</b><br><br>' +
        'f1: %{x}<br>' +
        'f2: %{y}<br>',
        text=plasmidsall
    ))


    # add a layout
    fig2.update_layout(
        title='Comparison of Labeled Fractions between Replicates',
        xaxis_title=f'Fraction of cells labeled ({replist[0]})',
        yaxis_title=f'Fraction of cells labeled ({replist[1]})',
        width=900,  # width in pixels
        height=600  # height in pixels
    )

    st.plotly_chart(fig2)
        



