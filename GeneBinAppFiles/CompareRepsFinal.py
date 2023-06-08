import os
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import numpy as np

def load_data(datapath):
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

    # Load data for both replicates and calculate fraction of nonzero entries
    nonzero_fractions_final = []
    plasmid_names = []
    # nonzero_fractions_initial = []
    for rep in replist:
        datapathfinal = os.path.join(datadir, f"{cond}_{rep}_Labeled.tsv")
        # datapathinitial = os.path.join(datadir, f"{cond}_{rep}_CellCounts.tsv")
        dffinal = load_data(datapathfinal)
        # dfinitial = load_data(datapathinitial)
        
        plasmid_cols_final = dffinal.columns[1:]
        plasmid_names.append(plasmid_cols_final)
        # plasmid_cols_initial = dfinitial.columns[1:]
        nonzero_fraction_final = dffinal[plasmid_cols_final].mean()
        nonzero_fractions_final.append(nonzero_fraction_final)
        # nonzero_fraction_initial = dfinitial[plasmid_cols_initial].apply(lambda col: (col != 0).sum() / len(col), axis=0)
        # nonzero_fractions_initial.append(nonzero_fraction_initial)


    # create a figure with two bar charts
    fig = go.Figure()

    # add bar chart for nonzero_fraction of the first replicate
    fig.add_trace(go.Bar(
        x=plasmid_names[0],
        y=nonzero_fractions_final[0],
        name=f'Fraction of cells labeled ({replist[0]})',
        marker_color='blue'
    ))

    # add bar chart for nonzero_fraction of the second replicate
    fig.add_trace(go.Bar(
        x=plasmid_names[1],
        y=nonzero_fractions_final[1],
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
    
    rep1fractions = dict(zip(plasmid_names[0], nonzero_fractions_final[0]))
    rep2fractions = dict(zip(plasmid_names[1], nonzero_fractions_final[1]))
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
    # fig2.add_trace(go.Scatter(
    #     x=nonzero_fractions_initial[0],
    #     y=nonzero_fractions_initial[1],
    #     mode='markers',
    #     name='Fraction of cells labeled',
    #     marker_color='blue',
    #     marker_opacity=0.2
    # ))


    # add a layout
    fig2.update_layout(
        title='Comparison of Labeled Fractions between Replicates',
        xaxis_title=f'Fraction of cells labeled ({replist[0]})',
        yaxis_title=f'Fraction of cells labeled ({replist[1]})',
        width=900,  # width in pixels
        height=600  # height in pixels
    )

    st.plotly_chart(fig2)


