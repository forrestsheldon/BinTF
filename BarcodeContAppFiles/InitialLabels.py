import os
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import numpy as np

def load_data(datapath):
    # Load the TSV file
    df = pd.read_csv(datapath, sep='\t')
    return df

def count_nonzeros(row, selected_plasmids):
    # Count non-zero entries for selected plasmids
    return np.count_nonzero(row[selected_plasmids])

def create_histogram(df, selected_plasmids, numempty):
    # Add a new column 'nonzero_count' to the DataFrame
    df['nonzero_count'] = df.apply(count_nonzeros, axis=1, args=(selected_plasmids,))

    # Create a new array that includes the non-zero counts and the specified number of zeros
    counts_with_zeros = np.concatenate([df['nonzero_count'].values, np.zeros(numempty)])

    # Create a histogram using Plotly
    fig = go.Figure(data=[go.Histogram(x=counts_with_zeros)])

    # Calculate and print the average
    avg = counts_with_zeros.mean()
    st.write(f"Mean Unique Plasmids from SCseq: {avg:.2f}")

    # Add a vertical line for the average
    fig.add_shape(type='line', 
                  line=dict(dash='dash'),
                  x0=avg, x1=avg, y0=0, y1=1, yref='paper')

    # Update layout
    fig.update_layout(xaxis_title='Number of Unique Plasmids',
                      yaxis_title='Number of Cells',
                      shapes=[dict(type= 'line',
                                   yref= 'paper', y0= 0, y1= 1,
                                   xref= 'x', x0= avg, x1= avg)])

    return fig


def app(datadir, rep, cond, contdist):

    # Use the title
    st.title('Naive Labeling')

    # Use a sidebar for user inputs
    datapath = os.path.join(datadir, f"{cond}_{rep}_CellCounts.tsv")

    with open(os.path.join(datadir, f"{rep}_ddPCR.txt"), 'r') as file:
        ddPCR = float(file.read())

    with open(os.path.join(datadir, f"{cond}_{rep}_noBCcells.txt"), 'r') as file:
        noBCcells = int(file.read())
    
    st.write(f"Mean number of insertions from ddPCR: {ddPCR}\n")

    # Load data
    df = load_data(datapath)

    

    # Get list of plasmids
    plasmid_cols = df.columns[1:]

    # Create histogram
    fig = create_histogram(df, plasmid_cols, noBCcells)

    # Display plot
    st.plotly_chart(fig)

    st.text("The number of unique plasmids from SCseq is much higher than the number\n"+
            "found by ddPCR, across replicates and sequencing preparations. The line\n"+
            "marks the average number of unque plasmids per cell")
