#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 09:57:51 2022

Produces a chemical space plot using tmap, molecules are plotted on x,y 
axis accoring to LSH Forest indexing.

Code adapted from https://tmap.gdb.tools/#lsh-forest and 
https://xinhaoli74.github.io/posts/2020/05/TMAP/

INPUTS: MHFP_BRICS_frags.csv; MHFP_reactants.csv
OUTPUTS: faerun plot to visualise chemical space

@author: gah
"""
import pandas as pd
import numpy as np
from mhfp.encoder import MHFPEncoder
import tmap
from faerun import Faerun

def main():

   
    reactants = input("Enter reactant fps file name: ")
    BRICS_frags = input("Enter BRICS frags fps file name: ")
    
    df = _combine_dataframes(reactants, BRICS_frags)
    fingerprints =  _np_to_vectorUintd(df)
    x, y, s, t = LSH_forest_index(fingerprints)
    
    faerun = Faerun(view="front", coords=False)
    faerun.add_scatter(
        "BRICS_fragments_chemical_space_plot",
        {   "x": x, 
            "y": y, 
            "c": df['label'].to_list(), 
            "labels": df["Smiles"]},
        point_scale=3,
        colormap = 'Set1_r',
        has_legend=True,
        legend_title = 'Molecule / Fragment',
        series_title = 'Reactant vs Fragments',
        legend_labels= [[(0, "Reactant"),
                         (1, "Unique BRICS fragment")]],
        categorical=True
        #shader = 'smoothCircle'
    )
    
    faerun.add_tree("BRICS_fragments_Basic_tree",
                    {"from": s, "to": t},
                    point_helper="BRICS_fragments_chemical_space_plot")
    
    faerun.plot("BRICS_fragments_chemical_space_plot", template="smiles")
    

def _combine_dataframes(infile1, infile2):
    """takes input file names and returns combined dataframe """
    df1 = pd.read_csv(infile1)
    df1['label'] = 0
    df2 = pd.read_csv(infile2)
    df2['label'] = 1
    return pd.concat([df1,df2])
    
def _np_to_vectorUintd(dataframe):
    """ takes the numpy array of a mhfp fingerprint and converts it to
        the tmap datatype tmap.VectorUint for use in LSH forest indexing"""
    array = dataframe.iloc[:,2:].to_numpy()
    fingerprints = [tmap.VectorUint(array[i,:]) for i in range(array.shape[0])]
    return fingerprints

def LSH_forest_index(fingerprints):
    
    #set # permutations
    perm = 512
    
    # Initialize the LSH Forest
    lf1 = tmap.LSHForest(perm)

    # Add the Fingerprints to the LSH Forest and index
    lf1.batch_add(fingerprints)
    lf1.index()

    # Get the coordinates
    x, y, s, t, _ = tmap.layout_from_lsh_forest(lf1)
    return x, y, s, t

if __name__ == "__main__":
    main()