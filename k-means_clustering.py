#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:17:11 2022

Does k-means clustering analysis of reactant molecule and BRICS fragment
SMILES to determine the chemical similarity of the two sets.



INPUTS: ECFP4_BRICS_frags.csv; ECFP4_reactants.csv
OUTPUTS: csv files of reactant and fragment counts and % of total count in
         each cluster.

@author: gah
"""
import numpy as np
import pandas as pd
from sklearn.cluster import MiniBatchKMeans
import math
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt


def main():

    num_clusters = int(input("number of clusters: "))
    reactants = input("Enter reactant ecfp file name: ")
    BRICS_frags = input("Enter BRICS frags ecfp file name: ")

    data = _combine_dataframes(reactants, BRICS_frags)

    arr = np.array(data.values[0::, 2::], dtype=np.float16)
    km = MiniBatchKMeans(n_clusters=num_clusters,
                         random_state=0,
                         batch_size=3 * num_clusters)
    km.fit(arr)

    chunk_size = 500
    all_data = np.array(data.values[0::, 2::], dtype=bool)
    chunks = math.ceil(all_data.shape[0] / chunk_size)

    # It looks like the predict method chokes if you send too much data,
    # chunking to 500 seems to work
    cluster_id_list = []
    for row, names in zip(np.array_split(all_data, chunks),
                          np.array_split(data['Smiles'].values,
                                         chunks)):
        p = km.predict(row)
        cluster_id_list += list(p)
    data.insert(2, "Cluster", cluster_id_list)
    center_list = find_cluster_centers(data, km.cluster_centers_)
    data.insert(3, "Center", center_list)
    out_df = data[["Smiles", "Cluster", "Center", "label"]]

    count_df = pd.DataFrame(out_df.groupby(['label', 'Cluster']).size())
    count_df.reset_index(inplace=True)
    count_df.rename(columns={0: "Count"}, inplace=True)

    frag_counts_df = count_df[count_df["label"] == 1]
    frag_counts_df['% of all frags'] = (frag_counts_df['Count']
                                        / frag_counts_df['Count'].sum())*100
    reactant_counts_df = count_df[count_df["label"] == 0]
    reactant_counts_df['% of all reactants'] = (reactant_counts_df['Count']
                                               / reactant_counts_df['Count'].sum())*100
    
    frag_counts_df.plot(kind="bar",
                    x="Cluster",
                    y="% of all frags",
                    title="Fragment distribution",
                    ylabel="% of total fragments")
    plt.savefig("k-means_30_frag_plot.png")
    
    reactant_counts_df.plot(kind="bar",
                    x="Cluster",
                    y="% of all reactants",
                    title="Reactant distribution",
                    ylabel="% of total reactants")
    plt.savefig("k-means_30_reactant_plot.png")
    
    return( frag_counts_df.to_csv(
        'frag_counts_k-means_c{}.csv'
                                 .format(
                                     str(
                                     num_clusters)), index=False), reactant_counts_df.to_csv(
                                         'reactant_counts_k-means_c{}.csv'
                                         .format(
                                             str(
                                             num_clusters)), index=False))


def _combine_dataframes(infile1, infile2):
    """takes input file names and returns combined dataframe """
    df1=pd.read_csv(infile1)
    df1['label']=0
    df2=pd.read_csv(infile2)
    df2['label']=1
    return pd.concat([df1, df2])


def find_cluster_centers(df, centers):
    center_set=set()
    for k, v in df.groupby("Cluster"):
        fp_list=v.values[0::, 3::]
        XA=np.array([centers[k]]).astype(float)
        XB=np.array(fp_list).astype(float)
        dist_list=cdist(XA, XB)
        min_idx=np.argmin([dist_list])
        center_set.add(v['Smiles'].values[min_idx])
    return ["Yes" if x in center_set else "No" for x in df['Smiles'].values]

if __name__ == "__main__":
    main()
