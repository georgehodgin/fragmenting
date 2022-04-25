#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 09:25:33 2022
Takes SMILES codes for unique reactant molecules from OPRD dataset and
computes their ecfp4 and mhfp (512) fingerprints for use in plotting chemical
space maps with tmap.

INPUT: unique_smiles.csv
OUTPUT: ECFP_reactants.csv; MHFP_reactants.csv

@author: gah
"""
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
from mhfp.encoder import MHFPEncoder
import tmap

def main():

    # Import csv file and input name for output file
    input_file = input("Enter input file name: ")
    
    df = pd.read_csv(input_file)
    reactant_smiles = df['Smiles'].to_list()
    
    ecfp_df = pd.DataFrame(compute_ecfp_descriptors(reactant_smiles))
    ecfp_df.to_csv("ECFP4_reactants.csv", index=False)
    
    mhfp_df = pd.DataFrame(compute_mhfp_descriptors(reactant_smiles))
    mhfp_df.to_csv("MHFP_reactants.csv", index=False)
    
def _compute_single_mhfp_descriptor(smiles):
    # The number of permutations used by the MinHashing algorithm
    perm = 512
    
    # Initializing the MHFP encoder with 512 permutations
    enc = MHFPEncoder(perm)
    # Create MHFP fingerprints from SMILES
    # The fingerprint vectors have to be of the tm.VectorUint data type
    # Will import them as numpy arrays for consistency and change type
    try:
        fp = tmap.VectorUint(enc.encode(smiles))
        fp_array = np.array(fp)
        return fp_array
    except Exception as E:
        return None
    
def compute_mhfp_descriptors(smiles_list):
    
    """ Computes ecfp descriptors  for list of SMILES"""

    keep_idx = []
    descriptors = []
    for i, smiles in enumerate(smiles_list):
        mhfp = _compute_single_mhfp_descriptor(smiles)
        if mhfp is not None:
            keep_idx.append(i)
            descriptors.append(mhfp)
    df = pd.DataFrame(smiles_list)
    kept_frags = df.iloc[keep_idx]
    fp_df = pd.DataFrame(np.vstack(descriptors))
    fp_df.insert(loc=0,
                 column='Smiles',
                 value=kept_frags)
    fp_df.insert(loc=1,
                   column='label',
                   value='Reactant')
    
    return fp_df  

def _compute_single_ecfp_descriptor(smiles):
    
    """ Computes the extended connectivity fingerprint for one molecule
        and returns an array of the bit vector """
    try:
        mol = Chem.MolFromSmiles(smiles)
    except Exception as E:
        return None

    if mol:
        fp = Chem.AllChem.GetMorganFingerprintAsBitVect(mol,
                                                        radius=2,
                                                        nBits=2048)
        return np.array(fp)

    return None

def compute_ecfp_descriptors(smiles_list):
    
    """ Computes ecfp descriptors  for list of SMILES"""

    keep_idx = []
    descriptors = []
    for i, smiles in enumerate(smiles_list):
        ecfp = _compute_single_ecfp_descriptor(smiles)
        if ecfp is not None:
            keep_idx.append(i)
            descriptors.append(ecfp)
    df = pd.DataFrame(smiles_list)
    kept_frags = df.iloc[keep_idx]
    fp_df = pd.DataFrame(np.vstack(descriptors))
    fp_df.insert(loc=0,
                 column='Smiles',
                 value=kept_frags)
    fp_df.insert(loc=1,
                   column='label',
                   value='Reactant')
    
    return fp_df

if __name__ == "__main__":
    main()
