#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 10:30:21 2022

Takes unique reactant SMILES and generates BRICS fragments.

Returns CSV files containing BRICS frags as SMILES and their MHFP (512)
and ECFP4 fingerprints in an array.

Also returns the molecules which couldn't be broken down by the 
BRICS algorithm. - not yet

INPUTS: Reactant SMILES csv file
OUTPUTS: ECFP4_BRICS_frags.csv - cols['Smiles', 'label', fp]
         MHFP_BRICS_frags.csv - cols['Smiles', 'label', fp]
         
@author: George Hodgin
"""
from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
import re
from mhfp.encoder import MHFPEncoder
import tmap

def main():

    # Import csv file and input name for output file
    input_file = input("Enter input file name: ")
    
    fragment_smiles = get_BRICS_frags(input_file)
    ecfp_df = pd.DataFrame(compute_ecfp_descriptors(fragment_smiles))
    ecfp_df.to_csv("ECFP4_BRICS_frags.csv", index=False)
    
    mhfp_df = pd.DataFrame(compute_mhfp_descriptors(fragment_smiles))
    mhfp_df.to_csv("MHFP_BRICS_frags.csv", index=False)
    
def _compute_single_mhfp_descriptor(smiles):
    # The number of permutations used by the MinHashing algorithm
    perm = 512
    
    # Initializing the MHFP encoder with 512 permutations
    enc = MHFPEncoder(perm)
    # Create MHFP fingerprints from SMILES
    # The fingerprint vectors have to be of the tm.VectorUint data type
    # Will import them as numpy arrays for consistency and change type
    fp = tmap.VectorUint(enc.encode(smiles))
    fp_array = np.array(fp)
    return fp_array

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
                   value='BRICS fragment')
    
    return fp_df  


def _prune_BRICS_frags(fragment_smiles):
    
    fragments_df = pd.DataFrame(fragment_smiles, columns=['SMILES'])
    
    # get rid of fragments which couldn't be decomposed
    weird_mols = []
    for i in fragment_smiles:
        z = re.compile("\[\d\*\]|\[\d\d\*\]")
        if not re.match(z, i):
            weird_mols.append(i)
            
    only_frags_df = fragments_df[~fragments_df['SMILES'].isin(weird_mols)]
    
    return only_frags_df['SMILES'].to_list()    

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
                   value='BRICS fragment')
    
    return fp_df

def get_BRICS_frags(input_file):
    
    """Gets the BRICS frags from the csv file input,
        returns them as smiles """
        
    # dataframe to store SMILES and Mols
    data = pd.read_csv(input_file)
    # convert SMILES to mols
    Chem.PandasTools.AddMoleculeColumnToFrame(data,
                                              'Smiles',
                                              'Molecule',
                                              includeFingerprints=False)
    # get rid of any that cant be converted
    data.dropna(axis='rows', inplace=True)

    # Get the BRICS fragments for each molecule
    fragments = []
    for m in data['Molecule']:
        try:
            res_all = list(Chem.BRICS.BRICSDecompose(m, silent=True,
                                                     keepNonLeafNodes=False,
                                                     returnMols=False))
            fragments += res_all
        except:
            pass
    return _prune_BRICS_frags(set(fragments))


    

if __name__ == "__main__":
    main()
