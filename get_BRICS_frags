#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 10:30:21 2022

Takes unique reactant SMILES and generates BRICS fragments.

Returns CSV files containing BRICS frags as SMILES and their MHFP and ECFP4
fingerprints in an array. -just ecfp4 in this version.

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


def main():

    # Import csv file and input name for output file
    input_file = input("Enter input file name")
    outfile = input("Enter output file name")
    fragment_smiles = get_BRICS_frags(input_file)
    ecfp_df = pd.DataFrame(compute_ecfp_descriptors(fragment_smiles))
    ecfp_df.insert(loc=0,
                   column='Smiles',
                   value=fragment_smiles)
    
    ecfp_df.insert(loc=1,
                   column='label',
                   value='BRICS fragment')
    
    ecfp.to_csv(outfile, index=False)

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

    return np.vstack(descriptors), keep_idx

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
    return fragments


    

if __name__ == "__main__":
    main()
