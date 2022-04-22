#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 13:36:21 2022

Takes the cleaned reaxys csv data and splits up reaction SMILES
into products and reactants.

Returns the unique reactants as a .csv file

@author: gah
"""
import pandas as pd

def main():
    
    infile = input("Please enter the name of the input file")
    outfile = input("Please enter the name of the output file")
    
    get_smiles(infile).to_csv(outfile)    
    
def get_smiles(infile):
    
    data = pd.read_csv(infile)
    
    # Put in reactants as individual column
    data['Reactants'] = data['Reaction'].str.split(">>",
                                                         n=1,
                                                         expand=True)[0]
    # Split the reactants into individual molecules
    split_smiles = [i.split('.')for i in data['Reactants']]
    
    # combine all the molecules and remove duplicates
    set_smiles = {n for row in split_smiles for n in row}
    
    unique_smiles = pd.DataFrame(set_smiles, columns=['Smiles'])
    
    return unique_smiles

if __name__ == "__main__":
    main()