# -*- coding: utf-8 -*-
"""
Created on Thu Nov 27 19:39:50 2025

@author: tblasco
"""

import os
import pandas as pd
import time
from multiprocessing import freeze_support
from persephone import runMgPipe

def main():
    # Load data
    counts = pd.read_table("/Users/josefinaarcagni/Downloads/persephonepy_toJose/counts.tsv", sep='\t', index_col=0)
    taxa = pd.read_table("/Users/josefinaarcagni/Downloads/persephonepy_toJose/taxa.tsv", sep='\t', index_col=0)

    # Select some diets
    diets = ['EUAverageDiet','HighFiberDiet','HighProteinDiet','VegetarianDiet','UnhealthyDiet']

    # Create microbiome community models and calculate fluxes on different diets
    for diet in diets:
        # Run function
        tic = time.time()
        rxnAb, rxnPre, subsAb = runMgPipe(
            counts=counts,
            taxa=taxa,
            compute_profiles=True,
            solver='highs',
            output_dir=os.path.join('.', 'output'),
            diet_file_name=diet
        )
        print(f'Elapsed time: {time.time()-tic} seconds.')

        # Rename files
        os.rename(os.path.join('.', 'output', 'netSecretionFluxes.csv'),
                  os.path.join('.', 'output', f'netSecretionFluxes_{diet}.csv'))
        os.rename(os.path.join('.', 'output', 'netUptakeFluxes.csv'),
                  os.path.join('.', 'output', f'netUptakeFluxes_{diet}.csv'))
        


    # Net secretion fluxes indicate the amount of each metabolite that can be produced

    # Heatmap of subsystem abundance ?

    # Run differential analysis to check differences between diets

    # Boxplots of some significant metabolites        

if __name__ == "__main__":
    freeze_support()  # ensures multiprocessing works on macOS and Windows
    main()
