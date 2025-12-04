# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 14:47:09 2025

@author: tblasco
"""

import os
import pandas as pd
import time

from persephone import runMARS, runMgPipe, persWBM, generateMWBM, analyzeMWBM, runPersephone

###############################################################################
# Running entire pipeline
###############################################################################

# Define configuration file path
configFile = 'test/configPersephone.json'

# Run persephone
timeLog = runPersephone(configFile)

###############################################################################
# Running submodules
###############################################################################

# Load data
counts = pd.read_table("/Users/josefinaarcagni/Downloads/persephonepy_toJose/counts.tsv", sep = '\t', index_col = 0)
taxa = pd.read_table("/Users/josefinaarcagni/Downloads/persephonepy_toJose/taxa.tsv", sep = '\t', index_col = 0)
meta = pd.read_table(os.path.join('.','test','metadata.tsv'), sep = '\t', index_col = 0)

# Run MARS module
tic = time.time()
countsMapped, taxaMapped = runMARS(counts = counts, taxa = taxa, database = 'AGORA2-APOLLO', output_dir = os.path.join('.','test','output','MARS'))
time_MARS = time.time() - tic
print(f'MARS computing time: {time_MARS}')

# Run MgPipe module
tic = time.time()
rxnAb, rxnPre, subsAb = runMgPipe(counts = countsMapped, taxa = taxaMapped, database = 'AGORA2-APOLLO', compute_profiles = True, 
                                  diet_file_name = 'EUAverageDiet',
                                  solver = 'cplex', threads = 8, output_dir = os.path.join('.','test','output','MgPipe'))
time_MgPipe = time.time()-tic
print(f'MgPipe computing time: {time_MgPipe}')

# Run WBM personalization
tic = time.time()
paramsConfig = persWBM(meta = meta, params_file_path = os.path.join('.','test','physParams.txt'), mets_file_path = os.path.join('.','test','metCon.txt'),
                       solver = 'cplex', threads = 8, pers_wbm_path = os.path.join('.','test','output','persWBM'))
time_persWBM = time.time()-tic
print(f'persWBM computing time: {time_persWBM}')

# Generate mWBM
tic = time.time()
stats = generateMWBM(meta = meta, micro_mod_path = os.path.join('.','test','output','MgPipe','models'), 
                     diet_file_name = 'EUAVerageDiet', solver = 'cplex', threads = 8,
                     output_dir = os.path.join('.','test','output','mWBM'))
time_MWBM = time.time()-tic
print(f'mWBM computing time: {time_MWBM}')

# Generate mWBM from pWBM
tic = time.time()
stats = generateMWBM(meta = meta, micro_mod_path = os.path.join('.','test','output','MgPipe','models'), 
                     diet_file_name = 'EUAverageDiet', solver = 'cplex', threads = 8,
                     use_pers_wbm = True, pers_wbm_path = os.path.join('.','test','output','persWBM','models'),
                     output_dir = os.path.join('.','test','output','mWBM_pWBM'))
time_MWBM2 = time.time()-tic
print(f'mWBM from pWBM computing time: {time_MWBM2}')

# Apply FBA
tic = time.time()
sol = analyzeMWBM(meta = meta, rxn_list_path = os.path.join('.','test','objRxns.tsv'), mwbm_path = os.path.join('.','test','output','mWBM','models'),
                  output_dir = os.path.join('.','test','output','FBA'), solver = 'cplex', threads = 8)
time_FBA = time.time()-tic
print(f'FBA computing time: {time_FBA}')

