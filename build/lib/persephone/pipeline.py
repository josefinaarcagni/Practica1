# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 16:28:45 2025

@author: tblasco
"""

import json
import pandas as pd
import os
import time

from .MARS import runMARS
from .MgPipe import runMgPipe
from .pWBM import persWBM
from .WBM import generateMWBM, analyzeMWBM

def runPersephone(config_file: str) -> (pd.DataFrame):
    ###########################################################################
    # 0 General settings
    ###########################################################################
    # Read configuration file
    with open(config_file, 'r', encoding = 'utf-8') as f:
        config = json.load(f)
        
    # Load metadata
    meta = pd.read_table(config['general_metadataPath'], sep = '\t', index_col = 0)
    
    # Read flags
    flagMARS = config['mars_flagMars']
    flagMicro = config['micro_flagMicro']
    flagPersWBM = config['pers_flagPersWBM']
    flagMWBM = config['mwbm_flagMWBM']
    flagFBA = config['fba_flagFBA']
    
    # Initialize log file
    timeLog = pd.DataFrame(0, index = ['MARS','MgPipe','WBM personalisation','mWBM generation','FBA'], columns = ['Time'])
    
    ###########################################################################
    # 1.0 SeqC
    ###########################################################################
    
    ###########################################################################
    # 1.5 MARS
    ###########################################################################
    tic = time.time()
    if flagMARS:
        # Define output directory
        config['mars_outputDir'] = os.path.join(config['results_path'], 'MARS')
        
        # Read input data 
        counts = pd.read_table(config['mars_readTablePath'], sep = '\t', index_col = 0)
        taxa = pd.read_table(config['mars_taxaTablePath'], sep = '\t', index_col = 0)
        
        # Run MARS
        countsMapped, taxaMapped = runMARS(counts = counts, taxa = taxa, database = config['general_database'], output_format = config['mars_outputExtensionMars'],
                                           taxa_delimiter = config['mars_taxaDelimiter'], flag_lone_species = config['mars_flagLoneSpecies'],
                                           cutoff = config['mars_cutoff'], sample_read_counts_cutoff = config['mars_sampleCutoff'],
                                           remove_clades = config['mars_removeClades'], output_dir = config['mars_outputDir'])
        
    timeLog.loc['MARS','Time'] = time.time() - tic
    ###########################################################################
    # 2 MgPipe
    ###########################################################################
    tic = time.time()
    if flagMicro:
        # Define output directory
        config['micro_outputDir'] = os.path.join(config['results_path'], 'MgPipe')
        
        # Read input data (if MARS not executed)
        if not flagMARS:
            print('TO DO')
            
        # Run MgPipe
        rxnAb, rxnPre, subsAb = runMgPipe(counts = countsMapped, taxa = taxaMapped, database = config['general_database'],
                                          compute_profiles = config['micro_computeProfiles'], solver = config['general_solver'], threads = config['general_threads'],
                                          diet_file_name = config['general_diet'], output_dir = config['micro_outputDir'])
        
    timeLog.loc['MgPipe','Time'] = time.time() - tic
    ###########################################################################
    # 3 WBM personalisation
    ###########################################################################
    tic = time.time()
    if flagPersWBM:
        # Define output directory
        config['pers_outputDir'] = os.path.join(config['results_path'], 'persWBM')
        
        # Run WBM personalisation
        persWBM(meta = meta, params_file_path = config['pers_physiologicalParamFile'], mets_file_path = config['pers_metaboliteConcentrationFile'],
                pers_wbm_path = config['pers_outputDir'], solver = config['general_solver'], threads = config['general_threads'])
        
    timeLog.loc['WBM personalisation','Time'] = time.time() - tic
    ###########################################################################
    # 4 Community microbiome WBM creation
    ###########################################################################
    tic = time.time()
    if flagMWBM:
        # Define output directory
        config['mwbm_outputDir'] = os.path.join(config['results_path'], 'mWBM')
        
        # Run mWBM creation
        if config['mwbm_usePersonalisedWBM']:
            stats = generateMWBM(meta = meta, micro_mod_path = os.path.join(config['micro_outputDir'], 'models'), diet_file_name = config['general_diet'],
                                 use_pers_wbm = config['mwbm_usePersonalisedWBM'], pers_wbm_path = os.path.join(config['pers_outputDir'], 'models'), 
                                 output_dir = config['mwbm_outputDir'], solver = config['general_solver'], threads = config['general_threads'])
        else:
            stats = generateMWBM(meta = meta, micro_mod_path = os.path.join(config['micro_outputDir'], 'models'), diet_file_name = config['general_diet'],
                                 output_dir = config['mwbm_outputDir'], solver = config['general_solver'], threads = config['general_threads'])
    
    timeLog.loc['mWBM generation'] = time.time() - tic
    ###########################################################################
    # 5 FBA
    ###########################################################################
    tic = time.time()
    if flagFBA:
        # Define output directory
        config['fba_outputDir'] = os.path.join(config['results_path'], 'FBA')
        
        # Run FBA
        res = analyzeMWBM(meta = meta, rxn_list_path = config['fba_rxnListFilePath'], mwbm_path = os.path.join(config['mwbm_outputDir'],'models'),
                          output_dir = config['fba_outputDir'], solver = config['general_solver'], threads = config['general_threads'])
    
    timeLog.loc['FBA','Time'] = time.time() - tic
    ###########################################################################
    # 6 Statistical analysis
    ###########################################################################
    
    ###########################################################################
    # Output variables
    ###########################################################################
    return(timeLog)