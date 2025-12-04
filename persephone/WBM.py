# -*- coding: utf-8 -*-
"""
Created on Tue Aug 26 09:27:26 2025

@author: tblasco
"""

import os
import re
import math
import numpy as np
import pandas as pd
import pickle
import copy
import importlib.resources as pkg_resources
from functools import partial
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import cobra

from .utils import _buildModel, _addCouplingConstraintsWBM, _addInfoToModel, _addCouplingConstraints, _setSolverParameters, _setDietConstraints

def _generateOriginalWBM(name: str, solver: str, threads: int = 0) -> (cobra.Model):
    # Load model information
    dataPath = pkg_resources.files('persephone').joinpath('data/models/WBM')
    metInfo = pd.read_csv(dataPath.joinpath(name + '.metInfo.csv'))
    rxnInfo = pd.read_csv(dataPath.joinpath(name + '.rxnInfo.csv'))
    constraints = pd.read_csv(dataPath.joinpath(name + '.consInfo.csv'))
    
    # Build model
    model = _buildModel(name, metInfo, rxnInfo)
    
    # Add coupling constraints
    print(f'Adding coupling constraints to sample {name}')
    model = _addCouplingConstraintsWBM(model, constraints)
    
    # Setting solver parameters
    model = _setSolverParameters(model, solver, threads)
    
    # Optimize model
    model.optimize()
    
    # Output variable
    return (model)

def _combineHarveyMicrobiota(name: str, wbm_model: cobra.Model, wbm_mets: pd.DataFrame, wbm_rxns: pd.DataFrame, micro_mets: pd.DataFrame,
                             micro_rxns: pd.DataFrame) -> (cobra.Model, pd.DataFrame, pd.DataFrame):
    # Remove exchange reactions:
    # 1. Diet exchange: 'Diet_EX_met_d': met_d <=>
    # 2. Fecal exchange: 'EX_met_fe': met_fe <=>
    micro_rxns = micro_rxns.loc[~(micro_rxns['rxnID'].str.startswith('Diet_EX_') | micro_rxns['rxnID'].str.startswith('EX_')),:].reset_index(drop = True)
    
    # Convert diet transport reactions.
    # 'DUt_met': 'met_d' -> 'met_u' ==> 'Micro_EX_met_luLI_luM': 'met_luLI' -> 'met_luM'
    ExR = micro_rxns['rxnID'].str.startswith('DUt_')
    micro_rxns.loc[ExR, 'rxnID'] = [re.sub(r'^DUt_','Micro_EX_', r) + '_luLI_luM' for r in micro_rxns.loc[ExR,'rxnID']]
    micro_rxns.loc[ExR, 'rxnName'] = micro_rxns.loc[ExR, 'rxnID']
    micro_rxns['rxnEq'] = micro_rxns['rxnEq'].str.rstrip()
    micro_rxns.loc[:,'rxnEq'] = [re.sub(r'_d;','_luLI;', r) for r in micro_rxns.loc[:,'rxnEq']]
    micro_rxns.loc[:,'rxnEq'] = [re.sub(r'_u;','_luM;', r) for r in micro_rxns.loc[:,'rxnEq']]
    micro_rxns.loc[:,'rxnEq'] = [re.sub(r'_u$','_luM', r) for r in micro_rxns.loc[:,'rxnEq']]
    
    # Make those reactions reversible
    micro_rxns.loc[ExR, 'lb'] = -1000
    micro_rxns.loc[ExR, 'ub'] = 1000
    
    # Get all _d metabolites
    
    # Rename _d and _u metabolites
    micro_mets['metID'] = [re.sub(r'_d$','_luLI', r) for r in micro_mets.loc[:,'metID']]
    micro_mets['metName'] = [re.sub(r'_d$','_luLI', r) for r in micro_mets.loc[:,'metName']]
    micro_mets['metComp'] = [re.sub(r'd$','luLI', r) for r in micro_mets.loc[:,'metComp']]
    micro_mets['metID'] = [re.sub(r'_u$','_luM', r) for r in micro_mets.loc[:,'metID']]
    micro_mets['metName'] = [re.sub(r'_u$','_luM', r) for r in micro_mets.loc[:,'metName']]
    micro_mets['metComp'] = [re.sub(r'u$','luM', r) for r in micro_mets.loc[:,'metComp']]
    
    # Remove fecal transport reactions
    micro_rxns = micro_rxns.loc[~micro_rxns['rxnID'].str.startswith('UFEt_'), :]
    
    # Rename microbiota biomass
    micro_mets.loc[micro_mets['metID'] == 'microbeBiomass_luM', 'metID'] = 'microbiota_LI_biomass_luM'
    micro_mets.loc[micro_mets['metName'] == 'microbeBiomass_luM', 'metName'] = 'microbiota_LI_biomass_luM'
    micro_rxns.loc[:,'rxnEq'] = [re.sub(r'microbeBiomass_luM','microbiota_LI_biomass_luM', r) for r in micro_rxns.loc[:,'rxnEq']]

    # Add biomass transport and excretion reactions
    data = [['LI_EX_microbiota_LI_biomass_luLI_fe','LI_EX_microbiota_LI_biomass_luLI_fe','Transport, biofluid', 0, 10000, 'microbiota_LI_biomass_luM;microbiota_LI_biomass_fe', '-1;1'],
            ['Excretion_EX_microbiota_LI_biomass_fe', 'Excretion_EX_microbiota_LI_biomass_fe', 'Exchange/demand reaction', 0, 10000, 'microbiota_LI_biomass_fe', '-1']]
    micro_rxns = pd.concat([micro_rxns, pd.DataFrame(data, columns = ['rxnID','rxnName','rxnSubs','lb','ub','rxnEq','rxnS'])]).reset_index(drop = True)
    
    # Rename fecal microbiota biomass
    micro_mets.loc[micro_mets['metID'] == 'microbeBiomass_fe', 'metID'] = 'microbiota_LI_biomass_fe'
    micro_mets.loc[micro_mets['metName'] == 'microbeBiomass_fe', 'metName'] = 'microbiota_LI_biomass_fe'
    
    # Adjust bounds
    micro_rxns.loc[micro_rxns['lb'] < 0, 'lb'] = -1000 * 1000
    micro_rxns.loc[micro_rxns['ub'] > 0, 'ub'] = 1000 * 1000
    
    # Remove overlapping reactions between WBM and microbiota model
    if any(micro_rxns['rxnID'].isin(wbm_rxns['rxnID'])):
        micro_rxns = micro_rxns.loc[~micro_rxns['rxnID'].isin(wbm_rxns['rxnID']),:]
        
    # Adjust communityBiomass to percentage
    bIdx = micro_rxns['rxnID'] == 'communityBiomass'
    micro_rxns.loc[bIdx,'rxnS'] = ';'.join([str(float(x) * 100) for x in micro_rxns.loc[bIdx, 'rxnS'].str.split(';').to_list()[0]])
    
    # Make sure that all new metabolites in luLI compartment can be excreted
      
    # Concatenate whole-body and microbiota models information
    mInfo = pd.concat([wbm_mets, micro_mets]).reset_index(drop = True)
    rInfo = pd.concat([wbm_rxns, micro_rxns]). reset_index(drop = True)
    
    # Build model
    model = _addInfoToModel(name, wbm_model, micro_mets, micro_rxns)
    
    # Detect metabolites involved in 0 reactions
    model_mets = [m.id for m in model.metabolites]
    mets_in_rxns = set()
    for rxn in model.reactions:
        mets_in_rxns.update([m.id for m in rxn.metabolites])
    orphan_mets = set(model_mets) - mets_in_rxns
   
    # Remove orphan metabolites
    mInfo = mInfo.loc[~mInfo['metID'].isin(orphan_mets),:]
    for met_id in orphan_mets:
        met = model.metabolites.get_by_id(met_id)
        model.remove_metabolites([met])
    
    # Add microbiota coupling constraints
    bioMets = rInfo.loc[rInfo['rxnID']=='communityBiomass', 'rxnEq'].str.split(';').to_list()[0]
    taxa = [re.sub(r'_biomass_c$','',x) for x in bioMets if x != 'microbiota_LI_biomass_luM']
    for taxon in taxa:
        model = _addCouplingConstraints(model, taxon, 400.0, 0.0)
    
    # Set objective
    model.objective = 'Whole_body_objective_rxn'
    
    # Additional constraints: Limit WBM and microbiota flux
    rxnWBMBio = model.reactions.get_by_id('Whole_body_objective_rxn')
    rxnWBMBio.lower_bound = 1.0
    rxnWBMBio.upper_bound = 1.0
    rxnMicroBio = model.reactions.get_by_id('Excretion_EX_microbiota_LI_biomass_fe')
    rxnMicroBio.lower_bound = 1.0
    rxnMicroBio.upper_bound = 1.0
    
    # Output variable
    return(model, mInfo, rInfo)

def _runGenerationMWBM(sample: str, meta: pd.DataFrame, use_pers_wbm: bool, pers_wbm_path: str, micro_mod_path: str, diet_file_name: str, output_dir: str, solver: str) -> (pd.DataFrame):
    if (not os.path.exists(os.path.join(output_dir, 'models', sample + '_mWBM.pkl'))):
        # Load WBM information
        if use_pers_wbm:
            try:
                wbm_metInfo = pd.read_csv(os.path.join(pers_wbm_path, sample + '_WBM.metInfo.csv'))
                wbm_rxnInfo = pd.read_csv(os.path.join(pers_wbm_path, sample + '_WBM.rxnInfo.csv'))
                with open(os.path.join(pers_wbm_path, sample + '.pkl'), 'rb') as f:
                    wbm_model = pickle.load(f)
            except:
                raise ValueError(f'Personalized WBM of sample {sample} can not be found at {pers_wbm_path}. Execute personalizeWBM function first.')
        else:
            dataPath = pkg_resources.files('persephone').joinpath('data/models/WBM')
            if meta.loc[sample]['sex'].lower() == 'male':
                wbm_metInfo = pd.read_csv(dataPath.joinpath('Harvey_1_03d' + '.metInfo.csv'))
                wbm_rxnInfo = pd.read_csv(dataPath.joinpath('Harvey_1_03d' + '.rxnInfo.csv'))
                with open(str(dataPath.joinpath('Harvey_1_03d.pkl')), 'rb') as f:
                    wbm_model = pickle.load(f)
            elif meta.loc[sample]['sex'].lower() == 'female':
                wbm_metInfo = pd.read_csv(dataPath.joinpath('Harvetta_1_03d' + '.metInfo.csv'))
                wbm_rxnInfo = pd.read_csv(dataPath.joinpath('Harvetta_1_03d' + '.rxnInfo.csv'))
                with open(str(dataPath.joinpath('Harvetta_1_03d.pkl')), 'rb') as f:
                    wbm_model = pickle.load(f)
            else: 
                raise ValueError('Invalid sex type. Options are [male, female].')
        
        # Load sample microbiota model
        try:
            micro_metInfo = pd.read_csv(os.path.join(micro_mod_path, sample + '.metInfo.csv'))
            micro_rxnInfo = pd.read_csv(os.path.join(micro_mod_path, sample + '.rxnInfo.csv'))
        except:
            raise ValueError(f'Microbiota model of sample {sample} can not be found at {micro_mod_path}. Execute runMgPipe function first.')
        
        # Define final model name
        if use_pers_wbm:
            name = 'miWBM_' + sample
        else:
            name = 'mWBM_' + sample
        
        # Merge WBM and microbiota models
        model, mInfo, rInfo = _combineHarveyMicrobiota(name, wbm_model, wbm_metInfo, wbm_rxnInfo, micro_metInfo, micro_rxnInfo)
        
        # Set diet constraints
        diet = pd.read_table(pkg_resources.files('persephone').joinpath('data/diet/' + diet_file_name + 'WBM.txt'), sep = '\t')
        model = _setDietConstraints(model, diet)
            
        # Export model information
        mInfo.to_csv(os.path.join(output_dir, 'models', sample + '_mWBM.metInfo.csv'), index = False)
        rInfo.to_csv(os.path.join(output_dir, 'models', sample + '_mWBM.rxnInfo.csv'), index = False)
        with open(os.path.join(output_dir, 'models', sample + '_mWBM.pkl'), "wb") as f:
            pickle.dump(model, f)
    else:
        with open(os.path.join(output_dir, 'models', sample + '_mWBM.pkl'), "rb") as f:
            model = pickle.load(f)
           
    # Extract model statistics
    stats = {'Sample': sample, 'Metabolites': len(model.metabolites), 'Reactions': len(model.reactions), 'Constraints': len(model.constraints) - len(model.metabolites)}       
            
    # Output variable
    return (stats)

def generateMWBM(meta: pd.DataFrame, micro_mod_path: str = os.path.join('.','output','MgPipe','models'), diet_file_name: str = 'EUAverageDiet', use_pers_wbm: bool = False,
                pers_wbm_path: str = os.path.join('.','output','persWBM'), output_dir: str = os.path.join('.','output','MWBM'), solver: str = 'glpk',
                threads: int = round(os.cpu_count())) -> (pd.DataFrame):
    # Set COBRApy default solver
    try:
        cobra.Configuration().solver = solver
    except:
        print(f'COBRApy could not found installation for {solver}. Simulations will be done with default solver.')
        
    # Generate output folders
    os.makedirs(os.path.join(output_dir, 'models'), exist_ok = True)
    
    # Generate original WBM if required
    if not use_pers_wbm:
        # Generate Harvey model
        if (not os.path.exists(os.fspath(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvey_1_03d.pkl')))):
            print('Generating Harvey_1_03d original model. This could take a while...')
            male = _generateOriginalWBM('Harvey_1_03d', solver, threads)
            with open(str(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvey_1_03d.pkl')), "wb") as f:
                pickle.dump(male, f)
        
        # Generate Harvetta model
        if (not os.path.exists(os.fspath(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvetta_1_03d.pkl')))):
            print('Generating Harvetta_1_03d original model. This could take a while...')
            female = _generateOriginalWBM('Harvetta_1_03d', solver, threads)
            with open(str(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvetta_1_03d.pkl')), "wb") as f:
                pickle.dump(female, f)
    
    # Initialize output
    statsDict = {}
    
    # Create microbiota whole body models
    parModelsFunc = partial(_runGenerationMWBM, meta = meta, use_pers_wbm = use_pers_wbm, pers_wbm_path = pers_wbm_path, micro_mod_path = micro_mod_path, diet_file_name = diet_file_name, output_dir = output_dir, solver = solver)
    with ProcessPoolExecutor(max_workers = threads) as executor:
        # Run function
        futures = {executor.submit(parModelsFunc, sample): sample for sample in meta.index}
        
        # Progress bar
        for future in tqdm(as_completed(futures), total=len(futures), desc="Generating microbiota whole body models"):
            sample = futures[future]
            statsDict[sample] = future.result()
            
    # Create dataframe
    stats = pd.DataFrame.from_dict(statsDict, orient = "index")
    stats.index = stats.loc[:,'Sample']
    stats = stats.select_dtypes(include = ['number'])
    
    # Run model optimization with all available threads
    for sample in tqdm(meta.index, desc = "Optimizing microbiota whole body models"):
        # Load model
        with open(os.path.join(output_dir, 'models', sample + '_mWBM.pkl'), "rb") as f:
            model = pickle.load(f)
            
        # Setting solver parameters
        #model = _setSolverParameters(model, solver, threads)
        print(model.solver.problem.parameters.threads.get())
        
        # Optimize model
        model.optimize()
        
        # Save model
        with open(os.path.join(output_dir, 'models', sample + '_mWBM.pkl'), "wb") as f:
            pickle.dump(model, f)
    
    # Output variable
    return (stats)

def _prepareModel(model: cobra.Model, rxnList: list) -> (cobra.Model):
    # Extract model information
    modMets = [x.id for x in model.metabolites]
    modRxns = [x.id for x in model.reactions]
    
    # Find demand reactions not yet in the model
    dmReactions = [x for x in rxnList if x not in modRxns]
    
    # Add demand reactions for missing metabolites
    if any([x.startswith('DM_') for x in dmReactions]):
        # Extract demand reaction
        demandRxns = [x for x in dmReactions if x.startswith('DM_')]
        demandMets = [re.sub('DM_','',x) for x in demandRxns]
        n = len(demandRxns)
        
        # Check missing metabolites
        missingMets = [x for x in demandMets if x not in modMets]
        
        # Prepare metabolite data TODOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
        metInfo = pd.DataFrame()
        
        # Prepare reaction data
        rxnInfo = pd.DataFrame({'rxnID': rxnList, 'rxnName': rxnList, 'rxnSubs': ['Demand']*n, 'lb': [0]*n, 'ub': [100000]*n, 'rxnS': ['-1']*n, 'rxnEq': demandMets})
        
        # Add information to the model
        model = _addInfoToModel(model.id, model, metInfo, rxnInfo)
        
    else:
        print('The following reactions are not present in the model')
        print(dmReactions)
        
    # Output variable
    return (model, demandRxns)

def _runOptimizationWBM(sample: str, rxn_list_path: str, mwbm_path: str, output_dir: str, solver: str):
    # Load model
    with open(os.path.join(mwbm_path, sample + '_mWBM.pkl'), "rb") as f:
         model = pickle.load(f)
    
    # Set solver parameters
    model = _setSolverParameters(model, solver)
    
    # Load reaction list
    rxnList = pd.read_table(rxn_list_path).iloc[:,0].to_list()
    
    # Add missing demand reactions to the model
    model, rxnList = _prepareModel(model, rxnList)
    
    # Extract model metabolites and reactions
    modMets = [x.id for x in model.metabolites]
    modRxns = [x.id for x in model.reactions]
    
    # Define sample sex
    if any('Testis' in rxn for rxn in modRxns):
        sex = 'male'
    else:
        sex = 'female'
    
    # Identify taxon names and relative abundances
    comBio = model.reactions.get_by_id('communityBiomass')
    taxonNames = [re.sub(r'_biomass_c','',t) for t in [x.id for x in comBio.metabolites] if t != 'microbiota_LI_biomass_luM']
    relAbundances = [-x for x in list(comBio.metabolites.values())[:-1]]
    
    # Preallocate solution
    solution = {'ID': model.id, 'Sex': sex, 'taxonNames': taxonNames, 'relAbundances': relAbundances, 'rxnList': rxnList, 'modelMets': modMets,
                'modelRxns': modRxns, 'Stat': pd.DataFrame(data = np.nan, index = ['0'], columns = rxnList),
                'Obj': pd.DataFrame(np.nan, index = ['0'], columns = rxnList), 'Flux': pd.DataFrame(np.nan, index = modRxns, columns = rxnList),
                'RedCost': pd.DataFrame(np.nan, index = modRxns, columns = rxnList), 'Shadow': pd.DataFrame(np.nan, index = modMets, columns = rxnList),
                'ShadowBIO': pd.DataFrame(np.nan, index = taxonNames, columns = rxnList)}
    
    # Make simulations
    for rxn in rxnList:
        print(f'Investigate reaction {rxn} in sample {sample}')
        
        # Change the objective
        model.objective = rxn
        
        try:
            # Perform FBA
            sol = model.optimize()
            
            # Status
            solution['Stat'].loc[:,rxn] = 1
            
            # Objective
            solution['Obj'].loc[:,rxn] = sol.objective_value
            
            # Flux vector, shadow prices and reduced costs
            solution['Flux'].loc[:,rxn] = sol.fluxes
            solution['Shadow'].loc[:,rxn] = sol.shadow_prices
            solution['RedCost'].loc[:,rxn] = sol.reduced_costs
            
            # Microbial biomass shadow prices
            solution['ShadowBIO'].loc[:,rxn] = sol.shadow_prices.loc[[x + '_biomass_c' for x in taxonNames]].values
            
        except:
            # Save status
            solution['Stat'].loc[:,rxn] = 0
    
    # Export solutions
    with open(os.path.join(output_dir, 'FBA_sol_' + model.id + '.pkl'), 'wb') as f:
        pickle.dump(solution, f)  
    
    # Output variable
    solution['Obj'].index = [sample]
    return solution['Obj']

def _describeFluxes(fluxes: pd.DataFrame, cutoff_type: str, cutoff_val: float, roundingFactor: float, analyse_gf: bool) -> (dict):
    # Remove sex information
    sexInfo = fluxes.loc[:,['ID','Sex']]
    fluxes.drop(['ID','Sex'], axis = 1, inplace = True)
    
    # Create temporary variable for statistics
    fluxesForStats = copy.deepcopy(fluxes)
    
    # Remove germfree samples
    
    # Preallocate table
    varNames = ['Reaction','No results','Unique results','Mean','Variance','Standard deviation','Removed']
    fluxStats = pd.DataFrame(index = fluxesForStats.columns, columns = varNames)
    
    # Add reaction names
    fluxStats.loc[:,'Reaction'] = fluxesForStats.columns
    
    # Find the number of samples with microbial flux contribution (previously zeros)
    fluxStats.loc[:,'No results'] = fluxesForStats.isna().sum()
    
    # Find the number of unique flux results per reaction
    fluxStats.loc[:,'Unique results'] = fluxesForStats.nunique(dropna = True)
    
    # Calculate mean of fluxes
    fluxStats.loc[:,'Mean'] = fluxesForStats.mean(skipna = True)
    
    # Calculate reaction variance
    fluxStats.loc[:,'Variance'] = fluxesForStats.var(skipna = True)
    
    # Calculate the standard deviation
    fluxStats.loc[:,'Standard deviation'] = fluxesForStats.std(skipna = True)
    
    # Remove metabolites based on user defined threshold value
    if cutoff_type == 'fraction':
        fluxStats.loc[:,'Removed'] = (fluxStats.loc[:,'Unique results'] / len(fluxStats.index)) < cutoff_val
    elif cutoff_type == 'SD':
        fluxStats.loc[:,'Removed'] = fluxStats.loc[:,'Standard deviation'] < cutoff_val
    elif cutoff_type == 'count':
        fluxStats.loc[:,'Removed'] = fluxStats.loc[:,'Unique results'] < cutoff_val
        
    # Remove reactions according to 'Removed'
    fluxes_rm = copy.deepcopy(fluxes)
    fluxes_rm = fluxes_rm.loc[:,~fluxStats['Removed'].astype(bool).values]
    
    # Scale fluxes
    
    # Define output variable
    stats = {'Flux_summary_statistics': fluxStats, 'Fluxes': fluxes, 'Fluxes_removed_reactions': fluxes_rm}
    return(stats)

def _analyzeFBAsol(meta: pd.DataFrame, output_dir: str, num_rounding: float = 1e-06, cutoff_type: str = 'fraction',
                   cutoff_val: float = 0.1, analyse_gf: bool = False):
    # Extract the number of digits to round from numerical rounding
    roundingFactor = -int(math.floor(math.log10(abs(num_rounding))))
    
    # PART 1: Load FBA solutions
    print('Loading the FBA solutions...')
    
    # Find paths to FBA solutions
    solPaths = [x for x in os.listdir(output_dir) if '.pkl' in x]
    modelNames = [x.replace('FBA_sol_','').removesuffix('.pkl') for x in solPaths]
    
    # Get number of metabolites investigated
    with open(os.path.join(output_dir, solPaths[0]), 'rb') as f:
        sol = pickle.load(f)
    reactions = sol['rxnList']
    
    # Find duplicate metabolites and remove where needed
    reactions = pd.Series(reactions).drop_duplicates().to_list()
    
    # Preallocate tables for FBA results
    fluxes = pd.DataFrame(np.nan, index = modelNames, columns = reactions)
    metadata = pd.DataFrame(index = modelNames, columns = ['ID','Sex'])
    fbaStats = pd.DataFrame(np.nan, index = modelNames, columns = reactions)
    
    # Load results and produce tables for the fluxes
    for mod in modelNames:
        # Load FBA solution
        with open(os.path.join(output_dir, 'FBA_sol_' + mod + '.pkl'), 'rb') as f:
            sol = pickle.load(f)
        
        # Add solution to the metadata table
        metadata.loc[mod,'ID'] = sol['ID']
        metadata.loc[mod,'Sex'] = sol['Sex']
        
        # Set flux results to nan if stat is not equal to one
        sol['Obj'][sol['Stat'] != 1] = np.nan
        
        # Round the fluxes and store them in table
        fluxes.loc[mod, :] = np.round(sol['Obj'].values, roundingFactor)
        
        # Add FBA statistics to tables
        fbaStats.loc[mod,:] = sol['Stat'].values
    
    # PART 2: Scale fluxes and produce summary statistics
    print('Processing and analysing the flux results...')
    
    # Obtain summary statistics
    fluxes = pd.concat([metadata, fluxes], axis = 1)
    
    # Create statistics for fluxes and prune results
    stats = _describeFluxes(fluxes, roundingFactor, cutoff_type, cutoff_val, analyse_gf)
    
    # Get fluxes for further analysis
    fluxesPruned = stats['Fluxes_removed_reactions']
    
    # Are the models personalised with gut microbiota?
    microbiomePresent = False
    if any([x.startswith('mWBM_') for x in modelNames]) or any([x.startswith('miWBM_') for x in modelNames]):
        microbiomePresent = True
    
    # Correlate fluxes with microbial relative abundances
    if microbiomePresent:
        print('> Extract metagenomic relative abundances from mWBMs...')
        
        # Preallocate microbe info dictionary
        modelSP = dict.fromkeys(modelNames)
        
        # Load relative abundances and taxa from fba results
        for solution in solPaths:
            # Load FBA solution
            with open(os.path.join(output_dir, solution), 'rb') as f:
                sol = pickle.load(f)
                
            # Attach information
            modelSP[sol['ID']] = pd.DataFrame({'Taxon': sol['taxonNames'], 'Abundance': sol['relAbundances']})
            
        # Find all microbial species
        allSpecies = set()
        for df in modelSP.values():
            allSpecies.update(df['Taxon'].unique())
        allSpecies = list(sorted(allSpecies))
        
        # Find the microbial relative abundances from the models
        relAbun = pd.DataFrame(index = allSpecies, columns = modelNames)
        for model, df in modelSP.items():
            ab = df.set_index('Taxon')['Abundance'] / 100
            relAbun[model] = ab
        relAbun.index = [re.sub(r'^pan','',x) for x in relAbun.index]
        
        # Save microbiome relative abundances
        print('> Wrote metagenomic relative abundances from mWBMs to file...')
        relAbun.transpose().to_csv(os.path.join(output_dir,'FluxAnalysis','WBM_relative_abundances.csv'))
        
        # Process flux data from samples
        print('> Perform spearman correlations on flux results and relative abundances...')
        fluxesToCorrelate = stats['Fluxes']
        
        fluxesToCorrelate.index = fluxesToCorrelate.index.str.replace('|'.join(['^mWBM_','^miWBM_','^iWBM_']), '', regex = True)
    
def analyzeMWBM(meta: pd.DataFrame, rxn_list_path: str, mwbm_path: str = os.path.join('.','output','mWBM','models'), output_dir: str = os.path.join('.','output','mWBM','resultFlux'),
                solver: str = 'glpk', threads: int = round(os.cpu_count())) -> (pd.DataFrame):
    # Set COBRApy default solver
    try:
        cobra.Configuration().solver = solver
    except:
        print(f'COBRApy could not found installation for {solver}. Simulations will be done with default solver.')
        
    # Generate output folders
    os.makedirs(os.path.join(output_dir, 'FluxAnalysis'), exist_ok = True)
    
    # Initialize results dataframe
    res  = []
    
    # Perform FBA to input reactions
    parOptFunc = partial(_runOptimizationWBM, rxn_list_path = rxn_list_path, mwbm_path = mwbm_path, output_dir = output_dir, solver = solver)
    with ProcessPoolExecutor(max_workers = threads) as executor:
        # Run function
        futures = {executor.submit(parOptFunc, sample): sample for sample in meta.index}
        
        for future in tqdm(as_completed(futures), total=len(futures), desc="Performing FBA on input reaction list"):
            res.append(future.result())
    
    # Process FBA solutions
    
    # Output variable
    res = pd.concat(res, axis = 0).T
    return(res)