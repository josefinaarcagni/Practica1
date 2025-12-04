# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 13:17:46 2025

@author: tblasco
"""

import os
import re
import pickle
import pandas as pd
import numpy as np
import importlib.resources as pkg_resources
from functools import partial
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
import cobra
from cobra.flux_analysis import flux_variability_analysis
from .utils import _labelTaxonModel, _buildModel, _addCouplingConstraints, _adaptDiet, _useDiet

def _calculateRxnAbundance(abundMat: pd.DataFrame, taxInfo: pd.DataFrame, rxnInfo: pd.DataFrame, rxnTaxMat: pd.DataFrame, tol: float = 0.0000001) -> (pd.DataFrame, pd.DataFrame):
    # Extract taxa
    taxa = taxInfo['Taxon']

    # Remove column from abundance matrix
    ab = abundMat.set_index(taxa)
        
    # Select tax
    rxnTaxMat = rxnTaxMat[ab.index]
    
    # Calculate reaction abundance
    rxnAbundance = rxnTaxMat @ ab
    rxnAbundance = rxnAbundance.apply(pd.to_numeric, errors = 'coerce')
    
    # Add reaction IDs
    rxnAbundance.index = rxnInfo['rxnID']
    
    # Remove reactions not present in the models
    rxnAbundance = rxnAbundance.loc[~(rxnAbundance == 0).all(axis=1)]
    
    # Calculate reaction presence
    rxnPresence = rxnAbundance
    rxnPresence = pd.DataFrame(np.where(rxnPresence > tol, 1, 0))
    rxnPresence.columns = rxnAbundance.columns
        
    # Add reaction IDs
    rxnPresence.index = rxnAbundance.index
    
    # Output variable
    return (rxnAbundance, rxnPresence)

def _calculateSubsAbundance(rxnAb: pd.DataFrame, rxnInfo: pd.DataFrame) -> (pd.DataFrame):
    # Filter reactions present in the models
    rxnInfo = rxnInfo[rxnInfo['rxnID'].isin(rxnAb.index)]
    
    # Define subsystems
    subs = sorted(rxnInfo['rxnSubs'].dropna().unique().tolist())
    
    # Remove biomass reaction
    rxnInfo = rxnInfo[~rxnInfo['rxnID'].str.contains('biomass', case=False, na=False)]
    
    # Define subsystems-reaction ratio
    subsRxnMat = pd.crosstab(rxnInfo['rxnSubs'], rxnInfo['rxnID'])
    subsRxnMat = subsRxnMat.div(subsRxnMat.sum(axis = 1), axis = 0).fillna(0)
    
    # Filter reaction abundance scores
    rxnAb = rxnAb.loc[subsRxnMat.columns]
    
    # Calculate subsystem abundance scores
    subsAbundance = subsRxnMat @ rxnAb
    subsAbundance = subsAbundance.apply(pd.to_numeric, errors = 'coerce')
    subsAbundance.index = subs
    
    # Output variable
    return (subsAbundance)

def _defineActiveExchangesMat(taxa: pd.DataFrame, database: str) -> (pd.DataFrame):
    # Load database information
    dataPath = pkg_resources.files('persephone').joinpath('data/models/' + database)
    metInfo = pd.read_csv(dataPath.joinpath(database + '.metInfo.csv'))
    rxnInfo = pd.read_csv(dataPath.joinpath(database + '.rxnInfo.csv'))
    metTaxMat = pd.read_csv(dataPath.joinpath(database + '.metTaxMat.csv'))
    rxnTaxMat = pd.read_csv(dataPath.joinpath(database + '.rxnTaxMat.csv'))
    
    # Inititalize output
    metInfoDataSet = sorted(metInfo.loc[metTaxMat.loc[:,taxa.iloc[:,0]].sum(axis = 1) > 0, 'metID'].tolist())
    activeExMat = pd.DataFrame(data = 0, index = [re.sub(r'_e$', '(e)', m) for m in ['EX_' + y for y in [x for x in metInfoDataSet if x.endswith('_e')]]], columns = taxa.iloc[:,0]).sort_index()
    
    # Check exchanges feasibility
    for taxon in taxa.iloc[:,0]:
        taxa_metInfo = metInfo.loc[metTaxMat.loc[:,taxon] == 1].reset_index()
        taxa_rxnInfo = rxnInfo.loc[rxnTaxMat.loc[:,taxon] == 1].reset_index()
        
        # Build model
        m = _buildModel(taxon, taxa_metInfo, taxa_rxnInfo)
        m.objective = 'biomassPan_' + taxon
        
        # Find exchange reactions in taxa
        exRxns = [x for x in taxa_rxnInfo['rxnID'] if x.startswith('EX_') and x.endswith('(e)') and 'biomass' not in x]
        
        # Apply flux variability analysis
        flux = flux_variability_analysis(m, exRxns, fraction_of_optimum = 0)
        
        # Add information to the matrix
        activeExMat.loc[flux[(flux["minimum"].abs() > 0.00000001) | (flux["maximum"].abs() > 0.00000001)].index, taxon] = 1
    
    # Return output
    return (activeExMat)

def _adaptTaxonTransports(taxon: str, mInfo: pd.DataFrame, rInfo: pd.DataFrame, activeExMat: pd.DataFrame) -> (pd.DataFrame):
    # Remove exchange reactions
    rInfo = rInfo[~rInfo.loc[:,'rxnID'].str.contains('_EX_')]

    # Define lumen reactions and metabolites
    lumenRxns = activeExMat.index.tolist()
    lumenMets = [re.sub(r'^EX_','',r) for r in lumenRxns]
    lumenMets = [re.sub(r'\(e\)$','_u',m) for m in lumenMets]
    
    # Identify lumen metabolites in the taxon
    taxMets = [re.sub(taxon + '_', '', x) for x in mInfo.loc[:,'metID']]
    lumenTaxMets = [x for x in lumenMets if x in taxMets]
    
    # Add transport reaction data between lumen and taxon
    rtList = []
    for met in lumenTaxMets:
        data = [[taxon + '_IEX_' + met, taxon + '_IEX_' + met, '', -1000, 1000, met + ';' + (taxon + '_' + met), '-1;1']]
        rtList.append(pd.DataFrame(data, columns = ['rxnID','rxnName','rxnSubs','lb','ub','rxnEq','rxnS']))
        
    lumenTrRxnsInfo = pd.concat(rtList).reset_index(drop = True)
    
    # Concatenate information
    rInfo = pd.concat([lumenTrRxnsInfo, rInfo]).reset_index(drop = True)
    
    # Output variable
    return (rInfo)

def _addDietMetabolites(taxa: list, mInfo: pd.DataFrame, rInfo: pd.DataFrame, activeExMat: pd.DataFrame, origMetInfo: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame):
    # Define active lumen reactions and metabolites
    lumenRxns = activeExMat.index.tolist()
    lumenMets = [re.sub(r'^EX_','',r) for r in lumenRxns]
    lumenMets = [re.sub(r'\(e\)$','_u',m) for m in lumenMets]
    
    # Define active lumen reactions and metabolites
    lumenActRxns = activeExMat.index[activeExMat.sum(axis = 1) > 0].tolist()
    lumenActMets = [re.sub(r'^EX_','',r) for r in lumenActRxns]
    lumenActMets = [re.sub(r'\(e\)$','_u',m) for m in lumenActMets]
    
    # Define lumen metabolites data
    origMetInfo['metID'] = origMetInfo['metID'].str.replace(r'_e$', '_u', regex=True)
    
    mList = []
    for met in lumenMets:
        if met in lumenActMets:
            data = [[re.sub(r'_u$', '_d', met), origMetInfo.loc[origMetInfo['metID']==met, 'metName'].iloc[0], origMetInfo.loc[origMetInfo['metID']==met, 'metFormula'].iloc[0], 'd'],
                    [re.sub(r'_u$', '_fe', met), origMetInfo.loc[origMetInfo['metID']==met, 'metName'].iloc[0], origMetInfo.loc[origMetInfo['metID']==met, 'metFormula'].iloc[0], 'fe'],
                    [met, origMetInfo.loc[origMetInfo['metID']==met, 'metName'].iloc[0], origMetInfo.loc[origMetInfo['metID']==met, 'metFormula'].iloc[0], 'u']]
        else:
            data = [[met, origMetInfo.loc[origMetInfo['metID']==met, 'metName'].iloc[0], origMetInfo.loc[origMetInfo['metID']==met, 'metFormula'].iloc[0], 'u']]
        
        mList.append(pd.DataFrame(data, columns = ['metID','metName','metFormula','metComp']))
    
    lumenMetsInfo = pd.concat(mList).reset_index(drop = True)
    
    # Define lumen reactions data
    rList = []
    for met in lumenActMets:
        data = [['Diet_EX_' + re.sub(r'_u$', '_d', met), 'EX_' + re.sub(r'_u$', '_d', met), '', -1000, 10000, re.sub(r'_u$', '_d', met), '-1'],
                ['DUt_' + re.sub(r'_u$', '', met), 'DUt_' + re.sub(r'_u$', '', met), '', 0, 10000, ';'.join([re.sub(r'_u$', '_d', met), met]), '-1;1'],
                ['UFEt_' + re.sub(r'_u$', '', met), 'UFEt_' + re.sub(r'_u$', '', met), '', 0, 10000, ';'.join([met, re.sub(r'_u$', '_fe', met)]), '-1;1'],
                ['EX_' + re.sub(r'_u$', '_fe', met), 'EX_' + re.sub(r'_u$', '_fe', met), '', -1000, 10000, re.sub(r'_u$', '_fe', met), '-1']]
        rList.append(pd.DataFrame(data, columns = ['rxnID','rxnName','rxnSubs','lb','ub','rxnEq','rxnS']))
    
    lumenRxnsInfo = pd.concat(rList).reset_index(drop = True)

    # Merge metabolites and reactions information
    mInfo = pd.concat([lumenMetsInfo, mInfo]).reset_index(drop = True)
    rInfo = pd.concat([lumenRxnsInfo, rInfo]).reset_index(drop = True)
    
    # Remove any metabolite not present in the reactions
    all_rxn_mets = ';'.join(rInfo['rxnEq'].astype(str).tolist())
    mInfo = mInfo[mInfo['metID'].apply(lambda x: x in all_rxn_mets)].copy()
    
    # Output variable
    return mInfo, rInfo

def _createPersonalizedModels(name: str, taxInfo: pd.DataFrame, abundMat: pd.DataFrame, database: str, activeExMat: pd.DataFrame) -> (cobra.Model, pd.DataFrame, pd.DataFrame):
    # Load database information
    dataPath = pkg_resources.files('persephone').joinpath('data/models/' + database)
    metInfo = pd.read_csv(dataPath.joinpath(database + '.metInfo.csv'))
    rxnInfo = pd.read_csv(dataPath.joinpath(database + '.rxnInfo.csv'))
    metTaxMat = pd.read_csv(dataPath.joinpath(database + '.metTaxMat.csv'))
    rxnTaxMat = pd.read_csv(dataPath.joinpath(database + '.rxnTaxMat.csv'))
    
    # Identify sample species and abundances
    idx = abundMat[name] > 0
    taxa = taxInfo['Taxon'][idx].tolist()
    ab = abundMat[name][idx].tolist()
    
    # Label reactions and metabolites
    mInfo = []
    rInfo = []
    for taxon in taxa:
        singMets, singRxns = _labelTaxonModel(taxon, metInfo, rxnInfo, metTaxMat, rxnTaxMat)
        mInfo.append(singMets)
        rInfo.append(singRxns)
        
    # Add transport reactions between the lumen and taxon
    for i in range(len(taxa)):
        rInfo[i] = _adaptTaxonTransports(taxa[i], mInfo[i], rInfo[i], activeExMat)
    
    mInfo = pd.concat(mInfo).reset_index(drop = True)
    rInfo = pd.concat(rInfo).reset_index(drop = True)
    
    # Add diet and fecal compartments
    mInfo, rInfo = _addDietMetabolites(taxa, mInfo, rInfo, activeExMat, metInfo)
    
    # Define biomass metabolites
    bioMets = [taxon + '_biomass_c' for taxon in taxa]
    
    # Add community biomass metabolites information
    data = [['microbeBiomass_u', 'microbeBiomass_u', '', 'u'],['microbeBiomass_fe', 'microbeBiomass_fe', '', 'fe']]
    mInfo = pd.concat([mInfo, pd.DataFrame(data, columns = ['metID', 'metName', 'metFormula', 'metComp'])]).reset_index(drop = True)
    
    # Add community biomass reactions information
    data = [['communityBiomass', 'communityBiomass', '', 0, 10000, ';'.join(bioMets + ['microbeBiomass_u']), ';'.join(map(str, [-x for x in ab]+[1]))],
            ['UFEt_microbeBiomass','UFEt_microbeBiomass','', 0, 10000, 'microbeBiomass_u;microbeBiomass_fe', '-1;1'],
            ['EX_microbeBiomass_fe', 'EX_microbeBiomass_fe', '', -10000, 10000, 'microbeBiomass_fe', '-1']]
    rInfo = pd.concat([rInfo, pd.DataFrame(data, columns = ['rxnID','rxnName','rxnSubs','lb','ub','rxnEq','rxnS'])]).reset_index(drop = True)
    
    # Build model
    model = _buildModel(name, mInfo, rInfo)
    
    # Define model objective
    model.objective = 'EX_microbeBiomass_fe'
    
    # Add coupling constraints
    print('Adding constraints to sample ' + name)
    for taxon in taxa:
        model = _addCouplingConstraints(model, taxon, 400.0, 0.0)
        
    # Output variable
    return (model, mInfo, rInfo)

def _runCreatePersonalizedModels(sample: str, taxa: pd.DataFrame, counts: pd.DataFrame, database: str, activeExMat: pd.DataFrame, output_dir: str):
    # Check if model exists
    if (not os.path.exists(os.path.join(output_dir, 'models', sample + '.pkl'))):
        # Create model
        model, mInfo, rInfo = _createPersonalizedModels(sample, taxa, counts, database, activeExMat)
        
        # Optimize model to avoid issues
        model.optimize()
        
        # Export model information
        mInfo.to_csv(os.path.join(output_dir, 'models', sample + '.metInfo.csv'), index = False)
        rInfo.to_csv(os.path.join(output_dir, 'models', sample + '.rxnInfo.csv'), index = False)
        with open(os.path.join(output_dir, 'models', sample + '.pkl'), "wb") as f:
            pickle.dump(model, f)
            
    # Return
    return f'Generated model for sample {sample}.'

def _computeProfiles(model: cobra.Model, rxnInfo: pd.DataFrame, diet: pd.DataFrame, objLB: float = 0.4, objUB: float = 1.0) -> (pd.DataFrame):
    # Adapt diet
    dietAdapted = _adaptDiet(diet)
    
    # Apply diet to the model
    modelM = _useDiet(model, dietAdapted)
    
    # Additional constraints: 1) Manage demand and sink reaction bounds
    rxnBio = modelM.reactions.get_by_id('communityBiomass')
    components = [re.sub(r'_biomass_c$','', comp) for comp in [x.id for x, y in rxnBio.metabolites.items()]]
    for sp in components:
        #  Demand reactions
        dmRxns = [x for x in rxnInfo.loc[:,'rxnID'] if x.startswith(sp +'_DM_')]
        for r in dmRxns:
            rxn = modelM.reactions.get_by_id(r)
            rxn.lower_bound = 0.0
        
        # Sink reactions
        sinkRxns = [x for x in rxnInfo.loc[:,'rxnID'] if x.startswith(sp +'_sink_')]
        for r in sinkRxns:
            rxn = modelM.reactions.get_by_id(r)
            rxn.lower_bound = -1.0
            
    # Additional contraints: 2) Microbe biomass lb = 0.4, ub = 1
    rxnBio.lower_bound = objLB
    rxnBio.upper_bound = objUB

    # Additional constraints: 3) UFEt, DUt, EX reactions upper bound
    rxns = [x for x in rxnInfo.loc[:,'rxnID'] if x.startswith('UFEt_') or x.startswith('DUt_') or x.startswith('EX_')]
    for r in rxns:
        rxn = modelM.reactions.get_by_id(r)
        rxn.upper_bound = 1000000.0
    
    # Identify secretion reactions
    objRxnInfo = rxnInfo[rxnInfo['rxnID'].str.endswith('_fe') & ~rxnInfo['rxnID'].str.startswith('EX_microbeBiomass')]
    secRxns = [model.reactions.get_by_id(r) for r in objRxnInfo['rxnID']]
    
    # Apply flux variability analysis
    fluxesSec = flux_variability_analysis(modelM, secRxns, fraction_of_optimum = 0.9999)
    
    # Identify uptake reactions
    objRxnInfo = rxnInfo[rxnInfo['rxnID'].str.endswith('_d')]
    upRxns = [model.reactions.get_by_id(r) for r in objRxnInfo['rxnID']]
    
    # Apply flux variability analysis
    fluxesUp = flux_variability_analysis(modelM, upRxns, fraction_of_optimum = 0.9999)
    
    # Calculate final fluxes
    fluxes = pd.DataFrame({'minimum': abs(fluxesUp['maximum'].to_numpy() + fluxesSec['minimum'].to_numpy()),
                           'maximum': abs(fluxesSec['maximum'].to_numpy() + fluxesUp['minimum'].to_numpy())},
                          index = fluxesSec.index)
    
    # Output variable
    return (fluxes)

def _runComputeProfiles(sample: str, diet_file_name: str, output_dir: str) -> (pd.DataFrame):
    # Read diet file
    diet = pd.read_table(pkg_resources.files('persephone').joinpath('data/diet/' + diet_file_name + '.txt'), sep = '\t')
    
    # Read model information
    rxnInfo = pd.read_csv(os.path.join(output_dir, 'models', sample + '.rxnInfo.csv'))
    with open(os.path.join(output_dir, 'models', sample + '.pkl'), 'rb') as f:
        model = pickle.load(f)
        
    # Optimize model to avoid issues
    model.optimize()
    
    # Calculate secretion and uptake fluxes
    flux = _computeProfiles(model, rxnInfo, diet)
    
    # Output variable
    return flux
    
def runMgPipe(counts: pd.DataFrame, taxa: pd.DataFrame, database: str = 'AGORA2-APOLLO', compute_profiles: bool = False, solver: str = 'glpk',
              threads: int = round(os.cpu_count() * 0.8), diet_file_name: str = 'EUAverageDiet', output_dir: str = os.path.join('.','output','MgPipe')) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame):
    # Set COBRApy default solver
    try:
        cobra.Configuration().solver = solver
    except:
        print(f'COBRApy could not found installation for {solver}. Simulations will be done with default solver.')    
    
    # Generate output folders
    os.makedirs(os.path.join(output_dir, 'models'), exist_ok = True)
    
    # Load database information
    dataPath = pkg_resources.files('persephone').joinpath('data/models/' + database)
    rxnInfo = pd.read_csv(dataPath.joinpath(database + '.rxnInfo.csv'))
    rxnTaxMat = pd.read_csv(dataPath.joinpath(database + '.rxnTaxMat.csv'))
    
    # Calculate reaction scores
    rxnAbundance, rxnPresence = _calculateRxnAbundance(counts, taxa, rxnInfo, rxnTaxMat)

    # Calculate subsystem scores
    subsAbundance = _calculateSubsAbundance(rxnAbundance, rxnInfo)
    
    # Export abundances to text files
    rxnAbundance.to_csv(os.path.join(output_dir, 'rxnAbundance.csv'), index = True)
    rxnPresence.to_csv(os.path.join(output_dir, 'rxnPresence.csv'), index = True)
    subsAbundance.to_csv(os.path.join(output_dir, 'subsAbundance.csv'), index = True)
    
    # Define active exchanges matrix
    if os.path.exists(os.path.join(output_dir,'activeExMets.pkl')):
        with open(os.path.join(output_dir,'activeExMets.pkl'), 'rb') as f:
            activeExMat = pickle.load(f)
    else:
        activeExMat = _defineActiveExchangesMat(taxa, database)
        with open(os.path.join(output_dir,'activeExMets.pkl'), 'wb') as f:
            pickle.dump(activeExMat, f)
    
    # Create personalized microbiota models
    parModelsFunc = partial(_runCreatePersonalizedModels, taxa = taxa, counts = counts, database = database, activeExMat = activeExMat, output_dir = output_dir)
    with ProcessPoolExecutor(max_workers = threads) as executor:
        # Run function
        futures = [executor.submit(parModelsFunc, sample) for sample in counts.columns]
        
        # Progress bar
        for future in tqdm(as_completed(futures), total=len(futures), desc="Generating microbiota models"):
            result = future.result()
            print(result)
        
    # Calculate uptake and secretion fluxes
    if compute_profiles:
        # Initialize output
        fluxesFVA = {}
        
        # Compute profiles
        parProfilesFunc = partial(_runComputeProfiles, diet_file_name = diet_file_name, output_dir = output_dir)
        with ProcessPoolExecutor(max_workers = threads) as executor:
            # Run function
            futures = {executor.submit(parProfilesFunc, sample): sample for sample in counts.columns}
            
            # Progress bar
            for future in tqdm(as_completed(futures), total=len(futures), desc="Computing flux profiles"):
                sample = futures[future]
                fluxesFVA[sample] = future.result()
                
        # Build net fluxes dataframes
        netUptakeFluxes = pd.concat([fluxesFVA[sample]['minimum'].rename(sample) for sample in counts.columns], axis = 1).sort_index()
        netSecretionFluxes = pd.concat([fluxesFVA[sample]['maximum'].rename(sample) for sample in counts.columns], axis = 1).sort_index()
        
        # Export results
        netUptakeFluxes.to_csv(os.path.join(output_dir, 'netUptakeFluxes.csv'), index = True)
        netSecretionFluxes.to_csv(os.path.join(output_dir, 'netSecretionFluxes.csv'), index = True)
          
    # Output variable
    return (rxnAbundance, rxnPresence, subsAbundance)