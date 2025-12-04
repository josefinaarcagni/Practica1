# -*- coding: utf-8 -*-
"""
Created on Mon Aug 11 13:11:08 2025

@author: tblasco
"""

import re
import pandas as pd
import cobra
import copy
from cobra import Model, Metabolite, Reaction

def _setSolverParameters(model: cobra.Model, solverName: str, threads: int = 0) -> (cobra.Model):
    if solverName == 'cplex':
        # No scaling problem matrix
        model.solver.problem.parameters.read.scale.set(-1)
        
        # Emphasize precision in numerically unstable or difficult problems
        model.solver.problem.parameters.emphasis.numerical.set(1)
        
        # Time limit in seconds
        model.solver.problem.parameters.timelimit.set(100000)
        
        # Threads
        if threads > 0:
            model.solver.problem.parameters.threads.set(threads)
        
    # Output variable
    return (model)

def _buildModel(name: str, metInfo: pd.DataFrame, rxnInfo: pd.DataFrame):
    # Number of metabolites and reactions
    nM = metInfo.shape[0]
    nR = rxnInfo.shape[0]
    
    # Define model
    model = Model(name)
    
    # Add metabolites to the model
    print('Adding metabolites to sample ' + name)
    for i in range(nM):
        model.add_metabolites([Metabolite(metInfo['metID'][i], formula = metInfo['metName'][i],
                                          name = metInfo['metFormula'][i], compartment= metInfo['metComp'][i])])
    
    # Add reactions to the model
    mets = {met.id: met for met in model.metabolites}
    print('Adding reactions to sample ' + name)
    for i in range(nR):
        rxn = Reaction(rxnInfo.iloc[i]['rxnID'], name = rxnInfo.iloc[i]['rxnName'], subsystem = rxnInfo.iloc[i]['rxnSubs'],
                       lower_bound = rxnInfo.iloc[i]['lb'], upper_bound = rxnInfo.iloc[i]['ub'])
        rxn.add_metabolites({mets[m]: s for m, s in zip(rxnInfo['rxnEq'][i].split(';'), [float(x) for x in rxnInfo['rxnS'][i].split(';')])})
        model.add_reactions([rxn])

    # Return output variable
    return(model)

def _addInfoToModel(name: str, model: cobra.Model, metInfo: pd.DataFrame, rxnInfo: pd.DataFrame) -> (cobra.Model):
    # Rename model
    model.id = name
    
    # Number of metabolites and reactions
    nM = metInfo.shape[0]
    nR = rxnInfo.shape[0]
    
    # Add metabolites to the model
    print('Adding metabolites to sample ' + name)
    for i in range(nM):
        model.add_metabolites([Metabolite(metInfo.loc[i,'metID'], formula = metInfo.loc[i,'metName'],
                                          name = metInfo.loc[i,'metFormula'], compartment= metInfo.loc[i,'metComp'])])
    
    # Add reactions to the model
    mets = {met.id: met for met in model.metabolites}
    print('Adding reactions to sample ' + name)
    for i in range(nR):
        rxn = Reaction(rxnInfo.iloc[i]['rxnID'], name = rxnInfo.iloc[i]['rxnName'], subsystem = rxnInfo.iloc[i]['rxnSubs'],
                       lower_bound = rxnInfo.iloc[i]['lb'], upper_bound = rxnInfo.iloc[i]['ub'])
        rxn.add_metabolites({mets[m]: s for m, s in zip(rxnInfo['rxnEq'][i].split(';'), [float(x) for x in rxnInfo['rxnS'][i].split(';')])})
        model.add_reactions([rxn])

    # Return output variable
    return(model)

def _labelTaxonModel(sp: str, metInfo: pd.DataFrame, rxnInfo: pd.DataFrame, metTaxMat: pd.DataFrame, rxnTaxMat: pd.DataFrame):
    # Filter metabolites and reactions present in the taxon
    metInfo = metInfo[metTaxMat[sp] == 1]
    rxnInfo = rxnInfo[rxnTaxMat[sp] == 1]
    
    # Label metabolite abbreviation
    metInfo.loc[:,'metID'] = sp + '_' + metInfo.loc[:,'metID'].astype(str)
    
    # Rename extracellular compartment to lumen compartment
    metInfo.loc[:,'metID'] = metInfo.loc[:,'metID'].str.replace(r'_e$', '_u', regex = True)
    metInfo.loc[:,'metComp'] = metInfo.loc[:,'metComp'].str.replace(r'e', 'u', regex = True)
    
    # Label reaction abbreviation
    rxnInfo.loc[:,'rxnID'] = sp + '_' + rxnInfo.loc[:,'rxnID'].astype(str)
    
    # Rename reaction equation
    rxnInfo.loc[:,'rxnEq'] = rxnInfo.loc[:,'rxnEq'].apply(lambda x: ';'.join([sp + '_' + (m[:-2] + '_u' if m.endswith('_e') else m) for m in x.split(';')]))

    # Rename biomass reaction
    rxnInfo.loc[rxnInfo['rxnID'] == (sp + '_biomassPan_' + sp),'rxnID'] = sp + '_biomassPan'

    # Output variable
    return metInfo, rxnInfo

def _addCouplingConstraints(model, taxon, c: float = 1000.0, u: float = 0.01):
    # Identify objective reaction
    objRxn = [r.id for r in model.reactions if 'biomass' in r.id and taxon in r.id]
    objRxnFlux = model.reactions.get_by_id(objRxn[0]).flux_expression
    
    # Identify reactions to be constrained
    rxnList = [r.id for r in model.reactions if r.id.startswith(taxon) and r.id not in objRxn]
    
    # Add coupling constraints
    for r in rxnList:
        # Extract reaction information
        rxn = model.reactions.get_by_id(r)
        rxnFlux = rxn.flux_expression
        
        # Add forward direction constraint
        cons = model.problem.Constraint(rxnFlux - c * objRxnFlux, ub = u)
        model.add_cons_vars(cons)
        
        # Add backward direction constraint
        if rxn.lower_bound < 0:
            revCons = model.problem.Constraint(rxnFlux + c * objRxnFlux, lb = -u)
            model.add_cons_vars(revCons)
            
    # Output variable
    return(model)

def _addCouplingConstraintsWBM(model: cobra.Model, constraints: pd.DataFrame) -> (cobra.Model):
    # Number of constraints
    nC = constraints.shape[0]
    
    # Constrain reactions of each organ based on its biomass
    for i in range(nC):
        # Extract reaction IDs and coefficients
        consRxns = [model.reactions.get_by_id(x) for x in constraints.loc[i,'rxnID'].split(';')]
        consRxnFlux = [r.flux_expression for r in consRxns]
        consCoef = [float(x) for x in constraints.loc[i,'rxnCoef'].split(';')]
        
        # Build constraint linear expression
        constraint_expr = consCoef[0] * consRxnFlux[0]
        for coef, flux in zip(consCoef[1:], consRxnFlux[1:]):
            constraint_expr += coef * flux
        
        # Add coupling constraint
        if constraints.loc[i,'conSense'] == 'L':
            finalCons = model.problem.Constraint(constraint_expr, ub = constraints.loc[i, 'conBound'])
        elif constraints.loc[i,'conSense'] == 'G':
            finalCons = model.problem.Constraint(constraint_expr, lb = constraints.loc[i, 'conBound'])
        model.add_cons_vars(finalCons)
        
    # Output variable
    return (model)

def _adaptDiet(dietTable: pd.DataFrame, includeHumanMets: bool = True) -> (pd.DataFrame):
    # Remove all 0 values in the diet
    dietTable = dietTable[dietTable.iloc[:,1]!=0]
    
    # Transform uptakes to negative
    dietTable.iloc[:,1] = [-x for x in dietTable.iloc[:,1]]
    
    # Save original diet information
    origDiet = copy.deepcopy(dietTable)
    
    # Fix exchange nomenclature issues
    replacements = {'\[e\]': '(e)','EX_adpcbl\(e\)': 'EX_adocbl(e)','EX_glc\(e\)': 'EX_glc_D(e)','EX_sbt-d\(e\)': 'EX_sbt_D_(e)'}
    for pattern, repl in replacements.items():
        dietTable.iloc[:, 0] = dietTable.iloc[:, 0].replace(pattern, repl, regex=True)
    
    # Add missing essential metabolites
    essentialMetabolites = ['EX_12dgr180(e)', 'EX_26dap_M(e)', 'EX_2dmmq8(e)', 'EX_2obut(e)', 'EX_3mop(e)', 'EX_4abz(e)',
                            'EX_4hbz(e)', 'EX_ac(e)', 'EX_acgam(e)', 'EX_acmana(e)', 'EX_acnam(e)', 'EX_ade(e)', 'EX_adn(e)',
                            'EX_adocbl(e)', 'EX_ala_D(e)', 'EX_ala_L(e)', 'EX_amet(e)', 'EX_amp(e)', 'EX_arab_D(e)',
                            'EX_arab_L(e)', 'EX_arg_L(e)', 'EX_asn_L(e)', 'EX_btn(e)', 'EX_ca2(e)', 'EX_cbl1(e)',
                            'EX_cgly(e)', 'EX_chor(e)', 'EX_chsterol(e)', 'EX_cit(e)', 'EX_cl(e)', 'EX_cobalt2(e)',
                            'EX_csn(e)', 'EX_cu2(e)', 'EX_cys_L(e)', 'EX_cytd(e)', 'EX_dad_2(e)', 'EX_dcyt(e)',
                            'EX_ddca(e)', 'EX_dgsn(e)', 'EX_fald(e)', 'EX_fe2(e)', 'EX_fe3(e)', 'EX_fol(e)', 'EX_for(e)',
                            'EX_gal(e)', 'EX_glc_D(e)', 'EX_gln_L(e)', 'EX_glu_L(e)', 'EX_gly(e)', 'EX_glyc(e)',
                            'EX_glyc3p(e)', 'EX_gsn(e)', 'EX_gthox(e)', 'EX_gthrd(e)', 'EX_gua(e)', 'EX_h(e)', 'EX_h2o(e)',
                            'EX_h2s(e)', 'EX_his_L(e)', 'EX_hxan(e)', 'EX_ile_L(e)', 'EX_k(e)', 'EX_lanost(e)', 'EX_leu_L(e)',
                            'EX_lys_L(e)', 'EX_malt(e)', 'EX_met_L(e)', 'EX_mg2(e)', 'EX_mn2(e)', 'EX_mqn7(e)', 'EX_mqn8(e)',
                            'EX_nac(e)', 'EX_ncam(e)', 'EX_nmn(e)', 'EX_no2(e)', 'EX_ocdca(e)', 'EX_ocdcea(e)', 'EX_orn(e)',
                            'EX_phe_L(e)', 'EX_pheme(e)', 'EX_pi(e)', 'EX_pnto_R(e)', 'EX_pro_L(e)', 'EX_ptrc(e)', 'EX_pydx(e)',
                            'EX_pydxn(e)', 'EX_q8(e)', 'EX_rib_D(e)', 'EX_ribflv(e)', 'EX_ser_L(e)', 'EX_sheme(e)', 'EX_so4(e)',
                            'EX_spmd(e)', 'EX_thm(e)', 'EX_thr_L(e)', 'EX_thymd(e)', 'EX_trp_L(e)', 'EX_ttdca(e)', 'EX_tyr_L(e)',
                            'EX_ura(e)', 'EX_val_L(e)', 'EX_xan(e)', 'EX_xyl_D(e)', 'EX_zn2(e)', 'EX_glu_D(e)', 'EX_melib(e)',
                            'EX_chtbs(e)', 'EX_metsox_S_L(e)', 'EX_hdca(e)', 'EX_gam(e)', 'EX_indole(e)', 'EX_glcn(e)',
                            'EX_coa(e)', 'EX_man(e)', 'EX_fum(e)', 'EX_succ(e)', 'EX_no3(e)', 'EX_ins(e)', 'EX_uri(e)', 
                            'EX_drib(e)', 'EX_pime(e)', 'EX_lac_L(e)', 'EX_glypro(e)', 'EX_urea(e)', 'EX_duri(e)', 'EX_h2(e)',
                            'EX_mal_L(e)', 'EX_tre(e)', 'EX_orot(e)', 'EX_glymet(e)', 'EX_glyleu(e)', 'EX_pydx5p(e)',
                            'EX_so3(e)', 'EX_nh4(e)']
    missingMetabolites = [x for x in essentialMetabolites if x not in dietTable.iloc[:,0].values]
    diet = pd.concat([dietTable, pd.DataFrame([[m, -0.1] for m in missingMetabolites], columns=dietTable.columns)]).reset_index(drop = True)
    
    # Allow uptake of certain dietary compounds that are currently not mapped in the Diet Designer
    unmappedMets = ['EX_asn_L(e)','EX_gln_L(e)','EX_crn(e)','EX_elaid(e)','EX_hdcea(e)','EX_dlnlcg(e)','EX_adrn(e)',
                    'EX_hco3(e)','EX_sprm(e)', 'EX_carn(e)','EX_7thf(e)','EX_Lcystin(e)','EX_hista(e)','EX_orn(e)',
                    'EX_ptrc(e)','EX_creat(e)','EX_cytd(e)','EX_so4(e)']
    
    # Remove any missing metabolite in the original diet file
    # missingMetabolites = [x for x in unmappedMets if x not in dietTable.iloc[:,0].values]
    missingMetabolites = unmappedMets
    diet = diet[~diet.iloc[:, 0].isin(missingMetabolites)]
    diet = pd.concat([diet, pd.DataFrame([[m, -50] for m in missingMetabolites], columns=dietTable.columns)]).reset_index(drop = True)
    
    # Map choline flux
    if ~(diet.iloc[:,0] == 'EX_chol(e)').any():
        diet = pd.concat([diet, pd.DataFrame([['EX_chol(e)', -41.251]], columns = dietTable.columns)]).reset_index(drop = True)
        
    # Increase the uptake rate of micronutrients with too low defined uptake
    micronutrients = ['EX_adocbl(e)','EX_vitd2(e)','EX_vitd3(e)','EX_psyl(e)','EX_gum(e)','EX_bglc(e)','EX_phyQ(e)','EX_fol(e)',
                      'EX_5mthf(e)','EX_q10(e)','EX_retinol_9_cis(e)','EX_pydxn(e)','EX_pydam(e)','EX_pydx(e)','EX_pheme(e)',
                      'EX_ribflv(e)','EX_thm(e)','EX_avite1(e)','EX_pnto_R(e)','EX_na1(e)','EX_cl(e)','EX_k(e)','EX_pi(e)','EX_zn2(e)',
                      'EX_cu2(e)','EX_btn(e)']
    idx = (diet.iloc[:,0].isin(micronutrients)) & (abs(diet.iloc[:,1]) <= 0.1)
    diet.loc[idx,dietTable.columns[1]] = diet.loc[idx,dietTable.columns[1]] * 100
    
    # Pantothenate uptake needs to be at least 0.1
    idx = (diet.iloc[:,0].isin(['EX_pnto_R(e)'])) & (abs(diet.iloc[:,1]) < 0.1)
    diet.loc[idx,dietTable.columns[1]] = -0.1
    
    # Folate, L-arabinose, D-xylose, AMP, NH4 uptake need to be at least 1
    nutrients = ['EX_fol(e)','EX_arab_L(e)','EX_xyl_D(e)','EX_amp(e)','EX_nh4(e)','EX_cobalt2(e)']
    idx = (diet.iloc[:,0].isin(nutrients)) & (abs(diet.iloc[:,1]) < 1)
    diet.loc[idx,dietTable.columns[1]] = -1
    
    # Include human metabolites
    if includeHumanMets:
        humanMets = ['Diet_EX_gchola_d', 'Diet_EX_tdchola_d', 'Diet_EX_tchola_d', 'Diet_EX_dgchol_d', 'Diet_EX_34dhphe_d', 'Diet_EX_5htrp_d', 'Diet_EX_Lkynr_d',
                     'Diet_EX_f1a_d', 'Diet_EX_gncore1_d', 'Diet_EX_gncore2_d', 'Diet_EX_dsT_antigen_d', 'Diet_EX_sTn_antigen_d', 'Diet_EX_core8_d', 'Diet_EX_core7_d',
                     'Diet_EX_core5_d', 'Diet_EX_core4_d', 'Diet_EX_ha_d', 'Diet_EX_cspg_a_d', 'Diet_EX_cspg_b_d', 'Diet_EX_cspg_c_d', 'Diet_EX_cspg_d_d',
                     'Diet_EX_cspg_e_d', 'Diet_EX_hspg_d']
        humanValues = [-10, -10, -10, -10, -10, -10, -10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
        humanDiet = pd.DataFrame(list(zip(humanMets, humanValues)), columns=diet.columns)

        # Identify human metabolites not included
        idx = ~humanDiet.iloc[:, 0].isin(diet.iloc[:, 0])
        
        # Append only missing ones
        diet = pd.concat([diet, humanDiet[idx]], ignore_index=True)

    # Set upper bound of original diet information
    diet['ub'] = 0.0
    for ind, met in enumerate(origDiet.iloc[:,0]):
        diet.loc[diet.iloc[:,0] == met, 'ub'] = origDiet.iloc[ind,1] * 0.8    

    # Set upper bound of human metabolites
    diet.loc[diet.iloc[:,0].isin(humanMets),'ub'] = 10000
    
    # Adapt exchange IDs to microbiota models
    diet.iloc[:,0] = diet.iloc[:,0].str.replace(r'^EX_','Diet_EX_', regex = True)
    diet.iloc[:,0] = diet.iloc[:,0].str.replace(r'\(e\)$','_d', regex = True)
    
    # Output variable
    return (diet)

def _setDietConstraints(model: cobra.Model, diet: pd.DataFrame, factor: int = 1) -> (cobra.Model):
    # Define model reactions
    modRxns = [x.id for x in model.reactions]
    
    # Load AGORA essential metabolites
    essentialMetabolites = ['Diet_EX_12dgr180_d', 'Diet_EX_26dap_M_d', 'Diet_EX_2dmmq8_d', 'Diet_EX_2obut_d', 'Diet_EX_3mop_d', 'Diet_EX_4abz_d', 'Diet_EX_4hbz_d',
                            'Diet_EX_5HPET_d', 'Diet_EX_ac_d', 'Diet_EX_acgam_d', 'Diet_EX_ach_d', 'Diet_EX_acmana_d', 'Diet_EX_ade_d', 'Diet_EX_adn_d', 'Diet_EX_adocbl_d',
                            'Diet_EX_adpcbl_d', 'Diet_EX_ala_D_d', 'Diet_EX_ala_L_d', 'Diet_EX_alaasp_d', 'Diet_EX_alagln_d', 'Diet_EX_alahis_d', 'Diet_EX_alathr_d',
                            'Diet_EX_amet_d', 'Diet_EX_amp_d', 'Diet_EX_appnn_d', 'Diet_EX_arab_D_d', 'Diet_EX_arbt_d', 'Diet_EX_arg_L_d', 'Diet_EX_asn_L_d', 'Diet_EX_btn_d',
                            'Diet_EX_ca2_d', 'Diet_EX_cbl1_d', 'Diet_EX_cgly_d', 'Diet_EX_chor_d', 'Diet_EX_chsterol_d', 'Diet_EX_chtbs_d', 'Diet_EX_cit_d', 'Diet_EX_cl_d',
                            'Diet_EX_coa_d', 'Diet_EX_cobalt2_d', 'Diet_EX_csn_d', 'Diet_EX_cu2_d', 'Diet_EX_cys_L_d', 'Diet_EX_cytd_d', 'Diet_EX_dad_2_d', 'Diet_EX_dcyt_d',
                            'Diet_EX_ddca_d', 'Diet_EX_dextrin_d', 'Diet_EX_dgsn_d', 'Diet_EX_fald_d', 'Diet_EX_fe2_d', 'Diet_EX_fe3_d', 'Diet_EX_fol_d', 'Diet_EX_for_d',
                            'Diet_EX_gal_d', 'Diet_EX_gam_d', 'Diet_EX_gd1c_hs_d', 'Diet_EX_glc_D_d', 'Diet_EX_glcn_d', 'Diet_EX_gln_L_d', 'Diet_EX_glu_L_d', 'Diet_EX_gly_d',
                            'Diet_EX_glyasn_d', 'Diet_EX_glyc3p_d', 'Diet_EX_glyc_d', 'Diet_EX_glygn5_d', 'Diet_EX_glyleu_d', 'Diet_EX_glymet_d', 'Diet_EX_glytyr_d',
                            'Diet_EX_gsn_d', 'Diet_EX_gthox_d', 'Diet_EX_gthrd_d', 'Diet_EX_gua_d', 'Diet_EX_h2_d', 'Diet_EX_h2o_d', 'Diet_EX_h2s_d', 'Diet_EX_h_d',
                            'Diet_EX_his_L_d', 'Diet_EX_hxan_d', 'Diet_EX_ile_L_d', 'Diet_EX_indole_d', 'Diet_EX_k_d', 'Diet_EX_ksi_d', 'Diet_EX_lanost_d', 'Diet_EX_leu_L_d',
                            'Diet_EX_lys_L_d', 'Diet_EX_malt_d', 'Diet_EX_mantr_d', 'Diet_EX_melib_d', 'Diet_EX_met_L_d', 'Diet_EX_metala_d', 'Diet_EX_metsox_S_L_d',
                            'Diet_EX_mg2_d', 'Diet_EX_mn2_d', 'Diet_EX_mnl_d', 'Diet_EX_mqn7_d', 'Diet_EX_mqn8_d', 'Diet_EX_nac_d', 'Diet_EX_ncam_d', 'Diet_EX_nmn_d',
                            'Diet_EX_no2_d', 'Diet_EX_no3_d', 'Diet_EX_ocdca_d', 'Diet_EX_ocdcea_d', 'Diet_EX_orn_d', 'Diet_EX_phe_L_d', 'Diet_EX_pheme_d', 'Diet_EX_pi_d',
                            'Diet_EX_pime_d', 'Diet_EX_pnto_R_d', 'Diet_EX_pro_L_d', 'Diet_EX_prostge1_d', 'Diet_EX_ptrc_d', 'Diet_EX_pydx5p_d', 'Diet_EX_pydx_d',
                            'Diet_EX_pydxn_d', 'Diet_EX_q8_d', 'Diet_EX_rib_D_d', 'Diet_EX_ribflv_d', 'Diet_EX_ser_L_d', 'Diet_EX_sheme_d', 'Diet_EX_so4_d', 'Diet_EX_spmd_d',
                            'Diet_EX_sucr_d', 'Diet_EX_thm_d', 'Diet_EX_thr_L_d', 'Diet_EX_thymd_d', 'Diet_EX_tre_d', 'Diet_EX_trp_L_d', 'Diet_EX_ttdca_d', 'Diet_EX_tyr_L_d',
                            'Diet_EX_ura_d', 'Diet_EX_val_L_d', 'Diet_EX_xan_d', 'Diet_EX_xyl_D_d', 'Diet_EX_zn2_d']
    
    # Block all diet uptake reactions
    dietRxns = [x for x in modRxns if x.startswith('Diet_EX_')]
    for reaction in dietRxns:
        rxn = model.reactions.get_by_id(reaction)
        rxn.upper_bound = 0
        rxn.lower_bound = 0
    
    # Ensure uptake of metabolites required for microbiota
    missingUptakes = [x for x in essentialMetabolites if x not in diet.iloc[:,0].values]
    for reaction in missingUptakes:
        rxn = model.reactions.get_by_id(reaction)
        rxn.lower_bound = -0.1
        
    # Add more required compounds
    missingDietCompounds = ['Diet_EX_asn_L_d', 'Diet_EX_gln_L_d', 'Diet_EX_chol_d', 'Diet_EX_crn_d', 'Diet_EX_elaid_d', 'Diet_EX_hdcea_d', 'Diet_EX_dlnlcg_d',
                            'Diet_EX_adrn_d', 'Diet_EX_hco3_d', 'Diet_EX_sprm_d', 'Diet_EX_carn_d', 'Diet_EX_7thf_d', 'Diet_EX_Lcystin_d', 'Diet_EX_hista_d', 'Diet_EX_orn_d',
                            'Diet_EX_ptrc_d', 'Diet_EX_creat_d', 'Diet_EX_cytd_d', 'Diet_EX_so4_d']
    for reaction in missingDietCompounds:
        rxn = model.reactions.get_by_id(reaction)
        rxn.lower_bound = -50
        
    # Adjust choline uptake
    rxn = model.reactions.get_by_id('Diet_EX_chol_d')
    rxn.lower_bound = -41.251
    
    # Define additional groups of metabolites
    microNutrients = ['Diet_EX_adocbl_d', 'Diet_EX_vitd2_d', 'Diet_EX_vitd3_d', 'Diet_EX_psyl_d', 'Diet_EX_gum_d', 'Diet_EX_bglc_d', 'Diet_EX_phyQ_d', 'Diet_EX_fol_d',
                      'Diet_EX_5mthf_d', 'Diet_EX_q10_d', 'Diet_EX_retinol_9_cis_d', 'Diet_EX_pydxn_d', 'Diet_EX_pydam_d', 'Diet_EX_pydx_d', 'Diet_EX_pheme_d',
                      'Diet_EX_ribflv_d', 'Diet_EX_thm_d', 'Diet_EX_avite1_d', 'Diet_EX_pnto_R_d']
    ions = ['Diet_EX_na1_d', 'Diet_EX_cl_d', 'Diet_EX_k_d', 'Diet_EX_pi_d', 'Diet_EX_zn2_d', 'Diet_EX_cu2_d']
    so4 = ['Diet_EX_so4_d']
    
    # Add diet uptakes
    for reaction, flux in zip(diet.iloc[:,0], diet.iloc[:,1]):
        rxn = model.reactions.get_by_id(reaction)
        if reaction in microNutrients and flux <= 10:
            rxn.lower_bound = -10 * factor
        elif reaction in microNutrients and flux > 0.1:
            rxn.lower_bound = -1.2 * 100 * factor * flux
        elif reaction in ions:
            rxn.lower_bound = -1.2 * 100 * factor * flux
        elif reaction in so4:
            rxn.lower_bound = -1000 * factor
        else:
            rxn.lower_bound = -1.2 * factor * flux
            
    # Force diet uptakes
    for reaction, flux in zip(diet.iloc[:,0], diet.iloc[:,1]):
        rxn = model.reactions.get_by_id(reaction)
        if reaction in microNutrients:
            rxn.upper_bound = -0.8 * factor * flux
        elif reaction in ions:
            rxn.upper_bound = -0.8 * factor * flux
        else:
            rxn.upper_bound = -0.8 * factor * flux
            
    # Output variable
    return (model)
             
def _changeRxnBoundsToModel(model: cobra.Model, rxnList: list, rxnFlux: list, boundType: str = 'lb') -> (cobra.Model):
    # Check bound type value
    if not boundType in ['lb', 'ub']:
        raise ValueError(f'Invalid bound: {boundType}. Must be one of [lb, ub].')
    
    # Filter reactions
    modelRxns = [x.id for x in model.reactions]
    filtRxns, filtFlux = zip(*[(r, f) for r, f in zip(rxnList, rxnFlux) if r in modelRxns])
    
    # Check if any reaction is present
    if not filtRxns:
        print('There are no common reactions between the input and the database reactions')
        return (model)
    
    # Change the bound of each reaction in the model
    for i in range(len(filtRxns)):
        rxn = model.reactions.get_by_id(filtRxns[i])
        if boundType == 'lb':
            rxn.lower_bound = filtFlux[i]
        elif boundType == 'ub':
            rxn.upper_bound = filtFlux[i]
    
    # Output variable
    return (model)

def _useDiet(model: cobra.Model, dietTable = pd.DataFrame) -> (cobra.Model):
    # Extract diet information
    rxnList = dietTable.iloc[:,0].to_list()
    rxnLB = dietTable.iloc[:,1].to_list()
    
    # Identify diet exchange reactions
    modelRxns = [x.id for x in model.reactions]
    dietRxns = [x for x in modelRxns if x.startswith('Diet_EX_')]
    
    # Set diet reactions uptake to 0
    model = _changeRxnBoundsToModel(model, dietRxns, [0] * len(dietRxns), 'lb')
    
    # Change input diet reactions uptake
    model = _changeRxnBoundsToModel(model, rxnList, rxnLB)
    
    # Change input diet reactions secretion
    if dietTable.shape[1] > 2:
        model = _changeRxnBoundsToModel(model, rxnList, dietTable.iloc[:,2].to_list(), 'ub')
    
    # Output variable
    return (model)
