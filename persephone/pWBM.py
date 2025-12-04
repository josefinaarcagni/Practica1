# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 23:26:45 2025

@author: tblasco
"""

import pandas as pd
import os
import re
import cobra
import pickle
import importlib.resources as pkg_resources
from functools import partial
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

from .WBM import _generateOriginalWBM

def _defaultIndividualParameters(sex: str) -> (dict, pd.DataFrame, list):

    if sex == 'female':
        # Define individual parameters
        indParam = {'ID': 'Default',
            'bodyWeight': 58,
            'Height': 160,
            'sex': 'female',
            'HeartRate': 67,
            'StrokeVolume': 80,
            'CardiacOutput': 5360,
            'Hematocrit': 0.4,
            'MConUrCreatinineMax': 1.2,
            'MConUrCreatinineMin': 0.5,
            'MConDefaultBc': 20,
            'MConDefaultCSF': 20,
            'MConDefaultUrMax': 20,
            'MConDefaultUrMin': 0,
            'CSFFlowRate': 0.35,
            'CSFBloodFlowRate': 0.52,
            'UrFlowRate': 2000,
            'GlomerularFiltrationRate': 90}
        
        # Define blood flow data
        bloodData = pd.DataFrame({'Organ': ['Agland', 'Adipocytes', 'Bones', 'Brain', 'Breast', 'Cervix', 'Colon', 'Esophagus', 'Gall', 'Heart', 'Kidney',
                                            'Liver', 'Lung', 'Muscle', 'Ovary', 'Pthyroidgland', 'Pancreas', 'Prostate', 'Rectum', 'Retina', 'Scord',
                                            'sIEC', 'Skin', 'Spleen', 'Stomach', 'Testis', 'Thyroidgland', 'Urinarybladder', 'Uterus', 'Total blood flow'],
                                  'Blood Flow Rate Min': [None, None, None, 0.63, None, None, 0.273, None, None, None, None, None, None, None, None, 0.045,
                                                          0.05, None, 0.121333333, None, None, 0.273, None, 0.2, 0.121333333, None, 0.045, None, None, 5.8],
                                 'Blood Flow Percentage': [0.000909091, 0.05, 0.05, 0.10862069, 0.000909091, 0.000909091, 0.047068966, 0.005, 0.000909091,
                                                           0.05, 0.2, 0.1, 0.025, 0.15, 0.000909091, 0.007758621, 0.00862069, 0.000909091,
                                                           0.02091954, 0.000909091, 0.000909091, 0.047068966, 0.05, 0.034482759, 0.02091954, 0.000909091,
                                                           0.007758621, 0.000909091, 0.000909091, 0.9923093]})
        
        # Define organ list
        organsListExt = ['Heart', 'Muscle', 'Lung', 'Skin', 'Stomach', 'sIEC', 'Colon', 'Urinarybladder', 'Retina', 'Scord', 'Brain', 'Adipocytes', 'Liver', 'Gall',
                        'Kidney', 'Pancreas', 'Spleen', 'Agland', 'Thyroidgland', 'Pthyroidgland', 'Ovary', 'Uterus', 'Breast', 'Cervix', 'Bcells', 'CD4Tcells',
                        'Nkcells', 'Monocyte', 'Platelet', 'RBC', 'BBB', 'Diet', 'SI', 'GI', 'LI', 'BileDuct', 'Excretion']
    
    elif sex == 'male':
        # Define individual parameters
        indParam = {'ID': 'Default',
            'bodyWeight': 70,
            'Height': 170,
            'sex': 'male',
            'HeartRate': 67,
            'StrokeVolume': 80,
            'CardiacOutput': 5360,
            'Hematocrit': 0.4,
            'MConUrCreatinineMax': 1.2,
            'MConUrCreatinineMin': 0.5,
            'MConDefaultBc': 20,
            'MConDefaultCSF': 20,
            'MConDefaultUrMax': 20,
            'MConDefaultUrMin': 0,
            'CSFFlowRate': 0.35,
            'CSFBloodFlowRate': 0.52,
            'UrFlowRate': 2000,
            'GlomerularFiltrationRate': 90}
        
        # Define blood flow data
        bloodData = pd.DataFrame({'Organ': ['Agland', 'Adipocytes', 'Bones', 'Brain', 'Breast', 'Cervix', 'Colon', 'Esophagus', 'Gall', 'Heart', 'Kidney',
                                            'Liver', 'Lung', 'Muscle', 'Ovary', 'Pthyroidgland', 'Pancreas', 'Prostate', 'Rectum', 'Retina', 'Scord',
                                            'sIEC', 'Skin', 'Spleen', 'Stomach', 'Testis', 'Thyroidgland', 'Urinarybladder', 'Uterus', 'Total blood flow'],
                                  'Blood Flow Rate Min': [None, None, None, 0.68, None, 0, 0.279, None, None, None, 1.17, None, None, None, None, 0.05,
                                                          0.06, None, 0.124, None, None, 0.279, None, 0.2, 0.124, None, 0.05, None, None, 5.8],
                                 'Blood Flow Percentage': [0.000909091, 0.05, 0.05, 0.117241379, 0.000909091, 0.000909091, 0.048103448, 0.005, 0.000909091,
                                                           0.05, 0.201724138, 0.1, 0.025, 0.15, 0.000909091, 0.00862069, 0.010344828, 0.000909091,
                                                           0.02137931, 0.000909091, 0.000909091, 0.048103448, 0.05, 0.034482759, 0.02137931, 0.000909091,
                                                           0.00862069, 0.000909091, 0.000909091, 1.009090909]})
        
        # Define organ list
        organsListExt = ['Heart', 'Muscle', 'Lung', 'Skin', 'Stomach', 'sIEC', 'Colon', 'Urinarybladder', 'Retina', 'Scord', 'Brain', 'Adipocytes', 'Liver', 'Gall',
                         'Kidney', 'Pancreas', 'Spleen', 'Agland', 'Thyroidgland', 'Pthyroidgland', 'Testis', 'Prostate', 'Bcells', 'CD4Tcells', 'Nkcells', 'Monocyte',
                         'Platelet', 'RBC', 'BBB', 'Diet', 'SI', 'GI', 'LI', 'BileDuct', 'Excretion']
        
    # Output variable
    return (indParam, bloodData, organsListExt)

def _defaultOrganWeights(sex: str) -> (list, dict, dict):
    
    if sex == 'male':
        organNames = ['Adipocytes', 'Agland', 'Brain', 'Breast', 'Colon', 'CSF', 'Esophagus', 'Heart', 'Kidney', 'Liver', 'Lung', 'Muscle', 'Ovary', 'Pancreas', 'Placenta',
                      'Prostate', 'Pthyroidgland', 'Rectum', 'Retina', 'Scord', 'sIEC', 'Skin', 'Spleen', 'Stomach', 'Testis', 'Thyroidgland', 'Urinarybladder', 'Urine in bladder',
                      'Uterus', 'Blood', 'WBC', 'Lymphocytes', 'Bcells', 'CD4Tcells', 'CD8Tcells', 'Nkcells', 'Monocyte', 'Platelet', 'RBC', 'Skeleton (bone, bone marrow other tissue)',
                      'Tendrons etc', 'Cervix', 'Gall', 'Gut', 'Salvary glands', 'Ureter', 'Urethra', 'Teeth', 'Nails', 'Nose mucosa', 'Hair',
                      'Total body weigth captured by Harvey when excluding bones', 'Total body weigth captured by Harvey with bones']
        
        organWeight = [15000, 14, 1400, 26, 300, 120, 40, 331, 310, 1800, 536, 28000, 0, 100, 0, 16, 0.12, 70, 0.326, 30, 640, 2600, 180, 150, 35, 20, 45, 350, 0, 5500, 55, 16.5,
                       1.485, 8.25, 4.4, 2.475, 2.75, 20, 2475, 10000, 850, 0, 10, 1005, 85, 16, 10, 46, 3, 32, 20, 58628.446, 69690.446]
        organWeight = dict(zip(organNames, organWeight))
        
        organWeightFract = [0.21429, 0.0002, 0.02, 0.00037143, 0.0042857, 0.0017143, 0.00057143, 0.0047286, 0.0044286, 0.025714, 0.0076571, 0.4, 0, 0.0014286, 0, 0.00022857, 1.7143e-06,
                           0.001, 4.6571e-06, 0.00042857, 0.0091429, 0.037143, 0.0025714, 0.0021429, 0.0005, 0.00028571, 0.00064286, 0.005, 0, 0.078571, 0.00078571, 0.00023571,
                           2.1214e-05, 0.00011786, 6.2857e-05, 3.5357e-05, 3.9286e-05, 0.00028571, 0.035357, 0.14286, 0.012143, 0, 0.00014286, 0.014357, None, None, None, None, None,
                           None, None, 0.83755, 0.99558]
        organWeightFract = dict(zip(organNames, organWeightFract))
        
    elif sex == 'female':
        organNames = ['Adipocytes', 'Agland', 'Brain', 'Breast', 'Colon', 'CSF', 'Esophagus', 'Heart', 'Kidney', 'Liver', 'Lung', 'Muscle', 'Ovary', 'Pancreas', 'Placenta',
                      'Prostate', 'Pthyroidgland', 'Rectum', 'Retina', 'Scord', 'sIEC', 'Skin', 'Spleen', 'Stomach', 'Testis', 'Thyroidgland', 'Urinarybladder', 'Urine in bladder',
                      'Uterus', 'Blood', 'WBC', 'Lymphocytes', 'Bcells', 'CD4Tcells', 'CD8Tcells', 'Nkcells', 'Monocyte', 'Platelet', 'RBC', 'Skeleton (bone, bone marrow other tissue)',
                      'Tendrons etc', 'Cervix', 'Gall', 'Gut', 'Salvary glands', 'Ureter', 'Urethra', 'Teeth', 'Nails', 'Nose mucosa', 'Hair',
                      'Total body weigth captured by Harvey when excluding bones', 'Total body weigth captured by Harvey with bones']
        
        organWeight = [19000, 14, 1200, 360, 290, 100, 34, 253, 275, 1400, 536, 17000, 11, 85, 528, None, 0.14, 70, 0.326, 28, 600, 1790, 150, 140, 0, 17, 45, 350, 80, 4100, 41, 12.3,
                       1.107, 6.15, 3.28, 1.845, 2.05, 20, 2050, 6800, 700, 50, 8, 1005, 70, 15, 3, 41, 3, 27, 300, 49519.466, 57478.466]

        organWeight = dict(zip(organNames, organWeight))
        
        organWeightFract = [0.32759, 0.00024138, 0.02069, 0.0062069, 0.005, 0.0017241, 0.00058621, 0.0043621, 0.0047414, 0.024138, 0.0092414, 0.2931, 0.00018966, 0.0014655, 0.0091034,
                            0, 2.4138e-06, 0.0012069, 5.6207e-06, 0.00048276, 0.010345, 0.030862, 0.0025862, 0.0024138, 0, 0.0002931, 0.00077586, 0.0060345, 0.0013793, 0.07069,
                            0.0007069, 0.00021207, 1.9086e-05, 0.00010603, 5.6552e-05, 3.181e-05, 3.5345e-05, 0.00034483, 0.035345, 0.11724, 0.012069, 0.00086207, 0.00013793,
                            0.017328, None, None, None, None, None, None, None, 0.85378, 0.99101]
        organWeightFract = dict(zip(organNames, organWeightFract))
    
    # Output variable
    return (organNames, organWeight, organWeightFract)

def _calcOrganFract(model: cobra.Model, indParam: dict) -> (dict):
    # Obtain reference organ weights
    refOrgans, refOrganWeight, refOrganWeightFract = _defaultOrganWeights(indParam['sex'])
    
    # Define sample body weight
    wt = indParam['bodyWeight']
    
    # Define polynomials for organ weights
    if indParam['sex'] == 'male':
        OF = {'Brain': [1.41e-01, -5.54e-06, 9.30e-11, -6.83e-16, 1.80e-21, 0.0], 'Heart': [6.32e-03, -1.67e-08, 0.0, 0.0, 0.0, 0.0], 'Kidney': [7.26e-03, -6.69e-08, 3.33e-13, 0.0, 0.0, 0.0],
              'Liver': [4.25e-02, -1.01e-06, 1.99e-11, -1.66e-16, 4.83e-22, 0.0], 'Lung': [1.86e-02, -4.55e-08, 0.0, 0.0, 0.0, 0.0], 'Spleen': [3.12e-03, -5.57e-09, 0.0, 0.0, 0.0, 0.0],
              'Agland': [8.04e-04, -1.98e-08, 2.01e-13, -6.11e-19, 0.0, 0.0], 'Pancreas': [1.48e-03, 0.0, 0.0, 0.0, 0.0, 0.0], 'Thymus': [3.70e-03, -1.05e-07, 7.94e-13, 0.0, 0.0, 0.0],
              'Thyroidgland': [2.42e-04, 0.0, 0.0, 0.0, 0.0, 0.0], 'Adipocytes': [1.61e-01, -3.59e-06, 8.28e-11, -3.57e-16, 4.73e-22, 0.0], 'Muscle': [9.68e-02, -3.32e-06, 1.83e-10, -1.24e-15, 0.0, 0.0],
              'Skin': [1.03e-01, -2.56e-06, 3.68e-11, -2.58e-16, 8.62e-22, -1.10e-27], 'Blood': [8.97e-02, -3.50e-07, 6.54e-13, 0.0, 0.0, 0.0]}
    elif indParam['sex'] == 'female':
        OF = {'Brain': [1.12e-01, -3.33e-06, 4.30e-11, -2.45e-16, 5.03e-22, 0.0], 'Heart': [5.40e-03, -1.07e-08, 0.0, 0.0, 0.0, 0.0], 'Kidney': [7.56e-03, -5.58e-08, 1.54e-13, 0.0, 0.0, 0.0],
              'Liver': [3.34e-02, -1.89e-07, 5.34e-13, 0.0, 0.0, 0.0], 'Lung': [1.89e-02, -5.94e-08, 0.0, 0.0, 0.0, 0.0], 'Spleen': [2.96e-03, -7.72e-09, 0.0, 0.0, 0.0, 0.0],
              'Agland': [8.04e-04, -1.98e-08, 2.01e-13, -6.11e-19, 0.0, 0.0], 'Pancreas': [1.48e-03, 0.0, 0.0, 0.0, 0.0, 0.0], 'Thymus': [3.70e-03, -1.05e-07, 7.94e-13, 0.0, 0.0, 0.0],
              'Thyroidgland': [2.42e-04, 0.0, 0.0, 0.0, 0.0, 0.0], 'Adipocytes': [1.84e-01, -6.86e-06, 2.46e-10, -2.11e-15, 7.58e-21, -9.94e-27], 'Muscle': [3.65e-02, 7.91e-06, -5.74e-11, 0.0, 0.0, 0.0],
              'Skin': [9.81e-02, -2.28e-06, 2.74e-11, -1.58e-16, 4.30e-22, -4.43e-28], 'Blood': [8.97e-02, -3.50e-07, 6.54e-13, 0.0, 0.0, 0.0], 'Breast': [0.01, 0.0, 0.0, 0.0, 0.0, 0.0]}
    
    # Calculate individual specific organ weights
    organs = list(OF.keys())
    organWeight = {}
    organWeightFract = {}
    for org in organs:
        organWeightFract[org] = OF[org][0] + OF[org][1]*wt + OF[org][2]*(wt**2) + OF[org][3]*(wt**3) + OF[org][4]*(wt**4) + OF[org][5]*(wt**5)
        organWeight[org] = organWeightFract[org]*wt
        
    # Define blood volume
    indParam['BloodVolume'] = organWeight['Blood'] / 1.0506
    
    # Calculate individual blood cells
    WBCWeight = 0.01 * organWeight['Blood']; # Whole body cells
    LympWeight = 0.3 * WBCWeight; # Lymphocytes
    
    organs.append('Bcells')
    organWeight['Bcells'] = 0.09 * LympWeight
    organWeightFract['Bcells'] = organWeight['Bcells'] / wt
    
    organs.append('CD4Tcells')
    organWeight['CD4Tcells'] = 0.15 * WBCWeight
    organWeightFract['CD4Tcells'] = organWeight['CD4Tcells'] / wt
    
    organs.append('Nkcells')
    organWeight['Nkcells'] = 0.15 * LympWeight
    organWeightFract['Nkcells'] = organWeight['Nkcells'] / wt

    organs.append('Monocyte')
    organWeight['Monocyte'] = 0.05 * WBCWeight
    organWeightFract['Monocyte'] = organWeight['Monocyte'] / wt
    
    organs.append('Platelet')
    organWeight['Platelet'] = 4 * indParam['BloodVolume'] / 1000
    organWeightFract['Platelet'] = organWeight['Platelet'] / wt
    
    organs.append('RBC')
    organWeight['RBC'] = 495 * indParam['BloodVolume'] / 1000
    organWeightFract['RBC'] = organWeight['RBC'] / wt
    
    # Get all organs defined in the model
    modRxns = [x.id for x in model.reactions]
    objectiveComponents = [x for x in modRxns if '_biomass_maintenance' in x]
    objectiveComponents.append('sIEC_biomass_reactionIEC01b')
    organsInModel = [o.split('_')[0] for o in objectiveComponents]

    # Readjust organ weight fraction for missing organs
    missingOrgans = [x for x in organsInModel if x not in organs]
    for org in missingOrgans:
        organs.append(org)
        organWeight[org] = refOrganWeight[org]
        organWeightFract[org] = refOrganWeightFract[org] / wt
        
    # Set data into individual parameters variable
    indParam['organs'] = organs
    indParam['organWeight'] = organWeight
    indParam['organWeightFract'] = organWeightFract
    
    # Output variable
    return (indParam)

def _defaultMetabolitesConcentration(comp: str) -> (pd.DataFrame):
    
    if comp == 'blood':
        # Define metabolites concentration
        metConcentration = pd.DataFrame({'Metabolite': ['C09642', 'melatn', 'tym', 'CE1401', 'C05767', '34dhpha', 'CE2006', 'gal1p', '23dpg', '12ppd_R', 'xylu_L',
                                                        '2hyoxplac', 'xylu_D', 'ppbng', 'C05770', 'C05769', 'no', 'C06314', 'leuktrB4', 'estradiolglc', 'prostgf2',
                                                        'leuktrE4', 'prostgd2', 'leuktrC4', 'CE2537', 'C14769', 'aflatoxin', 'normete_L', 'C05957', 'estrones',
                                                        'vitd2', '15HPET', 'prostge2', '5fthf', 'nrpphr', '5adtststerone', 'CE7172', 'thf', 'eicostet', 'phyQ',
                                                        'C04717', 'C14826', 'C14768', 'fmn', 'C14825', '25hvitd3', 'prgstrn', 'vitd3', 'tetdece1crn', 'lnlccrn',
                                                        '25hvitd2', '34dhoxpeg', 'fad', 'c3dc', 'c4crn', 'debrisoquine', 'ptrc', 'homoval', 'thymd', 'c10crn',
                                                        'c81crn', 'retn', 'retinal', 'cytd', 'tchola', '34hpp', 'crtsl', 'urcan', 'strdnc', 'pydx5p', 'gchola',
                                                        'srtn', 'hexc', 'duri', 'gsn', 'gly', 'lac_L', 'clpnd', 'ile_L', '4hpro_LT', 'tsul', 'C11695', 'imp',
                                                        'glu_L', 'sprm', 'spmd', 'carn', 'crvnc', 'octa', 'bhb', 'chsterol', '4abut', 'C02528', 'xan', 'nrvnc',
                                                        'ala_B', 'bz', 'arach', 'csn', 'lthstrl', 'acetone', 'pyr', 'mthgxl', 'taur', 'acnam', 'ptdca', 'caro',
                                                        'CE2510', 'tyr_L', 'dlnlcg', 'lnlnca', 'gudac', 'ethamp', 'udpg', 'bgly', 'gal', 'ppa', 'lnlncg', 'sucr',
                                                        'alltn', 'C04805', 'tmndnc', 'no2', 'ump', 'gthrd', 'for', 'xylt', 'orot', 'hxan', 'ind3ac', 'atp', 'ttdca',
                                                        'co', 'ser_D', 'nh4', 'hdcea', 'tcynt', 'adn', 'etha', 'arachd', 'estradiol', 'CE0955', 'dcyt', 'gam',
                                                        '4aabutn', 'CE1297', 'icit', 'subeac', 'acac', 'Lpipecol', 'mal_L', 'glyc_R', '3mob', 'fuM', '3hpp', 'rib_D',
                                                        'rbt', 'L2aadp', 'arab_L', 'abt', 'glyc_S', '4mop', 'C08261', 'etoh', 'phyt', 'bvite', 'CE2028', 'succ',
                                                        'sarcs', 'ppp9', 'T4hcinnm', '4hbz', 'hdd2crn', 'mlthf', 'dctp', 'estrone', 'CE7081', 'andrstndn', 'CE2049',
                                                        'dhea', 'prgnlone', '12harachd', 'C13856', 'r5p', 'g6p', 'prist', 'frdp', 'HC02202', 'HC02203', 'dchac',
                                                        'doco13ac', 'tdechola', 'pcholn204_hs', 'pcholn24_hs', 'prostge1', 'thyox_L', 'aldstrn', 'C14771',
                                                        '1a25dhvitd3', '17ahprgstrn', 'leuktrD4', 'CE6205', 'aqcobal', 'CE7085', 'CE2445', 'CE7096', 'fna5moxam',
                                                        'CE7083', 'C05298', 'CE2211', 'CE4877', '6hoxmelatn', 'leuktrB4woh', 'HC02213', 'CE7090', 'andrstrnglc',
                                                        'leuktrF4', 'txa2', 'mepi', 'andrstandn', 'mhista', '18harachd', 'adrnl', '3ityr_L', 'CE1918', 'hista',
                                                        'CE2209', 'tststerone', 'andrstrn', 'pyam5p', 'C05300', 'xol24oh', 'leuktrB4wcooh', 'triodthy', 'btn',
                                                        'aprgstrn', 'xol25oh', 'anth', 'ahandrostan', '13_cis_retn', '17ahprgnlone', '35diotyr', 'CE1447', 'CE1243', 
                                                        '3moxtyr', '13_cis_retnglc', 'CE1617', 'estriol', 'ttdcrn', 'fol', '4pyrdx', 'dhf', 'kynate', 'cynt',
                                                        '34dhoxmand', 'N1aspmd', '11docrtsl', 'CE5072', 'dopasf', '35cgmp', 'C03681', 'nrpphrsf', 'HC02193',
                                                        'ddeccrn', 'gmp', 'CE1925', '34dhphe', 'ribflv', 'CE2705', 'c8crn', 'pydxn', '1mncam', 'CE2047', 'camp',
                                                        'thbpt', 'C05302', '13dampp', 'ncam', 'pcrn', 'q10', 'prgnlones', 'cpppg1', 'sphgn', 'xol27oh', 'sph1p',
                                                        '5htrp', 'hgentis', 'omeprazole', 'c6crn', '3hanthrn', 'mev_R', 'crtstrn', '3mox4hoxm', 'HC00900', 'C02470',
                                                        'stcrn', '56dura', 'limnen', '2mcit', 'n8aspmd', 'cortsn', 'CE2176', 'c51crn', 'ins', 'ivcrn', '11docrtstrn',
                                                        '5g2oxpt', '5hoxindoa', 'dgchol', 'HC02191', 'sphs1p', 'lipoate', 'sphings', 'adrn', 'galt', 'hpdca', '3mlda',
                                                        'thm', 'ura', 'cyst_L', 'adpac', '4hphac', 'pmtcrn', 'HC02172', 'pydam', 'c101crn', 'im4ac', '5mthf', 'ade',
                                                        'CE4968', 'ddca', 'retinol', 'tdchola', 'docosac', 'pydx', '34hpl', '5aop', '4hdebrisoquine', 'cholp', 'egme',
                                                        'minohp', 'C10164', 'lgnc', 'cholate', 'iodine', 'glx', 'lanost', 'HC00460', 'dsmsterol', 'phpyr', 'but',
                                                        'acorn', 'dudp', 'cys_L', '3aib', 'sbt_D', 'mma', 'quln', 'pcholhep_hs', 'fdp', 'ahcys', 'CE0737',
                                                        '4tmeabutn', 'uri', 'ttdcea', 'selmeth', 'HC02192', 'acald', '7dhchsterol', 'yvite', 'ppi', 'pnto_R', 'Lkynr',
                                                        'prpp', 'CE6031', 'glucys', 'gthox', 'dheas', 'dmgly', 'oxa', 'acrn', 'avite1', 'sql', 'CE1352', 'g3p',
                                                        'hcys_L', 'pep', 'pcollg5hlys', 'dhap', 'met_L', 'asp_D', 'pcholole_hs', 'glcur', 'ascb_L', 'pser_L', 'hom_L',
                                                        'ac', 'glc_D', 'co2', 'cdp', '5oxpro', 'ca2', 'dopa', 'fald', 'adp', 'na1', 'nadh', 'ala_L', 'citr_L', 'nac',
                                                        'C02356', 'hco3', 'fru', 'asn_L', 'glyald', 'udp', 'hdca', 'urea', 'HC02020', 'g3pc', 'dodecanac', 'glyc',
                                                        'gdp', 'pheacgln', 'amp', 'avite2', 'leu_L', 'glyb', 'aact', 'so4', 'cit', 'ocdca', 'nad', 'glyc3p', 'man',
                                                        'trp_L', 'crn', 'meoh', 'crtn', 'HC00250', 'gln_L', 'o2', '3pg', 'glyclt', 'h2o2', 'pcholmyr_hs', 'akg',
                                                        '2obut', 'chol', 'val_L', 'ser_L', 'orn_D', 'nadph', 'k', 'his_L', 'phe_L', 'glcn', '3hmp', 'CE1310',
                                                        'HC00319', 'chsterols', 'xtsn', 'gtp', 'elaid', 'thr_L', 'bilirub', 'dca', 'ch4s', 'f6p', 'creat', 'arg_L',
                                                        'cgly', 'lys_L', 'pro_L', 'pcholste_hs'],
                                         'Min': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
                                                 3.4e-06, 6.4e-06, 8e-06, 2.8e-05, 2.9e-05, 3e-05, 3.4e-05, 3.8e-05, 4e-05, 4.6e-05, 5e-05, 5.6e-05, 5.9e-05, 6.78e-05,
                                                 7e-05, 7.6e-05, 7.7e-05, 8e-05, 8e-05, 8e-05, 9e-05, 0.0001, 0.000167, 0.000168, 0.00018, 0.0002, 0.00021, 0.000226,
                                                 0.000248, 0.00029, 0.0003, 0.00031, 0.000352, 0.00041, 0.00053, 0.0006, 0.00076, 0.00076, 0.0008, 0.0008, 0.0008,
                                                 0.0008, 0.0008, 0.000908, 0.00094, 0.001, 0.001, 0.0014, 0.00144, 0.001474, 0.0016, 0.0016, 0.0016, 0.002, 0.002,
                                                 0.002, 0.0029, 0.003, 0.003, 0.003, 0.003, 0.003, 0.0031, 0.0032, 0.0045, 0.0049, 0.005, 0.005, 0.005, 0.005, 0.0051,
                                                 0.0051, 0.00529, 0.0054, 0.0054, 0.006, 0.007, 0.007, 0.0077, 0.008, 0.0083, 0.00868, 0.01, 0.01, 0.01, 0.01, 0.01,
                                                 0.01, 0.01067, 0.011, 0.012, 0.0126, 0.014, 0.014, 0.015, 0.015, 0.016, 0.017, 0.019, 0.019, 0.02, 0.02, 0.02, 0.02,
                                                 0.021, 0.022, 0.0222, 0.025, 0.027, 0.03, 0.032, 0.034, 0.035, 0.038, 0.04, 0.04, 0.042, 0.043, 0.0474, 0.05, 0.05,
                                                 0.05, 0.05729, 0.06, 0.06, 0.065, 0.07, 0.07, 0.073, 0.08087, 0.088, 0.09, 0.09, 0.094, 0.1, 0.11, 0.11, 0.12, 0.14,
                                                 0.146, 0.149, 0.16, 0.17, 0.175, 0.2, 0.207, 0.21, 0.231, 0.234, 0.24, 0.25, 0.25, 0.252, 0.27, 0.29, 0.3, 0.3, 0.3,
                                                 0.3, 0.32, 0.35, 0.35, 0.37, 0.376, 0.39, 0.4, 0.42, 0.43, 0.5, 0.5, 0.55, 0.56, 0.588, 0.6, 0.6, 0.63, 0.64, 0.66,
                                                 0.7, 0.7, 0.8, 0.9, 0.93, 1.07, 1.1, 1.11, 1.18, 1.2, 1.454, 1.56, 1.6, 1.8, 1.8, 1.8, 1.9, 10, 10, 10.7, 101, 11, 11,
                                                 11.3, 11.7, 1100, 11100, 12, 12.1, 1254, 13, 13.3, 132, 132600, 14, 143, 15, 15, 15, 15500, 16, 16.4, 166, 17, 17.9,
                                                 1766.2, 178, 19, 191.7, 2, 2, 2.72, 2.8, 2.83, 20, 20, 20, 208.9, 22, 22.03, 23, 24, 25, 25.5, 26, 26.6, 27, 27.4,
                                                 274, 2780, 29.4, 3.1, 3.3, 3.38, 3.4, 3.4, 3.9, 30, 31.1, 33, 34, 3440, 35, 37.1, 4, 4, 4.29, 4.36, 4.36, 4.48, 42,
                                                 42.2, 45.7, 5, 5, 5.1, 5.3, 5.5, 6, 6.4, 62.5, 68.7, 8.1],
                                         'Max': [0.0052, 0.00592, 0.0137, 0.01508, 0.024, 0.038, 0.23, 127, 2130, 25.6, 33, 6.463, 5, 0.12, 0.012, 0.022, 2.4e-06,
                                                 7.9e-05, 0.00012, 0.000164, 0.000659, 0.00072, 0.001016, 0.001359, 0.001632, 0.002047, 0.00284, 0.00285, 0.003,
                                                 0.00411, 0.0042, 0.00605, 0.00792, 0.0097, 0.011, 0.011, 0.011, 0.0113, 0.014, 0.017, 0.01701, 0.02481, 0.02664,
                                                 0.034, 0.03448, 0.0705, 0.089, 0.1057, 0.12, 0.12, 0.1339, 0.18, 0.186, 0.2, 0.293, 0.323, 0.374, 0.46, 0.47, 0.48,
                                                 0.5, 0.56, 0.619, 0.63, 0.78, 0.83, 1.06, 1.09, 1.208, 1.26, 1.4, 1.53, 1.73, 1.8, 1.8, 1010, 10400, 11.5, 120, 120,
                                                 120.422, 13.77, 145, 157, 16.49, 17.86, 18.1, 180.7, 19, 192, 196, 2, 2.6, 2.83, 2.848, 20, 213.016, 22.2, 23.4, 23.9,
                                                 2310, 258, 264, 282, 3.09, 3.712, 3.92, 3.924, 320, 33.205, 34.3, 37, 38.8, 381, 39.06, 394, 4, 4.08, 4.2, 4.3, 4.65,
                                                 4080, 42.4, 438, 46.55, 464.1, 5, 5.8, 54.8, 6.27, 6548, 70, 706, 8, 80, 85, 88.3, 9.08, 91.7, 92.5, 0.00018,
                                                 0.001249, 0.5, 0.6, 1, 1.3, 10, 10, 113.6, 20, 21, 24, 28, 4, 4, 5, 5, 5, 5, 5, 5, 58, 58, 80, 9.6, 0.48, 19, 32,
                                                 625, 2.23, 2.734, 57.775, 0.044, 0.21, 67, 0.000312, 0.00483, 0.01, 0.0287, 0.042, 0.1079, 0.85, 17, 4.2, 42.7, 3,
                                                 0.017, 0.00017, 0.0043, 1.27, 4.444, 0.177, 11.23, 0.77, 0.00012, 0.1322, 0.002, 0.001836, 0.00016, 0.023, 5e-05,
                                                 0.00023, 0.00075, 0.001671, 0.000307, 0.000561, 0.00189, 0.00012, 0.00022, 0.000324, 0.00036, 0.0004, 0.00012,
                                                 0.00012, 0.000186, 2.272, 0.000443, 0.000312, 0.0025, 0.0007, 0.00047, 0.00032, 0.00188, 0.00109, 0.0015, 0.0022,
                                                 0.000548, 0.0286, 0.148, 0.025, 0.00081, 0.053, 0.0012, 0.00328, 0.0051, 0.1263, 0.0485, 0.05, 0.00164, 0.005, 0.012,
                                                 0.0099, 0.00416, 0.002586, 0.0029, 0.0116, 0.0096, 0.01, 0.06, 14, 3.14, 0.007, 0.044, 0.087, 0.019, 0.011, 0.0047,
                                                 0.0072, 0.0485, 0.0061, 0.105, 0.009, 0.015, 0.16, 0.0147, 0.02, 0.00917, 1.77, 0.0066, 0.139, 0.06, 0.85, 0.6814,
                                                 0.009, 0.0103, 0.01132, 0.07, 0.4508, 1.159, 5.15, 0.25, 0.02, 0.01133, 0.618, 0.127, 0.0234, 0.071, 7.24, 0.053,
                                                 0.209, 0.069, 0.037, 0.051, 0.355, 0.024, 0.088, 1.63, 2.4, 0.26, 0.078, 0.0804, 0.153, 0.07, 6.6, 0.158, 0.11, 0.2,
                                                 0.0652, 0.08, 0.41, 0.514, 0.111, 0.0526, 1.97, 2.63, 3.73, 0.11185, 0.79, 4.14, 5.5, 0.11, 144.25, 0.227, 29.4, 0.24,
                                                 0.25, 0.11, 0.7, 0.94, 0.23, 37, 3.5, 0.46, 1.536, 0.353, 1.12, 0.53, 0.259, 4.2, 0.267, 0.46, 0.367, 1.318, 2.3,
                                                 0.57, 12.78, 0.708, 1.62, 4.45, 0.7, 1.5, 1.9, 0.7, 309, 1.71, 24, 4, 0.564, 4.71, 5.2, 0.5, 5.84, 20, 5.74, 3.482,
                                                 0.82, 1.814, 1.4, 9.4, 11.7, 2.96, 6.8, 3.4, 9.1, 21.4, 41, 2.45, 8.2, 3.7, 26.64, 9.78, 93.3, 2.346, 5.03, 8, 16.2, 
                                                 25, 3.2, 24.72, 51.4, 31, 52.93, 229, 171, 23, 17, 40, 20880, 28460, 60, 161, 2690, 31.5, 19.3, 188, 155400, 40, 629,
                                                 55, 85, 31, 28480, 80, 97, 2786, 65, 85.768, 10383, 518.48, 48, 247.7, 632.4, 34, 3.96, 92, 5.63, 261, 144, 100, 580,
                                                 400, 470, 25.6, 36, 88, 155.9, 79, 936.7, 135, 76.7, 1405, 7780, 87, 12.7, 17.7, 52.3, 25.7, 15.4, 34.1, 388.3, 275.8,
                                                 145, 81, 5680, 205.8, 172, 6, 48, 5, 16.16, 6.44, 5.68, 70, 128.08, 335, 21, 17, 6.3, 14.5, 130.4, 230, 91.6, 480,
                                                 380, 88.9]})
    elif comp == 'csf':
        # Define metabolites concentration
        metConcentration = pd.DataFrame({'Metabolite': ['34dhpha', 'prostgf2', 'leuktrB4woh', 'leuktrB4wcooh', 'leuktrB4', 'C14825', 'C14826', 'btn', '12harachd',
                                                        'CE2049', 'CE2047', 'CE2211', 'CE4890', 'mhista', 'C09642', 'melatn', 'xylu_D', 'iodine', '12ppd_R',
                                                        'normete_L', '5hoxindoa', '34dhphe', 'Lkynr', 'nrpphr', 'kynate', 'pydx5p', 'hista', 'srtn', 'imp',
                                                        '5htrp', 'rib_D', 'xylt', 'cala', 'HC02172', 'ppa', 'mal_L', '2hb', 'cit', 'ddca', 'citr_L', 'pro_L',
                                                        'CE2028', 'urea', 'man', 'C02356', 'octa', 'asp_L', 'for', '3mop', 'hdcea', 'oaa', 'ala_L', 'C08261',
                                                        'acac', 'tmlys', 'subeac', 'thym', 'glu_L', 'tyr_L', 'hdca', 'asn_L', 'val_L', 'succ', 'gal', '3hmp',
                                                        'trp_L', 'his_L', 'rbt', 'fru', 'leu_L', 'phe_L', 'ocdca', 'icit', 'xan', 'acnam', 'acetone', 'hcys_L',
                                                        '3mob', 'cholp', 'but', 'homoval', 'pyr', 'nac', 'ile_L', 'orn', 'glyclt', 'akg', 'ttdca', '4mop',
                                                        'met_L', 'meoh', 'fe2', 'ocdcea', '5oxpro', 'glyc_R', 'ins', 'cytd', 'Lcystin', 'ala_B', 'ura', 'CE2176',
                                                        'cl', 'ade', 'amp', 'bilirub', 'adn', 'gudac', '3uib', 'duri', 'crtsl', 'adrnl', 'pcrn', 'C11695',
                                                        '7dhchsterol', 'prostge2', 'dopa', 'thmmp', 'CE1918', 'thmpp', 'CE2705', 'C13856', '5g2oxpt', 'estradiol',
                                                        'leuktrC4', 'C14771', 'tststerone', 'prgstrn', 'andrstndn', 'dhea', 'CE0955', 'C04805', 'prgnlone',
                                                        'mepi', 'anth', '3moxtyr', 'xol27oh', '35cgmp', 'im4ac', 'xol24oh', 'cyst_L', 'c4crn', 'ahcys', 'quln',
                                                        '3mlda', 'ind3ac', 'thm', 'fol', 'pcholn204_hs', 'gthox', 'cys_L', 'ach', 'arachd', 'acrn', 'pcholste_hs',
                                                        '5mthf', 'spmd', '4hphac', 'Nacasp', '4abut', '34dhoxpeg', '56dura', 'gly', 'HC00900', 'bgly', 'pac',
                                                        'crvnc', 'sprm', 'ametam', 'ptrc', 'amet', 'CE5643', 'C10164', 'CE2006', '2mcit', 'ribflv', 'agm',
                                                        'adpac', 'pnto_R', 'no2', '56dthm', 'hxan', 'galt', 'pser_L', 'gthrd', 'taur', 'crn', 'cgly', 'fuM',
                                                        'chol', 'no', 'g3pc', 'CE0737', 'gtp', 'atp', 'gdp', 'adp', 'glyc', 'abt', 'k', 'ascb_L', 'creat', 'crtn',
                                                        'etha', 'bhb', 'gln_L', 'sucsal', 'ser_D', 'acald', 'lys_L', '2obut', 'ethamp', 'chsterol', 'xylu_L',
                                                        'thr_L', 'ac', 'inost', 'urate', 'hco3', 'h2o', 'lnlc', 'tcynt', 'alltn', 'lac_L', 'ser_L', 'arg_L', 
                                                        'ca2', 'pppi', 'glc_D', 'arab_L', 'sbt_D', 'na1', 'nh4'],
                                         'Min': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 5e-07, 8.314e-06, 1.5e-05, 2e-05, 2e-05, 2e-05, 2e-05, 6.996e-05, 7e-05, 8.2e-05, 0.00032,
                                                 0.00046, 0.001, 0.0011, 0.0014, 0.0023, 0.0033, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02,
                                                 0.02, 0.03, 0.03, 0.031, 0.04, 0.04, 0.04, 0.05, 0.05, 0.07, 0.1, 0.1, 0.1, 0.11, 0.11, 0.12, 0.12, 0.13, 0.13,
                                                 0.14, 0.15, 0.21, 0.24, 0.24, 0.25, 0.27, 0.32, 0.46, 0.49, 0.5, 0.6, 0.77, 0.8, 0.88, 0.9, 0.95, 0.99, 1, 1,
                                                 1.03, 1.5, 1.7, 1.77, 1.79, 1.81, 1.85, 11, 11.9, 1615, 17, 18, 19, 2.26, 207, 24, 24, 25, 27, 3, 3.4, 3.8, 3.9,
                                                 4, 4, 4, 4, 4.2, 440, 44000000, 5, 5.2, 5.6, 50, 6, 6.5, 7, 7, 740, 8, 8, 89692, 9.52],
                                         'Max': [0.006, 0.004, 1.5e-05, 1.5e-05, 0.00459, 0.00393, 0.00357, 0.001171, 0.00153, 0.0015, 0.00137, 0.00038, 0.01,
                                                 0.01, 0.01, 0.06, 5, 32, 133, 0.01, 0.8, 0.01, 0.07, 0.03, 0.08, 0.1, 0.38, 0.01, 19.7, 0.01, 5, 5, 0.5, 9.3, 9.2,
                                                 9, 88, 830, 8.4, 8.27, 8, 8, 7640, 72, 7.9, 7.5, 7.45, 64, 6, 6, 57, 51.4, 51, 50.8, 5.5, 5.1, 5, 46.4, 42.3, 42,
                                                 41.2, 38.3, 370, 364, 36, 36.8, 35.7, 35.6, 342, 34, 32.4, 30, 28, 27, 26.1, 247, 24.1, 21.8, 2.84, 2.8, 2.22,
                                                 194, 18.42, 18.2, 17.7, 154, 15.7, 15, 13.8, 13, 116, 11.57, 108, 107, 102, 1.34, 1.08, 0.83, 0.6, 0.5, 0.3, 0.22,
                                                 0.2, 0.2, 0.2, 0.14, 0.11, 0.11, 0.1, 0.1, 0.06, 0.04, 0.03, 0.02, 0.02, 0.01, 0.01, 0.44, 0.01, 0.01, 0.38, 0.25,
                                                 1.1e-05, 0.00020208, 0.000255, 0.002, 0.002, 0.002, 0.002, 0.00223, 0.00071, 0.00019, 0.00048, 0.00154, 0.004524,
                                                 0.0028, 0.0034, 0.0031, 0.0036, 0.04, 0.04, 0.04, 0.1, 0.12, 0.49, 0.25, 0.1, 0.038, 0.06, 0.11, 0.05, 0.66, 0.62,
                                                 0.107, 0.2, 0.2, 0.06, 1.57, 0.45, 0.34, 4.1, 19, 0.77, 0.12, 0.31, 0.18, 0.16, 0.28, 0.33, 0.38, 0.19, 0.46, 1.7,
                                                 1.07, 0.45, 0.45, 2.16, 1.14, 18.51, 1.7, 6.68, 2.53, 7.6, 1.52, 18.8, 9, 3.62, 57, 13, 10.8, 6.3, 8.3, 1.89,
                                                 1.91, 1.91, 1.95, 67, 45.9, 3300, 251, 70, 92.4, 20.1, 365, 1050, 33, 105, 51.3, 55, 6.2, 9.6, 8.76, 6, 70.5,
                                                 536, 250.94, 46.7, 12000, 66000000, 17, 50, 13, 4150, 78, 34.4, 2087, 8.76, 8690, 38, 30.73, 153000, 14.28]})
    elif comp == 'urine':
        # Define metabolites concentration
        metConcentration = pd.DataFrame({'Metabolite': ['C05769', 'xylu_L', 'iodine', 'HC00460', 'glcn', 'C05767', 'C05770', 'C05302', 'trypta', 'ppbng',
                                                        '13dampp', 'mhista', 'glucys', 'tym', '2hyoxplac', 'mma', 'leuktrB4', 'estriol', 'estradiol',
                                                        'tdchola', 'argsuc', 'dgchol', 'ddca', 'btn', 'ttdca', '5mta', 'normete_L', 'fe2', 'thm', 'ribflv',
                                                        'urcan', 'hom_L', 'octa', 'gchola', 'xtsn', '2hb', '56dthm', 'CE4968', 'Lcystin', 'urate', 'glx', 
                                                        'rbt', 'quln', 'mal_L', 'gly', 'ocdcea', '4pyrdx', 'cys_L', 'inost', 'rib_D', 'glyald', 'ac', 'gudac',
                                                        'thymd', '4hphac', '4abut', '5hoxindoa', 'ppa', '4hpro_LT', 'leu_L', 'carn', 'for', '2oxoadp', 'dgsn',
                                                        'din', 'duri', 'thym', 'pcollg5hlys', 'forglu', 'lnlc', 'hgentis', 'cyst_L', 'bhb', 'gal', '34hpp',
                                                        'malt', 'glu_L', 'slfcys', 'orn', 'subeac', 'arg_L', 'ura', 'hxan', 'orot', '56dura', '3mox4hoxm',
                                                        'adpac', 'taur', 'ala_L', 'asp_L', 'icit', 'val_L', 'pnto_R', 'fru', 'glyc', 'trp_L', 'bgly', 'gln_L',
                                                        'CE4970', 'HC00319', 'met_L', 'sbt_D', 'ind3ac', 'crn', 'cit', 'pro_L', 'crtn', 'ascb_L', 'sebacid',
                                                        'C02470', 'phpyr', 'cytd', 'C02356', '3mop', 'hcys_L', 'alltn', 'man', 'glyclt', 'lac_L', 'glcur',
                                                        'adn', 'homoval', 'sucr', 'arachd', 'asn_L', 'fum', '5oxpro', 'tyr_L', 'phe_L', 'ser_L', 'ins',
                                                        'acgam', 'Lkynr', 'Nacasp', 'creat', 'CE2028', 'acac', 'citr_L', 'HC00900', 'pyr', 'ala_B', 'but', 
                                                        'galt', 'akg', 'lys_L', 'thr_L', 'ile_L', 'sarcs', 'glyc_R', 'glyb', 'glc_D', 'fol', 'prostge2',
                                                        'thyox_L', 'adrnl', 'prostgf2', 'crtsl', 'bilirub', 'pmtcrn', 'dheas', 'nrpphr', '34dhphe', 'his_L',
                                                        'hista', 'srtn', 'gthrd', '3moxtyr', 'aldstrn', '5htrp', '7dhchsterol', 'etoh', 'pcrn', 'gsn', '5aop',
                                                        '4mop', 'uri', 'dad_2', 'ocdca', 'dopa', '2mcit', 'mev_R', 'glyleu', 'acald', 'gua', 'HC02191',
                                                        'camp', '34hpl', 'C08261', 'hdca', '3hpp', 'dcyt', 'andrstrn', 'ppi', 'anth', 'cgly', '34dhoxpeg', 
                                                        '1mncam', 'pac', 'glypro', 'acrn', 'tsul', 'oaa', 'hyptaur', 'ade', 'pser_L', 'c10crn', 'lcts',
                                                        'kynate', 'tststerone', 'dmgly', 'cala', 'chol', 'acglu', 'pydxn', 'arab_L', 'csn', 'xan', 'ethamp', 
                                                        'ca2', 'na1', '3hmp', 'acnam', 'L2aadp', 'etha', 'cl', 'k', 'fuc_L', 'urea', 'nh4', 'glcr', '4hbz',
                                                        '2obut', 'CE4969', 'nmn', 'sphgn', 'estrone', 'C05298', 'C05301', 'CE5072', 'andrstndn', 'eandrstrn',
                                                        'sphings', '17ahprgstrn', 'C05299', '11docrtsl', 'dhea', 'ahandrostan', 'N1sprm', 'c8crn', 'dchac',
                                                        'pcholhep_hs', 'pcholste_hs', 'pcholn204_hs'],
                                         'Min': [0, 0, 0, 0, 0, 0.00028, 0.000914, 0.0022, 0.01, 0.01, 0.015, 0.032, 0.18, 0.2, 0.41, 1.2, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.3e-05, 1.7e-05, 2e-05, 2e-05, 5.88e-05,
                                                 0.0019, 0.0019, 0.002, 0.0026, 0.00342, 0.0036, 0.0036, 0.0036, 0.004, 0.005, 0.0052, 0.006, 0.00624, 0.0068,
                                                 0.008, 0.01, 0.01, 0.01, 0.02, 0.03, 0.04, 0.04, 0.043, 0.05, 0.0515, 0.056, 0.065, 0.085, 0.1, 0.1, 0.1, 0.1,
                                                 0.1, 0.1, 0.11, 0.11, 0.12, 0.15, 0.18, 0.2, 0.24, 0.3, 0.32, 0.4, 0.41, 0.51, 0.52, 0.52, 0.53, 0.62, 0.63,
                                                 0.66, 0.69, 0.7, 0.704, 0.72, 0.8, 0.8, 0.8, 1.1, 1.3, 1.644, 16.9, 1863.49, 2, 2.5, 2.5, 4.91, 5263.2,
                                                 553.33, 6.38, 6.579, 724, 0, 0, 0.04, 0.5, 66.67, 0, 0, 0, 0, 0, 0, 0, 3.75e-05, 0.00078, 0.004, 0.0091,
                                                 0.016, 0.11, 0, 0.03, 0.05, 0.0008, 0.0013, 0.0005],
                                         'Max': [0.0066, 0.5, 1.65, 2.93, 38.8, 0.00156, 0.0092, 0.0058, 0.15, 0.6948, 1.44, 0.98, 0.74, 0.78, 54.98, 19.8,
                                                 5.8e-05, 0.0021, 0.00232, 0.00445, 0.0065, 0.0065, 0.07, 0.119, 0.12, 0.35, 0.3864, 0.4579, 0.6, 0.72, 0.78,
                                                 1.03, 1.13, 1.14, 1.31, 10, 10, 10, 100, 1023.89, 11, 11, 115.3, 121.6, 1234.76, 13, 13.49, 147.6, 148.2,
                                                 15.53, 150, 151.99, 159, 16.59, 1600.6, 169.51, 18, 18.3, 187.35, 19.07, 19.97, 195.63, 2, 2, 2, 2, 2, 2.26,
                                                 2.32, 2.41, 2.8, 20.7, 200, 2026.02, 21.3, 21.4, 217.61, 22, 22.1, 22.31, 22.97, 23.85, 24.1, 24.1, 24.6,
                                                 25.9, 250, 251, 252.97, 26.8, 265.3, 27.68, 28.92, 285.46, 288.5, 29.4, 2921.04, 298.84, 3.3, 3.5, 3.947,
                                                 30.46, 31, 31.5, 3168.4, 34.76, 34337, 349.4, 35.1, 4, 4, 4.18, 4.4, 4.8, 4.89, 40.7, 42.49, 425, 444.29,
                                                 453.51, 46.7, 47.3, 486.2, 5.71, 53.01, 53.7, 54.2, 55.31, 55.63, 57.2, 6.5, 6.66, 6.73, 65.5, 654, 66.4, 67,
                                                 74.25, 79.4, 79.53, 8.4, 81, 81, 836.3, 86, 86.96, 9.25, 9.4, 91.1, 92.7, 942.7, 0.012, 0.019, 0.00062,
                                                 0.0099, 0.00055, 0.021, 5.13, 0.6, 2.52, 0.057, 0.04, 212.81, 37.44, 0.12, 0.085, 0.72, 0.014, 0.00936,
                                                 0.051, 3.1, 0.2, 0.23, 4.4, 0.45, 1.351, 0.4, 7.7, 5.221, 28.6, 0.6, 0.084, 43.4, 0.57, 0.2, 0.66, 11.5, 
                                                 17.6, 26.1, 48, 1.19, 96.48, 5, 1, 1.54, 0.3, 15, 5.98, 0.48, 7.5, 3.77, 5.94, 0.78, 5.1, 10.4, 15.2, 848.41,
                                                 4.2, 1.45, 11.2, 1.056, 28.5, 1.2, 1.2, 22, 10.7, 4.6, 12.2, 641, 37249.07, 59.8, 11.76, 16.4, 56.2, 17763.2,
                                                 8078.58, 26, 38812.4, 4704, 14.97, 29, 5.4, 3.8, 600, 0.00063, 0.003, 0.0041, 0.0044, 0.064, 0.206, 0.83,
                                                 3.947e-05, 0.00182, 0.006, 178.9, 0.244, 0.67, 0.08, 0.09, 0.06, 0.0087, 0.01, 0.0038]})
        
    # Output variable
    return(metConcentration)

def _adjustWholeBodyRxnCoeff(model: cobra.Model, indParam: dict) -> (cobra.Model):
    # Define model metabolites
    modMets = [x.id for x in model.metabolites]
    
    # Get model objective reaction
    rxn = model.reactions.get_by_id('Whole_body_objective_rxn')
    
    # Define dummy objectives
    dummObjs = [x for x in modMets if '_dummy_objective' in x]
    
    # Adjust whole body biomass coefficient
    for org in indParam['organs']:
        # Find dummy reaction for organ
        organID = [x for x in modMets if (org + '_biomass') in x]
        organID = [x for x in organID if x in dummObjs]
        
        # Set stoichiometric coefficient
        if len(organID) > 0:
            met = model.metabolites.get_by_id(organID[0])
            rxn.add_metabolites({met: -indParam['organWeightFract'][org] * 100}, combine = False)
            
    # Output variable
    return (model)

def _physiologicalConstraints(model: cobra.Model, indParam: dict, bloodData: pd.DataFrame, organsListExt: list, Type: str = 'HMDB', inputData: pd.DataFrame = [], bioFluid: str = 'all') -> (cobra.Model):
    # Define model reactions
    modRxns = [x.id for x in model.reactions]
    
    # Define concentration limits
    minConcConstraint = 5
    maxConcConstraint = 50
    
    # Lower bound option
    setLB = 0
    
    ########## Calculate GFR = Glomerular Filtration Rate (20% of the renal plasma flow) ##########
    GlomerularFiltrationRate = indParam['GlomerularFiltrationRate']

    ########## Read metabolite concentration data ##########
    if Type == 'HMDB':
        metConcDataBc = _defaultMetabolitesConcentration('blood')
        metConcDataCSF = _defaultMetabolitesConcentration('csf')
        metConcDataUr = _defaultMetabolitesConcentration('urine')
    elif Type == 'direct':
        if bioFluid == 'bc':
            metConcDataBc = inputData
        elif bioFluid == 'csf':
            metConcDataCSF = inputData
        elif bioFluid == 'u':
            metConcDataUr = inputData
    
    ########## Cardiac output and organ-specific blood flow rate ##########
    bloodFlowRate = [None] * len(organsListExt)
    plasmaFlowRate = [None] * len(organsListExt)
    for ind, org in enumerate(organsListExt):
        index = bloodData.loc[:,'Organ'] == org
        if any(index):
            bloodFlowRate[ind] = bloodData.loc[index, 'Blood Flow Percentage'].iloc[0] * indParam['CardiacOutput']
            plasmaFlowRate[ind] = bloodData.loc[index, 'Blood Flow Percentage'].iloc[0] * indParam['CardiacOutput'] * (1 - indParam['Hematocrit'])
        elif org == 'BBB':
            BScord = bloodData.loc[bloodData.Organ == 'Scord', 'Blood Flow Percentage'].iloc[0]
            BBrain = bloodData.loc[bloodData.Organ == 'Brain', 'Blood Flow Percentage'].iloc[0]
    
            bloodFlowRate[ind] = (BScord + BBrain) * indParam['CardiacOutput']
            plasmaFlowRate[ind] = (BScord + BBrain) * indParam['CardiacOutput'] * (1 - indParam['Hematocrit'])
        else:
            bloodFlowRate[ind] = 0.01 * indParam['CardiacOutput']
            plasmaFlowRate[ind] = 0.01 * indParam['CardiacOutput'] * (1 - indParam['Hematocrit'])
    
    # Organs excluded from secretion into bc
    exclOrgan = ['sIEC', 'Colon', 'Spleen', 'Pancreas', 'Gall', 'Brain']
    
    ########## Compute maximal possible uptake and secretion rate for each metabolite (BLOOD) ##########
    if bioFluid == 'all' or bioFluid == 'bc':
        # Excluding reactions
        R_Kidney = ['Kidney_EX_na1_e__bc', 'Kidney_EX_hco3_e__bc', 'Kidney_EX_urea_e__bc', 'Kidney_EX_k_e__bc', 'Kidney_EX_cl_e__bc', 'Kidney_EX_ca2_e__bc',
                    'Kidney_EX_HC02172_e__bc', 'Kidney_EX_avite1_e__bc']
        R = ['BBB_NH4_CSFupt', 'BBB_CHOL_CSFupt', 'BBB_PI_CSFupt', 'BBB_STRDNC_CSFupt', 'BBB_HC00250_CSFupt', 'BBB_PYDXN_CSFupt', 'BBB_5MTHF_CSFupt', 'BBB_SO3_CSFupt']
            
        for org in organsListExt:
            # Find objective reactions
            if org == 'BBB':
                ExR = [x for x in modRxns if x.endswith('CSFupt') and x not in R and not x.endswith('H2O_CSFupt')]
            else:
                ExR = [x for x in modRxns if x.startswith(org + '_EX_') and x.endswith('_bc') and '_o2s_e' not in x and '_h2o_e' not in x]
               
            # Filter objetive reactions if input data is provided
            if Type == 'direct':
                ExR = [rxn for rxn in ExR if any(pat in rxn for pat in metConcDataBc.loc[:,'Metabolite'].tolist())]
                
            for reaction in ExR:
                rxn = model.reactions.get_by_id(reaction)
                
                # Get metabolite ID associated with the reaction
                ExM = [met.id for met, coeff in rxn.metabolites.items() if coeff > 0][0]
                mName = re.sub(r'_bc$', '', ExM)
                if org != 'Kidney':
                    mName = re.sub(r'^RBC_', '', mName)
                
                if (metConcDataBc['Metabolite'] == mName).any():
                    MConMin = metConcDataBc.loc[metConcDataBc['Metabolite'] == mName,'Min'].iloc[0]
                    MConMax = metConcDataBc.loc[metConcDataBc['Metabolite'] == mName,'Max'].iloc[0]
                else:
                    MConMin = 0
                    MConMax = indParam['MConDefaultBc']
                    
                # Kidney is treated differently
                if org == 'Kidney':
                    MSecretRateKidneyLB = 1 * (MConMax/1000) * GlomerularFiltrationRate * 60 * 24 / 1000
                    MSecretRateKidneyUB = 1 * (MConMin/1000) * GlomerularFiltrationRate * 60 * 24 / 1000
            
                    if reaction not in R_Kidney:
                        # Temporaly change bounds (avoid errors)
                        rxn.lower_bound = -1e6
                        
                        # Set upper bound
                        if setLB == 1:
                            rxn.upper_bound = -MSecretRateKidneyUB
                        else:
                            if MConMax >= maxConcConstraint:
                                rxn.upper_bound = -MSecretRateKidneyUB
                            else:
                                rxn.upper_bound = 0
                    else:
                        rxn.upper_bound = 0
                        
                    # Set lower bound
                    rxn.lower_bound = -MSecretRateKidneyLB
                    
                else:
                    MuptakeRateBc = (MConMax/1000) * plasmaFlowRate[organsListExt.index(org)] * 60 * 24 / 1000
                    
                    # Set lower bound
                    if rxn.lower_bound < 0:
                        if MConMax > 1e-3:
                            rxn.lower_bound = -1*MuptakeRateBc
                        else:
                            MuptakeRateBc = ((1e-3)/1000) * plasmaFlowRate[organsListExt.index(org)] * 60 * 24 / 1000
                            rxn.lower_bound = -1*MuptakeRateBc
    
    ########## Compute maximal possible uptake and secretion rate for each metabolite (BRAIN) ##########
    if bioFluid == 'all' or bioFluid == 'csf':
        # Define excluding metabolites
        R = ['na1_csf', 'cl_csf', 'k_csf', 'h2o_csf', 'sucsal_csf', 'ca2_csf', 'ser_D_csf']
        
        # Define objective reactions
        ExR = [x for x in modRxns if 'BBB' in x and 'CSF' in x and 'exp' in x and '_o2_e' not in x and '_o2s_e' not in x and '_co2_e' not in x]
        
        # Filter objetive reactions if input data is provided
        if Type == 'direct':
            ExR = [rxn for rxn in ExR if any(pat in rxn for pat in metConcDataCSF.loc[:,'Metabolite'].tolist())]
        
        for reaction in ExR:
            rxn = model.reactions.get_by_id(reaction)
            
            # Get metabolite ID associated with the reaction
            ExM = [met.id for met, coeff in rxn.metabolites.items() if coeff < 0][0]
            mName = re.sub(r'_csf', '', ExM)
            
            # Get maximal concentration of the metabolite
            if any(metConcDataCSF['Metabolite'] == mName):
                MConMin = metConcDataCSF.loc[metConcDataCSF['Metabolite'] == mName,'Min'].iloc[0]
                MConMax = metConcDataCSF.loc[metConcDataCSF['Metabolite'] == mName,'Max'].iloc[0]
            else:
                MConMin = 0
                MConMax = indParam['MConDefaultCSF']
                
            # Set brain to blood secretion
            MSecretRateCSFLB = (MConMin/1000) * indParam['CSFBloodFlowRate'] * 60 * 24 / 1000
            MSecretRateCSFUB = (MConMin/1000) * indParam['CSFBloodFlowRate'] * 60 * 24 / 1000
            
            # Temporaly change bounds (avoid errors)
            rxn.upper_bound = 1e6
            
            if setLB == 1:
                rxn.lower_bound = MSecretRateCSFLB
            else:
                if MConMin >= minConcConstraint and MConMax >= maxConcConstraint and ExM not in R:
                    rxn.lower_bound = MSecretRateCSFLB
                else:
                    rxn.lower_bound = 0
            rxn.upper_bound = MSecretRateCSFUB
    
    ########## Compute maximal possible uptake and secretion rate for each metabolite (URINE) ##########
    if bioFluid == 'all' or bioFluid == 'u':
        # Convert creatinine from mg/dL into mmol/L
        MConDefaultUrCreatinineMin = indParam['MConUrCreatinineMin'] * 10 / 113.1179
        MConDefaultUrCreatinineMax = indParam['MConUrCreatinineMax'] * 10 / 113.1179
        
        # Define excluding reactions
        R = ['EX_na1_u', 'EX_cl_u', 'EX_k_u', 'EX_ca2_u', 'EX_C05767_u', 'EX_C05770_u', 'EX_C05302_u', 'EX_trypta_u', 'EX_ppbng_u', 'EX_13dampp_u', 'EX_mhista_u', 'EX_tym_u',
             'EX_2hyoxplac_u', 'EX_pmtcrn_u', 'EX_dheas_u', 'EX_34dhphe_u', 'EX_srtn_u', 'EX_gthrd_u', 'EX_pcholhep_hs_u', 'EX_pcholste_hs_u', 'EX_pcholn204_hs_u', 'EX_3moxtyr_u',
             'EX_aldstrn_u', 'EX_tststerone_u', 'EX_pydxn_u', 'EX_sphgn_u', 'EX_sphings_u', 'EX_csn_u', 'EX_arab_L_u', 'EX_tststerone_u', 'EX_tststerone_u', 'EX_pydxn_u', 'EX_mma_u',
             'EX_tsul_u', 'EX_5htrp_u', 'EX_7dhchsterol_u', 'EX_etoh_u', 'EX_gsn_u', 'EX_5aop_u', 'EX_uri_u', 'EX_dad_2_u', 'EX_ocdca_u', 'EX_gua_u', 'EX_dcyt_u', 'EX_glyleu_u',
             'EX_acald_u', 'EX_HC02191_u']
        
        # Define objective reactions
        ExR = [x for x in modRxns if x.startswith('EX_') and x.endswith('_u') and '_o2_e' not in x and '_o2s_e' not in x and '_co2_e' not in x]
        
        # Filter objetive reactions if input data is provided
        if Type == 'direct':
            ExR = [rxn for rxn in ExR if any(pat in rxn for pat in metConcDataUr.loc[:,'Metabolite'].tolist())]
        
        for reaction in ExR:
            rxn = model.reactions.get_by_id(reaction)
            
            # Get metabolite ID associated with the reaction
            ExM = [met.id for met, coeff in rxn.metabolites.items() if coeff < 0][0]
            mName = re.sub(r'_u$','', ExM)
            
            # Get maximal concentration of the metabolite
            if any(metConcDataUr['Metabolite'] == mName):
                MConMin = metConcDataUr.loc[metConcDataUr['Metabolite'] == mName,'Min'].iloc[0]
                MConMax = metConcDataUr.loc[metConcDataUr['Metabolite'] == mName,'Max'].iloc[0]
            else:
                MConMin = indParam['MConDefaultUrMin']
                MConMax = indParam['MConDefaultUrMax']
                
            # Set urine secretion
            MSecrRateUrLB = (MConMin/1000) * MConDefaultUrCreatinineMin * indParam['UrFlowRate'] / 1000
            MSecrRateUrUB = (MConMax/1000) * MConDefaultUrCreatinineMax * indParam['UrFlowRate'] / 1000
            
            if rxn.upper_bound > 0:
                # Temporaly change bounds (avoid errors)
                rxn.upper_bound = 1e6
                
                if setLB == 1:
                    if reaction not in R:
                        rxn.lower_bound = MSecrRateUrLB
                    else:
                        reaction.lower_bound = 0
                else:
                    if MConMax >=  maxConcConstraint and MConMin >= minConcConstraint and reaction not in R:
                        rxn.lower_bound = MSecrRateUrLB
                    else:
                        rxn.lower_bound = 0
                rxn.upper_bound = MSecrRateUrUB
            
    ########## Woman is not producing milk ##########
    tmpRxns = [x for x in modRxns if 'miB__mi' in x]
    for r in tmpRxns:
        rxn = model.reactions.get_by_id(r)
        rxn.lower_bound = 0
        rxn.upper_bound = 0
    
    if Type != 'direct':
        ########## Set o2_a and co2_a constraints ##########
        rxn = model.reactions.get_by_id('EX_o2_a')
        rxn.lower_bound = -1e6
        rxn.upper_bound = 1e6
        
        rxn.lower_bound = -25000
        rxn.upper_bound = -15000
        
        rxn = model.reactions.get_by_id('EX_co2_a')
        rxn.lower_bound = -1e6
        rxn.upper_bound = 1e6
        
        rxn.lower_bound = 15000 * 0.8
        rxn.upper_bound = 25000
        
        ########## Set water constraints ##########
        # Air
        rxn = model.reactions.get_by_id('EX_h2o_a')
        rxn.lower_bound = -1e6
        rxn.upper_bound = 1e6
        
        rxn.lower_bound = 47182 * 0.8
        rxn.upper_bound = 47182 * 1.2
        
        # Sweating
        rxn = model.reactions.get_by_id('EX_h2o_sw')
        rxn.lower_bound = -1e6
        rxn.upper_bound = 1e6
        
        rxn.lower_bound = 36080 * 0.8
        rxn.upper_bound = 36080 * 1.2
        
        # Urine
        rxn = model.reactions.get_by_id('EX_h2o_u')
        rxn.lower_bound = -1e6
        rxn.upper_bound = 1e6
        
        rxn.lower_bound = 77711 * 0.8
        rxn.upper_bound = 77711 * 1.2
        
        # Feces
        rxn = model.reactions.get_by_id('Excretion_EX_h2o_fe')
        rxn.lower_bound = -1e6
        rxn.upper_bound = 1e6
        
        rxn.lower_bound = 5550 * 0.8
        rxn.upper_bound = 5550 * 1.2
        
        # Limit water secretion into bile duct
        rxn = model.reactions.get_by_id('Gall_H2Ot_bdG')
        rxn.upper_bound = 1000
        
        rxn = model.reactions.get_by_id('Liver_H2Ot_bdL')
        rxn.upper_bound = 1000
        
        ########## Muscle can only take up glucose ##########
        rxn = model.reactions.get_by_id('Muscle_EX_glc_D_e__bc')
        if rxn.upper_bound >= -0.01 * 1000 and rxn.lower_bound <= -0.01*1000:
            rxn.upper_bound = -0.01 * 1000
        elif rxn.lower_bound > rxn.upper_bound:
            rxn.upper_bound = 0
        
        ########## Muscle_EX_ala_l(e)_[bc] Muscle_ala_L[e]  <=> ala_L[bc] 12.5 mg alanine/min/person (65 kg) ##########
        # Calculate mmol per day per person
        met = (12.5 * 60 * 24 * (indParam['bodyWeight'] / 65)) / 89.09
        
        # Set constraint
        rxn = model.reactions.get_by_id('Muscle_EX_ala_L_e__bc')
        if rxn.lower_bound < met * 0.8 and rxn.upper_bound >= met * 0.8:
            rxn.lower_bound = met * 0.8
            rxn.upper_bound = met * 1.2
        elif rxn.lower_bound > rxn.upper_bound:
            rxn.lower_bound = 0
        
        ########## RBC_EX_glc(e)_[bc]	RBC_glc_D[e]  <=> glc_D[bc]	RBC	 25 mg glc/min/person (65kg) ##########
        # Calculate mmol per day per person
        met = (25 * 60 * 24 * (indParam['bodyWeight'] / 65)) / 180.16
        
        # Set constraint
        rxn = model.reactions.get_by_id('RBC_EX_glc_D_e__bc')
        if rxn.lower_bound < -met * 1.2 and rxn.upper_bound >= -met * 1.2:
            rxn.lower_bound = -met * 1.2
            rxn.upper_bound = -met * 0.8
        elif rxn.lower_bound > rxn.upper_bound:
            rxn.upper_bound = 0
        
        ########## Brain_EX_glc(e)_[csf]	Brain_glc_D[e]  <=> glc_D[csf] brain	80 mg glc/min/person ##########
        # Calculate mmol per day per person
        met = (80 * 60 * 24 * (indParam['bodyWeight'] / 65)) / 180.16
        
        # Set constraint
        rxn = model.reactions.get_by_id('Brain_EX_glc_D_e__csf')
        if rxn.lower_bound < -met * 1.2 and rxn.upper_bound >= -met * 1.2:
            rxn.lower_bound = -met * 1.2
            rxn.upper_bound = -met * 0.8
        elif rxn.lower_bound > rxn.upper_bound:
            rxn.upper_bound = 0
        
        ######### Brain water exchange ##########
        # Calculate oxygen amount
        brain_o2 = (156 * 60 * 24 * indParam['organWeight']['Brain'] / 100) / 1000
        
        # Set constraint
        rxn = model.reactions.get_by_id('Brain_EX_o2_e__csf')
        if rxn.lower_bound < -(brain_o2*1.2) and rxn.upper_bound >= -(brain_o2*1.2):
            rxn.lower_bound = -brain_o2 * 1.2
            rxn.upper_bound = -brain_o2 * 0.7
        elif rxn.lower_bound > rxn.upper_bound:
            rxn.upper_bound = 0
        
        ######### Liver_EX_ala_l(e)_[bc] Liver_ala_L[e]  <=> ala_L[bc] liver	12.5 mg alanine/min/person (65 kg) ##########
        # Calculate mmol per day per person
        met = (12.5 * 60 * 24 * (indParam['bodyWeight'] / 65)) / 89.09
        
        # Set constraint
        rxn = model.reactions.get_by_id('Liver_EX_ala_L_e__bc')
        if rxn.lower_bound < -met * 1.2 and rxn.upper_bound >= met * 1.2:
            rxn.lower_bound = -met * 1.2
            rxn.upper_bound = -met * 0.8
        elif rxn.lower_bound > rxn.upper_bound:
            rxn.upper_bound = 0
        
        ########### Liver_EX_glc(e)_[bc] Liver_glc_D[e]  <=> glc_D[bc] liver	130 mg glc/min/person (65 kg) ###########
        # Calculate mmol per day per person
        met = (130 * 60 * 24 * (indParam['bodyWeight'] / 65)) / 180.16
        
        # Set constraint
        rxn = model.reactions.get_by_id('Liver_EX_glc_D_e__bc')
        if rxn.lower_bound < met * 0.8 and rxn.upper_bound >= met * 0.8:
            rxn.lower_bound = met * 0.8
            rxn.upper_bound = met * 1.2
        elif rxn.lower_bound > rxn.upper_bound:
            rxn.lower_bound = 0
            
        ########## Adipocytes_EX_glyc(e)_[bc] Adipocytes_glyc[e]  <=> glyc[bc] adipocytes	12 mg glycerol/min/person (65 kg) ##########    
        # Calculate mmol per day per person
        met = (12 * 60 * 24 * (indParam['bodyWeight'] / 65)) / 92.09
        
        # Set constraint
        rxn = model.reactions.get_by_id('Adipocytes_EX_glyc_e__bc')
        if rxn.lower_bound < met * 0.8 and rxn.upper_bound >= met * 0.8:
            rxn.upper_bound = met * 1.2
        elif rxn.lower_bound > rxn.upper_bound:
            rxn.lower_bound = 0
        
        ########## Brain and liver can do co2 fixation ##########    
        rxns = [x for x in modRxns if x.endswith('EX_co2_e__bc')]
        for r in rxns:
            rxn = model.reactions.get_by_id(r)
            rxn.lower_bound = 0
        
        ########## Co2 can cross BBB ##########    
        rxns = ['Brain_EX_co2_e__csf', 'Liver_EX_co2_e__bc', 'Lung_EX_co2_e__bc', 'Kidney_EX_co2_e__bc']
        for r in rxns:
            rxn = model.reactions.get_by_id(r)
            rxn.lower_bound = -10000
        
        ########## Brain Energy Requirement ##########    
        rxn = model.reactions.get_by_id('Brain_DM_atp_c_')
        rxn.lower_bound = 0
        
        ########## Heart Energy Requirement ##########
        rxn = model.reactions.get_by_id('Heart_DM_atp_c_')
        rxn.lower_bound = 6000
    
    # Output variable
    return (model)

def _persWBMmetabolomics(sample: str, model: cobra.Model, meta: pd.DataFrame, persMets: pd.DataFrame, indParam: dict, bloodData: pd.DataFrame, organsListExt: list) -> (cobra.Model):
    # Load VMH database
    dataPath = pkg_resources.files('persephone').joinpath('data/mets/')
    VMH = pd.read_csv(dataPath.joinpath('VMH_metabolites.csv'))
    
    # Check all metabolites present in VMH
    idx = persMets.iloc[:,0].isin(VMH.loc[:,'metID'])
    if any(idx == False):
        raise ValueError(f'The following metabolites were not present in VMH database: {persMets.loc[~idx,persMets.columns[0]].tolist()}. Please check metabolite ID matched name on vmh.life')
    
    # Check metabolites compartments
    idx = persMets.iloc[:,1].isin(['bc','csf','u'])
    if any(idx == False):
        raise ValueError('Some biofluids provided are incompatible. Available biofluids are: [bc, csf, u].')
    
    # Initialize biofluid data
    bioFluids = []
    metBcData = pd.DataFrame(columns=['Metabolite', 'Min', 'Max'])
    metCsfData = pd.DataFrame(columns=['Metabolite', 'Min', 'Max'])
    metUrData = pd.DataFrame(columns=['Metabolite', 'Min', 'Max'])
    
    # Build metabolite data for all the biofluids
    for i in range(len(persMets)):
        # Get molecular weight from VMH
        MW = VMH.loc[VMH.loc[:,'metID'] == persMets.iloc[i,0],'MW'].iloc[0]
        
        # Get the concentration and the unit of the metabolite
        Conc = meta.loc[sample, persMets.iloc[i,0]]
        unit = persMets.iloc[i,2]
        
        # Unit conversion to micro mol per litre
        if unit in ['µmol/L', 'umol/L', 'uM']:
            Conc = Conc
        elif unit == 'mg/dL':
            Conc = Conc * (10**4) / MW
        elif unit == 'g/dL':
            Conc = Conc * (10**7) / MW
        elif unit == 'pg/mL':
            Conc = Conc / MW
        elif unit == 'mmol/L':
            Conc = Conc * (10**3)
        elif unit == 'ng/dL':
            Conc = Conc * (10**2) / MW
        else:
            raise ValueError(f'Unrecognized unit: {unit}. Available units for metabolite concentrations are: [µmol/L, umol/L, uM, mg/dL, g/dL, pg/mL, mmol/L, ng/dL].')
            
        # Add information to specific biofluid
        row = {'Metabolite': persMets.iloc[i,0], 'Min': Conc * 0.8, 'Max': Conc * 1.2}
        if persMets.iloc[i,1] == 'bc':
            metBcData = pd.concat([metBcData, pd.DataFrame([row])], ignore_index=True)
            bioFluids.append('bc')
        elif persMets.iloc[i,1] == 'csf':
            metCsfData = pd.concat([metCsfData, pd.DataFrame([row])], ignore_index=True)
            bioFluids.append('csf')
        elif persMets.iloc[i,1] == 'u':
            metUrData = pd.concat([metBcData, pd.DataFrame([row])], ignore_index=True)
            bioFluids.append('u')
    
    # Set metabolic constraints
    pWBM = model
    if 'bc' in bioFluids:
        pWBM = _physiologicalConstraints(pWBM, indParam, bloodData, organsListExt, 'direct', metBcData, 'bc')
    elif 'csf' in bioFluids:
        pWBM = _physiologicalConstraints(pWBM, indParam, bloodData, organsListExt, 'direct', metCsfData, 'csf')
    elif 'u' in bioFluids:
        pWBM = _physiologicalConstraints(pWBM, indParam, bloodData, organsListExt, 'direct', metUrData, 'u')
    
    # Output variable
    return (pWBM)
    
def _generatePersWBM(sample: str, meta: pd.DataFrame, personalisingPhys: bool, persPhysiology: pd.DataFrame, personalisingMets: bool, persMets: pd.DataFrame, pers_wbm_path: str):
    if (not os.path.exists(os.path.join(pers_wbm_path, 'models', sample +'.pkl'))):
        # Load original WBM
        if meta.loc[sample,'sex'] == 'male':
            metInfo = pd.read_csv(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvey_1_03d.metInfo.csv'))
            rxnInfo = pd.read_csv(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvey_1_03d.rxnInfo.csv'))
            with open(str(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvey_1_03d.pkl')), "rb") as f:
                wbm_model = pickle.load(f)
        elif meta.loc[sample,'sex'] == 'female':
            metInfo = pd.read_csv(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvetta_1_03d.metInfo.csv'))
            rxnInfo = pd.read_csv(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvetta_1_03d.rxnInfo.csv'))
            with open(str(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvetta_1_03d.pkl')), "rb") as f:
                wbm_model = pickle.load(f)
        
        # Change model ID
        wbm_model.id = 'iWBM_' + sample
        
        # Initialize parameter log file
        paramSource = {}
        paramSource['sex'] = 'User defined'
        paramSource['bloodVolume'] = 'Estimated based on sex, height and weight'
        
        # Personalising physiological parameters
        if personalisingPhys:
            # Load default parameters
            indParam, bloodData, organsListExt = _defaultIndividualParameters(meta.loc[sample, 'sex'])
            
            # BODY WEIGHT
            if any(persPhysiology.iloc[:,0] == 'body weight'):
                idxWt = meta.columns == 'body weight'
                if sum(idxWt) > 1:
                    raise ValueError('More than one column in the metadata contains the parameter body weight!')
                    
                # Extract value and unit
                wt = meta.loc[sample,'body weight']
                currentUnit = persPhysiology.loc[persPhysiology.iloc[:, 0] == 'body weight', persPhysiology.columns[1]].iloc[0]
                
                # Unit transformation
                if currentUnit == 'kg':
                    wt = wt
                elif currentUnit == 'lb':
                    wt = wt * 0.453592
                else:
                    raise ValueError('Parameter body weight not provided in a valid unit. Available units are [kg, lb]')
                
                # Save parameter
                indParam['bodyWeight'] = wt
                paramSource['bodyWeight'] = 'User defined'
                
            else:
                paramSource['bodyWeight'] = 'Default'
                
            # AGE
            if any(persPhysiology.iloc[:,0] == 'age'):
                idxAge = meta.columns == 'age'
                if sum(idxAge) > 1:
                    raise ValueError('More than one column in the metadata contains the parameter age!')
                    
                # Extract value and unit
                age = meta.loc[sample,'age']
                currentUnit = persPhysiology.loc[persPhysiology.iloc[:, 0] == 'age', persPhysiology.columns[1]].iloc[0]
                
                # Unit transformation
                if currentUnit == 'days':
                    age = age / 365
                elif currentUnit == 'months':
                    age = age / 12
                elif currentUnit == 'years':
                    age = age
                else:
                    raise ValueError('Parameter age not provided in a valid unit. Available units are [days, months, years]')
                
                # Save parameter
                indParam['age'] = age
                paramSource['age'] = 'User defined'
                
            else:
                paramSource['age'] = 'Default'
                
            # HEIGHT
            if any(persPhysiology.iloc[:,0] == 'height'):
                idxHt = meta.columns == 'height'
                if sum(idxHt) > 1:
                    raise ValueError('More than one column in the metadata contains the parameter age!')
                    
                # Extract value and unit
                height = meta.loc[sample,'height']
                currentUnit = persPhysiology.loc[persPhysiology.iloc[:, 0] == 'height', persPhysiology.columns[1]].iloc[0]
                
                # Unit transformation
                if currentUnit == 'cm':
                    height = height
                elif currentUnit == 'ft':
                    height = height * 30.48
                elif currentUnit == 'in':
                    height = height * 2.54
                else:
                    raise ValueError('Parameter height not provided in a valid unit. Available units are [days, months, years]')
                
                # Save parameter
                indParam['Height'] = height
                paramSource['height'] = 'User defined'
                
            else:
                paramSource['height'] = 'Default'
            
            # HEART RATE
            if any(persPhysiology.iloc[:,0] == 'heart rate'):
                idxHr = meta.columns == 'heart rate'
                if sum(idxHr) > 1:
                    raise ValueError('More than one column in the metadata contains the parameter heart rate!')
                
                # Extract value and unit
                Hr = meta.loc[sample, 'heart rate']
                currentUnit = persPhysiology.loc[persPhysiology.iloc[:, 0] == 'heart rate', persPhysiology.columns[1]].iloc[0]
        
                # Unit transformation
                if currentUnit == 'bpm':
                    indParam['HeartRate'] = Hr
                    paramSource['heartRate'] = 'User defined'
                else:
                    raise ValueError('Heart rate must be provided in bpm')
                    
            else:
                paramSource['heartRate'] = 'Default'
            
            # STROKE VOLUME
            if any(persPhysiology.iloc[:,0] == 'stroke volume'):
                idxSV = meta.columns == 'stroke volume'
                if sum(idxSV) > 1:
                    raise ValueError('More than one column in the metadata contains the parameter stroke volume!')
                
                # Extract value and unit
                SV = meta.loc[sample, 'stroke volume']
                currentUnit = persPhysiology.loc[persPhysiology.iloc[:, 0] == 'stroke volume', persPhysiology.columns[1]].iloc[0]
        
                # Unit transformation
                if currentUnit == 'mL':
                    indParam['StrokeVolume'] = SV
                    paramSource['strokeVolume'] = 'User defined'
                else:
                    raise ValueError('Stroke volume must be provided in mL')
            
            else:
                paramSource['strokeVolume'] = 'Default'
            
            # HEMATOCRIT
            if any(persPhysiology.iloc[:,0] == 'hematocrit'):
                idxHmt = meta.columns == 'hematocrit'
                if sum(idxHmt) > 1:
                    raise ValueError('More than one column in the metadata contains the parameter hematocrit!')
                
                # Extract value
                Hmt = meta.loc[sample, 'hematocrit']
                
                # Unit transformation (proportion)
                if Hmt > 1:
                    Hmt = Hmt / 100
                    
                # Save parameter
                indParam['Hematocrit'] = Hmt
                paramSource['hematocrit'] = 'User defined'
                
            else:
                paramSource['hematocrit'] = 'Default'
                
            # CREATININE
            if any(persPhysiology.iloc[:,0] == 'creatinine'):
                idxCn = meta.columns == 'creatinine'
                if sum(idxCn) > 1:
                    raise ValueError('More than one column in the metadata contains the parameter creatinine!')
                
                # Extract value and unit
                Cn = meta.loc[sample, 'creatinine']
                currentUnit = persPhysiology.loc[persPhysiology.iloc[:, 0] == 'creatinine', persPhysiology.columns[1]].iloc[0]
        
                # Unit transformation
                if currentUnit == 'mg/dL':
                    indParam['MConUrCreatinineMin'] = Cn * 0.9
                    indParam['MConUrCreatinineMax'] = Cn * 1.1
                    paramSource['mConUrCreatinineMin'] = 'User defined'
                    paramSource['mConUrCreatinineMax'] = 'User defined'
                else:
                    raise ValueError('Creatinine must be provided in mg/dL')
            
            else:
                paramSource['mConUrCreatinineMin'] = 'Default'
                paramSource['mConUrCreatinineMax'] = 'Default'
            
            # CARDIAC OUTPUT
            if any(persPhysiology.iloc[:,0] == 'cardiac output'):
                idxCO = meta.columns == 'cardiac output'
                if sum(idxCO) > 1:
                    raise ValueError('More than one column in the metadata contains the parameter cardiac output!')
                    
                # Extract value and unit
                CO = meta.loc[sample,'cardiac output']
                currentUnit = persPhysiology.loc[persPhysiology.iloc[:, 0] == 'cardiac output', persPhysiology.columns[1]].iloc[0]
                
                # Unit transformation
                if currentUnit == 'L/min':
                    CO = CO * 1000
                elif currentUnit == 'mL/min':
                    CO = CO
                else:
                    raise ValueError('Parameter cardiac output not provided in a valid unit. Available units are [L/min, mL/min]')
                    
                # Save parameter
                indParam['CardiacOutput'] = CO
                paramSource['cardiacOutput'] = 'User defined'
                
            else:
                # Calculate cardiac output
                optionCardiacOutput = 1
                
                if optionCardiacOutput == 1:
                    indParam['CardiacOutput'] = indParam['HeartRate'] * indParam['StrokeVolume']
                    paramSource['cardiacOutput'] = 'Calculated from stroke volume and heart rate'
            
            # GLOMERULAR FILTRATION RATE
            if any(persPhysiology.iloc[:,0] == 'glomerular filtration rate'):
                idxGFR = meta.columns == 'glomerular filtration rate'
                if sum(idxGFR) > 1:
                    raise ValueError('More than one column in the metadata contains the parameter glomerular filtration rate!')
                
                # Extract value and unit
                GFR = meta.loc[sample, 'glomerular filtration rate']
                currentUnit = persPhysiology.loc[persPhysiology.iloc[:, 0] == 'glomerular filtration rate', persPhysiology.columns[1]].iloc[0]
        
                # Unit transformation
                if currentUnit == 'mL/min' or currentUnit == 'mL/min/1.73m^2':
                    indParam['GlomerularFiltrationRate'] = GFR
                    paramSource['glomerularFiltrationRate'] = 'User defined'
                else:
                    raise ValueError('Glomerular Filtration rate must be provided in mL/min/1.73m^2 or mL/min')
                    
            else:
                # Calculate glomerular filtration rate
                RenalFiltrationFactor = 0.2
                
                # Blood flow percentage that Kidney gets
                BK = bloodData.loc[bloodData.loc[:,'Organ'] == 'Kidney', 'Blood Flow Percentage'].iloc[0]
                
                # Save parameter
                RenalFlowRate = BK * indParam['CardiacOutput'] * (1 - indParam['Hematocrit'])
                indParam['GlomerularFiltrationRate'] = RenalFlowRate * RenalFiltrationFactor
                paramSource['glomerularFiltrationRate'] = 'Calculated based on sex, cardiac output and hematocrit'
                
            # Estimate blood volume
            if meta.loc[sample,'sex'] == 'male':
                indParam['bloodVolume'] = (0.3669 * (indParam['Height']/100)**3 + 0.03219 * indParam['bodyWeight'] + 0.6041) * 1000
            elif meta.loc[sample,'sex'] == 'female':
                indParam['bloodVolume'] = (0.3561 * (indParam['Height']/100)**3 + 0.03308 * indParam['bodyWeight'] + 0.1833) * 1000
            
            # Update organ weight and the biomass
            indParam = _calcOrganFract(wbm_model, indParam)
            wbm_model = _adjustWholeBodyRxnCoeff(wbm_model, indParam)
            
            # Update constraints based on new physiological parameters
            pWBM = _physiologicalConstraints(wbm_model, indParam, bloodData, organsListExt)
        
        # Personalising metabolic parameters
        if personalisingMets:
            # Set original parameters to the WBM
            if not personalisingPhys:
                indParam, bloodData, organsListExt = _defaultIndividualParameters(meta.loc[sample, 'sex'])
                indParam = _calcOrganFract(wbm_model, indParam)
                wbm_model = _physiologicalConstraints(wbm_model, indParam, bloodData, organsListExt)
                
                paramSource['bodyWeight'] = 'Default'
                paramSource['age'] = 'Default'
                paramSource['height'] = 'Default'
                paramSource['heartRate'] = 'Default'
                paramSource['strokeVolume'] = 'Default'
                paramSource['hematocrit'] = 'Default'
                paramSource['mConUrCreatinineMin'] = 'Default'
                paramSource['mConUrCreatinineMax'] = 'Default'
                paramSource['cardiacOutput'] = 'Default'
                paramSource['glomerularFiltrationRate'] = 'Default'
            else:
                wbm_model = pWBM
                
            # Personalize metabolite concentrations
            pWBM = _persWBMmetabolomics(sample, wbm_model, meta, persMets, indParam, bloodData, organsListExt)
            
        # Save model and parameters information
        metInfo.to_csv(os.path.join(pers_wbm_path, 'models', sample + '_WBM.metInfo.csv'), index = False)
        rxnInfo.to_csv(os.path.join(pers_wbm_path, 'models', sample + '_WBM.rxnInfo.csv'), index = False)
        with open(os.path.join(pers_wbm_path, 'models', sample + '.pkl'), 'wb') as f:
            pickle.dump(pWBM, f)
        with open(os.path.join(pers_wbm_path, 'params_' + sample + '.pkl'), 'wb') as f:
            pickle.dump(paramSource, f)
            
    else:
        # Load parameters information
        with open(os.path.join(pers_wbm_path, 'params_' + sample + '.pkl'), 'rb') as f:
            paramSource = pickle.load(f)
            
    # Output variable
    return (paramSource)
        
def persWBM(meta: pd.DataFrame, params_file_path: str = '', mets_file_path: str = '', pers_wbm_path: str = os.path.join('.','output','persWBM'),
            solver: str = 'glpk', threads: int = round(os.cpu_count())) -> (pd.DataFrame):
    # Set COBRApy default solver
    try:
        cobra.Configuration().solver = solver
    except:
        print(f'COBRApy could not found installation for {solver}. Simulations will be done with default solver.')
        
    # Generate output folders
    os.makedirs(os.path.join(pers_wbm_path, 'models'), exist_ok = True)

    # Check sex parameter in the meta data
    if any(meta.columns == 'sex'):
        if not meta['sex'].isin(["male", "female"]).all():
            raise ValueError('The possible values for the sex parameter are [male and female]. Check your metadata')
        elif meta['sex'].isin(['male']).all():
            print('All subjects are male')
        elif meta['sex'].isin(['female']).all():
            print('All subjects are female')
        else:
            print('Cohort includes both male and female subjects')
    else:
        raise ValueError('No parameter sex found in metadata')
    
    # Define all physiological parameters
    AllParams = ['age', 'modelID', 'body weight', 'body fat', 'lean body mass', 'height', 'sex', 'heart rate', 'stroke volume',
                 'cardiac output', 'hematocrit', 'creatinine', 'blood flow rate', 'glomerular filtration rate']
    
    # Read physiological parameters data
    if params_file_path != '':
        personalisingPhys = True
        persPhysiology = pd.read_table(params_file_path)
    else:
        personalisingPhys = False
        persPhysiology = []
    
    # Check that all parameters in the file are valid
    if personalisingPhys:
        for param in persPhysiology.iloc[:,0]:
            if param not in AllParams:
                raise ValueError(f'The parameter "{param}" does not match any parameter that can be personalised. Please check spelling and function documentation.')
                
    # Read metabolite concentration data
    if mets_file_path != '':
        personalisingMets = True
        persMets = pd.read_table(mets_file_path)
    else:
        personalisingMets = False
        persMets = []
        
    # If neither physiological or metabolic parameters were specified, extract available parameters in the metadata
    if params_file_path == '' and mets_file_path == '':
        persPhysiology = pd.DataFrame(data = [x for x in meta.columns if x in AllParams], columns = ['Param'])
        if persPhysiology.shape == (0,0):
            raise ValueError(f'No valid parameters found in the meta data. Consider checking the following valid list: {AllParams}')
        else:
            personalisingPhys = True
            
    # Filter physiological parameters with those present in the metadata
    if personalisingPhys:
        idx = persPhysiology.iloc[:,0].isin(meta.columns)
        if all(idx == False):
            personalisingPhys = False
            if personalisingMets:
                print('No physiological parameters found in the metadata. Consider checking spelling. Skipping physiological parameter personalisation.')
            else:
                raise ValueError('No physiological parameters found in the metadata. Consider checking spelling.')
        elif any(idx == False):
            print(f'The following physiological parameters were not found in the metadata and will be discarded: {persPhysiology.loc[~idx,persPhysiology.columns[0]].tolist()}')
            persPhysiology = persPhysiology.loc[idx, :]
    
    # Filter metabolite concentration data with those present in the metadata
    if personalisingMets:
        idx = persMets.iloc[:,0].isin(meta.columns)
        if all(idx == False):
            personalisingMets = False
            if personalisingPhys:
                print('No metabolite IDs found in the metadata. Consider checking spelling. Skipping metabolite concentration personalisation.')
            else:
                raise ValueError('No metabolite IDs found in the metadata. Consider checking spelling.')
        elif any(idx == False):
            print(f'The following metabolite concentrations were not found in the metadata and will be discarded: {persMets.loc[~idx,persMets.columns[0]].tolist()}')
            persMets = persMets.loc[idx, :]
        
    # Create male control model if not already exists and metadata contains male models
    if (not os.path.exists(os.fspath(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvey_1_03d.pkl')))) and any(meta['sex'] == 'Male'):
        print('Generating Harvey_1_03d original model. This could take a while...')
        male = _generateOriginalWBM('Harvey_1_03d', solver, threads)
        with open(str(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvey_1_03d.pkl')), "wb") as f:
            pickle.dump(male, f)
    
    # Create female control model if not already exists and metadata contains female models
    if (not os.path.exists(os.fspath(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvetta_1_03d.pkl')))) and any(meta['sex'] == 'Female'):
        print('Generating Harvetta_1_03d original model. This could take a while...')
        female = _generateOriginalWBM('Harvetta_1_03d', solver, threads)
        with open(str(pkg_resources.files('persephone').joinpath('data/models/WBM/Harvetta_1_03d.pkl')), "wb") as f:
            pickle.dump(female, f)
    
    # Initialize output
    paramInfo = {}
    
    # Run model personalization
    parModelsFunc = partial(_generatePersWBM, meta = meta, personalisingPhys = personalisingPhys, persPhysiology = persPhysiology, personalisingMets = personalisingMets, persMets = persMets, pers_wbm_path = pers_wbm_path)
    with ProcessPoolExecutor(max_workers = threads) as executor:
        # Run function
        futures = {executor.submit(parModelsFunc, sample): sample for sample in meta.index}
        
        # Progress bar
        for future in tqdm(as_completed(futures), total=len(futures), desc="Personalising models"):
            sample = futures[future]
            paramInfo[sample] = future.result()
            
    # Build output dataframe
    paramInfo = pd.DataFrame.from_dict(paramInfo, orient = 'index')
    paramInfo.to_csv(os.path.join(pers_wbm_path, 'persParams.csv'))
    
    # Run model optimization with all available threads
    for sample in tqdm(meta.index, desc = "Optimizing personalised whole body models"):
        # Load model
        with open(os.path.join(pers_wbm_path, 'models', sample + '.pkl'), "rb") as f:
            model = pickle.load(f)
            
        # Setting solver parameters
        #model = _setSolverParameters(model, solver, threads)
        print(model.solver.problem.parameters.threads.get())
        
        # Optimize model
        model.optimize()
        
        # Save model
        with open(os.path.join(pers_wbm_path, 'models', sample + '.pkl'), "wb") as f:
            pickle.dump(model, f)
    
    # Output variable
    return (paramInfo)