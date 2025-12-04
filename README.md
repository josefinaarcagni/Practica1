# persephonepy

Write a brief description of the package here: TO BE DONE

Cite: TO BE DONE

# Installation

You can install this package by cloning this repo and installing manually. To clone the repo:

```
git clone 
```

To install from this repo, check to be into the main directory (running the `ls` command the setup.py file must be displayed), and run:

```
pip install -e .
```

# Configuring solver

Although persephonepy can be used with glpk default solver, it is highly recommended to use more efficient linear optimization solvers such as cplex or gurobi.

## Configuring cplex into conda environment

In order to configure cplex into your conda environment:

```
cd CPLEX_INSTALLATION_PATH/ILOG/CPLEX_StudioXXXX/python
python setup.py install
```

Replace `CPLEX_INSTALLATION_PATH` with the one in your computer and `XXXX` to your cplex version.

## Configuring gurobi into conda environment

In order to configure gurobi into your conda environment:

```
conda install -c gurobi gurobi
export GRB_LICENSE_FILE=GUROBI_LICENSE_FILE_PATH
```

Replace `GUROBI_LICENSE_FILE_PATH` with the path to your gurobi license file.

# Modules 

## Running whole persephone pipeline

In order to run Persephone pipeline, a configuration JSON file is required. An example file can be found at **test** folder. Automatic generation of configuration files can be performed [here](https://vmh2.life/persephone). Then:

```
timeLog = runPersephone(configFile)
```

## Running persephone modules

### MARS

The MARS module of persephonepy translates OTU table and taxa information into suitable format for using different metabolic reconstructions.

**Inputs:**
- counts: OTU table in FeatureTable[Frequency] format.
- taxa: taxonomical assignment file in FeatureData[Taxonomy] format.

**Parameters:**
- database: Metabolic reconstruction. Options are: AGORA2, APOLLO, AGORA2-APOLLO. Default: AGORA2-APOLLO
- output-format: Format of the output files. Default: csv
- taxa-delimiter: Delimiter of the taxonomic assignment file lineage. Default: ;
- flag-lone-species: Boolean to indicate if the name of the genus is included in the species name. Default: False
- cutoff: Value under relative abundances are considered to be 0. Default: 0.000001
- sample-read-counts-cutoff: Value for exclusion of samples when the reads are lower than this value. Default: 1
- remove-clades: Boolean indicating if clade names should be removed from taxa names. Default: True
- output-dir: Path to output dir.

**Outputs:**
- counts-mapped: OTU table in FeatureTable[RelativeFrequency] format mapped to the selected database.
- taxa-mapped: taxonomical assignment file in FeatureData[Taxonomy] format mapped to the selected database.

### MgPipe

After running MARS, the MgPipe module of q2-persephone calculates reaction and subsystem abundances, generates microbial metabolic reconstructions and calculates maximum uptake/secretion rates of some samples.

**Inputs:**
- counts: mapped OTU table in FeatureTable[Frequency] format.
- taxa: mapped taxonomical assignment file in FeatureData[Taxonomy] format.

**Parameters:**
- database: Metabolic reconstruction. Options are: AGORA2, APOLLO, AGORA2-APOLLO. Default: AGORA2-APOLLO
- solver: Linear optimization solver name. Default: glpk
- threads: Number of threads to use for MgPipe. Default: 80% of available threads
- compute-profiles: Boolean to compute uptake/secretion fluxes
- diet-file-name: Name of the diet to contextualize microbial metabolic models during fluxes calculation. Default: EUAverageDiet
- output-dir: Path to output dir.

**Outputs:**
- rxn-ab: Reaction abundances across samples in FeatureTable[Frequency] format.
- rxn-pre: Reaction presence across samples in FeatureTable[PresenceAbsence] format.
- subs-ab: Subsystem abundances across samples in FeatureTable[Frequency] format.

### WBM personalization

This module personalizes whole body models by physiological and metabolomics data from a metadata file. For this function, users could set either physiological or metabolomics files that refer to the information stored in the metadata.

**Inputs:**
- meta: Sample phenotype data. REQUIRED

**Parameters:**
- params-file-path: Path to the file with physiological parameters information.
- mets-file-path: Path to the file with metabolomics information.
- solver: Linear optimization solver name. Default: glpk
- threads: Number of threads to use for MgPipe. Default: 80% of available threads
- pers-wbm-path: Path to the output folder where personalised WBMs will be saved.

**Outputs:**
- par-info: Information about the parameters applied to each sample.

Physiological parameter file should contain the following information in a tab separated format (assuming that parameter values for each sample are stored in the metadata).

| Param      | Unit  |
|------------|--------|
| sex        | nan    |
| age        | years  |
| hematocrit |        |
| height     | ft     |

The file with the metabolic concentrations should contain the following information in a tab separated format (assuming that metabolite concentration values for each sample are stored in the metadata).

| Metabolite | Biofluid | Unit     |
|------------|----------|----------|
| glc_D      | bc       | mg/dL    |
| crtn       | bc       | µmol/L   |
| lac_D      | csf      | mg/dL    |
| cl         | u        | mmol/L   |
| k          | u        | mmol/L   |
| cys_L      | u        | mmol/L   |

### mWBM generation

After generating the microbiome models with MgPipe, this module creates microbiome whole body models for a given list of samples.

**Inputs:**
- meta: Sample phenotype data. REQUIRED

**Parameters:**
- micro-mod-path: Path to the microbiota models generated by MgPipe.
- diet-file-name: Name of the diet to be used. Default: "EUAverageDiet"
- use-pers-wbm: Boolean to indicate if personalized WBMs should be used. Default: False
- pers-wbm-path: Path to the personalized WBMs.
- output-dir: Path to output dir.
- solver: Name of the solver to be used for simulations. Default: glpk
- threads: Number of threads to use for MgPipe. Default: 80% of available threads

**Outputs:**
- stats: Microbiome WBM statistics.

### FBA analysis

Finally, this module computes Flux Balance Analysis to each microbiome whole body model for a given reaction list. The reaction list should be provided in a file with the following format:

| Reaction        |
|------------------|
| DM_glc_D_bc      |
| DM_trp_L_bc      |
| DM_pcresol_bc    |
| DM_lys_L_bc      |
| DM_tmao_bc       |
| DM_tma_bc        |

**Inputs:**
- meta: Sample phenotype data. REQUIRED

**Parameters:**
- rxn-list-path: Path to the file with reaction IDs to be explored.
- mwbm-path: Name of the diet to be used. Default: "EUAverageDiet"
- output-dir: Path to output dir.
- solver: Name of the solver to be used for simulations. Default: glpk
- threads: Number of threads to use for MgPipe. Default: 80% of available threads

**Outputs:**
- fba-res: Results of the FBA.

# Testing

Please refer to **testing.py** script.