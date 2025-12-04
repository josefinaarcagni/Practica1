# -*- coding: utf-8 -*-
"""
Created on Fri Aug 29 20:00:55 2025

@author: tblasco
"""

import pandas as pd
import json
import os
import numpy as np
import re
import logging
import importlib.resources as pkg_resources

def setup_logger(name, log_file, level=logging.INFO):
    formatter = logging.Formatter('%(asctime)s\t%(name)s\t%(levelname)s\t%(filename)s:%(lineno)d\t%(funcName)s()\t%(message)s')

    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(file_handler)
    
    return logger

def read_file_as_dataframe(file_path, header=0, index_col=None):
    """
    Reads a file into a pandas DataFrame based on its extension.

    Args:
        file_path (str): Path to the file.
        header (int): Row number to use as column names (default is 0).
        index_col (int or None): Column to use as the row labels of the DataFrame (default is None).

    Returns:
        pd.DataFrame: A pandas DataFrame containing the data from the file.
    """
    logger = logging.getLogger('main.utils')
    
    # Extract the file extension
    _, extension = os.path.splitext(file_path)
    extension = extension.lower()

    # Map extensions to file types
    try:
        file_extension = file_path.type
    except AttributeError:
        ind = file_path.find('.')
        extension = file_path[ind+1:]
        if extension == 'txt':
            file_extension = 'text/plain'
        elif extension == 'csv':
            file_extension = 'text/csv'
        elif extension == 'xlsx':
            file_extension = 'spreadsheet'
        else:
            file_extension = 'not found'

    # Read the file based on its type
    if file_extension in ['text/plain', 'text/csv', 'text/tab-separated-values']:
        # Assuming text files are tab-separated by default; adjust delimiter as needed
        delimiter = '\t' if file_extension == 'text/plain' else ','
        return pd.read_csv(file_path, sep=delimiter, index_col=index_col, low_memory=False, header=header)
    elif file_extension == 'spreadsheet'or file_extension == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet':
        return pd.read_excel(file_path, engine='openpyxl', index_col=index_col, header=header)
    elif file_extension == '.parquet':
        return pd.read_parquet(file_path)
    else:
        raise ValueError(f"Unsupported file type: {file_extension}")


def merge_files(input1, input2, input_is_df=False):
    """
    Merge two files using pandas and return the merged DataFrame.

    Args:
        input1 (str or pd.DataFrame): The first input, either a file path or a DataFrame.
        input2 (str or pd.DataFrame): The second input, either a file path or a DataFrame, or None.
        input_is_df (bool): Whether the inputs are DataFrames instead of file paths.

    Returns:
        pd.DataFrame: The merged DataFrame.
    """
    logger = logging.getLogger('main.utils')
    
    # Read input files into pandas DataFrames  
    if input2 is None:
        if not input_is_df:
            df1 = read_file_as_dataframe(input1, header=0, index_col=0)
        else:
            df1 = input1

        merged_df = df1

        if merged_df.index.name != 'Taxon':
            if 'Taxon' in merged_df.columns:
                merged_df = merged_df.set_index('Taxon')
    else:
        if not input_is_df:
            df1 = read_file_as_dataframe(input1, header=0, index_col=0)
            df2 = read_file_as_dataframe(input2, header=1, index_col=0)
        else:
            df1 = input1
            df2 = input2

        # Merge DataFrames using their index values
        merged_df = pd.merge(df1, df2, left_index=True, right_index=True, how='inner')

        try:
            merged_df = merged_df.drop(columns=['Confidence'])
        except KeyError:
            pass

        merged_df = merged_df.reset_index(drop=True)

        if 'Taxon' in merged_df.columns:
            merged_df = merged_df.set_index('Taxon')

    return merged_df


def normalize_dataframe(dataframe, dfvalues_are_rel_abundances=False, cutoff=0.000001, dataframe_to_normalize_to=None):
    """
    Normalize a dataframe by first grouping and summing rows with the same name,
    and then calculating the relative abundances per taxa so that the sum of each sample (each column) is 1.
    Optionally, a cutoff can be provided to filter out low abundant taxa after normalization, as well as 
    an additional dataframe containing total read counts for all samples of the input dataframe, to which, in case
    being provided, will be normalized to instead.

    Args:
        dataframe (pd.Dataframe):                           The input DataFrame with taxonomic groups in the index and read counts per sample.
        dfvalues_are_rel_abundances (boolean, optional):    Flag indicating whether read counts are already normalized and therefore
                                                            only the cutoff will be applied.
                                                            Defaults to False.
        cutoff (float, optional):                           A cut-off value for filtering out low abundance taxa. Defaults to None.
        dataframe_to_normalize_to (pd.Dataframe, optional): Containing per-sample total read counts pre-mapping,
                                                            allowing for taxa abundance normalization against pre-mapped total read counts
                                                            instead of normalization against only present or absent total read counts.
                                                            Defaults to None.

    Returns:
        grouped_df_afterCutoff (pd.Dataframe): The input DataFrame with taxonomic groups in the index and total read counts per sample, after cutoff was applied.
        rel_abundances_df_afterCutoff (pd.Dataframe): The input DataFrame with taxonomic groups in the index and normalized relative abundances per sample, after cutoff was applied.
    """
    logger = logging.getLogger('main.utils')
    logger_taxa_below_cutoff = logging.getLogger('taxa_below_cutoff.utils')
    
    # Group by index and sum the rows with the same name
    grouped_df = dataframe.groupby(dataframe.index.name).sum()
    
    # Normalize each column so that the sum of each column is 1 (either
    # to pre-mapped total read counts, or to the subset read counts for 
    # the dataset with taxa present in model database - which needs to be done for modelling to work)
    if dataframe_to_normalize_to is not None:
        read_counts = dataframe_to_normalize_to.sum()
    else:
        read_counts = grouped_df.sum()
    
    # Normalize read counts to get relative abundances of taxa
    rel_abundances_df = grouped_df.div(read_counts)

    # Apply cut-off for low abundance taxa
    if cutoff != 0:
        # Identify which taxa are below relative abundance cutoff threshold & set their read counts to 0
        grouped_df[rel_abundances_df <= cutoff] = 0

        # Remove taxa which are non-abundant in any sample after cutoff has been applied
        grouped_df_afterCutoff = grouped_df[(grouped_df != 0).any(axis=1)]

        # Renormalize the dataframe after cutoff was applied so that relative abundances per sample sum to 1
        read_counts_afterCutoff = grouped_df_afterCutoff.sum()
        rel_abundances_df_afterCutoff = grouped_df_afterCutoff.div(read_counts_afterCutoff)

        # Log which taxa are below the cutoff in which samples and were set to 0
        entries_below_cutoff = rel_abundances_df[rel_abundances_df <= cutoff].stack().index.tolist()
        if entries_below_cutoff:
            logger.info(f"Taxa were below the cutoff & are listed in seperate log-file.")
            logger_taxa_below_cutoff.info(f"Taxa whose rel.abundance was below the cutoff & therefore set to 0: {entries_below_cutoff}")
        else:
            logger.info(f"No taxa were below the cutoff.")
    else:
        grouped_df_afterCutoff = grouped_df
        rel_abundances_df_afterCutoff = rel_abundances_df
        
    return grouped_df_afterCutoff, rel_abundances_df_afterCutoff


def combine_metrics(metrics1, metrics2, df_type="metrics", dfvalues_are_rel_abundances=False):
    """
    Combine the metrics from two different sets of taxonomic DataFrames into a single DataFrame for each level.

    Args:
        metrics1 (dict): A dictionary with keys as taxonomic levels and values as the calculated metrics for the first group.
        metrics2 (dict): A dictionary with keys as taxonomic levels and values as the calculated metrics for the second group.

    Returns:
        dict: A dictionary with keys as taxonomic levels and values as the combined DataFrames.
    """
    logger = logging.getLogger('main.utils')
    
    combined_metrics = {}

    for level in metrics1.keys():
        level_metrics_pre_mapping = metrics1[level]
        level_metrics_post_mapping = metrics2[level]

        if df_type == "metrics":
            combined_level_metrics = {}

            for metric_name in level_metrics_pre_mapping.keys():
                
                combined_metric = pd.DataFrame([level_metrics_pre_mapping[metric_name], level_metrics_post_mapping[metric_name]])
                combined_metric.index = ['Pre mapping', 'Post mapping']
                combined_level_metrics[metric_name] = combined_metric

            combined_metrics[level] = combined_level_metrics

        elif df_type == "summ_stats":
            combined_metric = pd.DataFrame({'Pre mapping': level_metrics_pre_mapping, \
                                            'Post mapping': level_metrics_post_mapping})
            
            combined_metric['Mapping coverage'] = np.nan
            combined_metric.loc[0, 'Mapping coverage'] = combined_metric.loc[0, 'Post mapping'] / combined_metric.loc[0, 'Pre mapping']
            combined_metric.loc[2, 'Mapping coverage'] = combined_metric.loc[2, 'Post mapping'] / combined_metric.loc[2, 'Pre mapping']
            combined_metric.loc[4, 'Mapping coverage'] = combined_metric.loc[4, 'Post mapping'] - combined_metric.loc[4, 'Pre mapping']
            combined_metric.loc[6, 'Mapping coverage'] = combined_metric.loc[6, 'Post mapping'] / combined_metric.loc[6, 'Pre mapping']
            
            if dfvalues_are_rel_abundances == False:
                combined_metric['Description'] = ['The number of taxa across all samples. MappingCoverage = post mapping/pre mapping.', \
                                                  'The estimated number of taxa following standard nomenclature. Estimated by excluding all taxa whose names contain "-" &/or multiple uppercase letters in a row.', \
                                                  'Mean number of taxa across samples (equals to species richness on species taxonomic level). MappingCoverage = post mapping/pre mapping.', \
                                                  'Standard deviation of taxa richness.', 'Mean alpha-diversity (calc. by pielous evenness) in samples. Towards 0: low diversity, towards 1: high diversity, with 1: complete evenness. MappingCoverage = post mapping - pre mapping.', \
                                                 'Standard deviation of pielous evenness.', 'Mean number of reads across samples. MappingCoverage = post mapping/pre mapping.', \
                                                 'Standard deviation of read counts.']
    
                combined_metric.index = ['Total number of taxa', 'Estimated total number of named taxa', 'Mean taxa richness', \
                                         'Std taxa richness', 'Mean pielous evenness', 'Std pielous evenness', \
                                         'Mean read counts', 'Std read counts']
            else:
                combined_metric['Description'] = ['The number of taxa across all samples. MappingCoverage = post mapping/pre mapping.', \
                                                  'The estimated number of taxa following standard nomenclature. Estimated by excluding all taxa whose names contain "-" &/or multiple uppercase letters in a row.', \
                                                  'Mean number of taxa across samples (equals to species richness on species taxonomic level). MappingCoverage = post mapping/pre mapping.', \
                                                  'Standard deviation of taxa richness.', 'Mean alpha-diversity (calc. by pielous evenness) in samples. Towards 0: low diversity, towards 1: high diversity, with 1: complete evenness. MappingCoverage = post mapping - pre mapping.', \
                                                 'Standard deviation of pielous evenness.', 'Mean relative abundance across samples. MappingCoverage = post mapping/pre mapping.', \
                                                 'Standard deviation of relative abundance.']

                combined_metric.index = ['Total number of taxa', 'Estimated total number of named taxa', 'Mean taxa richness', \
                                         'Std taxa richness', 'Mean shannon index', 'Std shannon index', \
                                         'Mean relative abundance', 'Std relative abundance']
            
            combined_metrics[level] = combined_metric

    return combined_metrics


def save_dataframes(dataframe_groups, output_path, output_format):
    """
    Save DataFrames to the specified output path in the given format.

    Args:
        dataframe_groups (dict): A dictionary of dictionaries containing DataFrames.
                                  Outer dictionary keys represent group names, and inner dictionaries contain
                                  DataFrames for each taxonomic level.
        output_path (str): The directory path where the output files will be saved.
        output_format (str): The format to save the DataFrames ('csv', 'excel', 'parquet', or 'json').
    """
    logger = logging.getLogger('main.utils')
    
    os.makedirs(output_path, exist_ok=True)

    for group_name, values in dataframe_groups.items():
        group_output_path = os.path.join(output_path, group_name)
        os.makedirs(group_output_path, exist_ok=True)
        if group_name == "metrics":
            for level, metrics_dataframes in values[0].items():
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                for metric_name, df in metrics_dataframes.items():
                    file_name = f"{metric_name}.{output_format}"
                    file_path = os.path.join(level_output_path, file_name)

                    if output_format == "csv":
                        df.to_csv(file_path)
                    elif output_format == "txt" or output_format == "tsv":
                        df.to_csv(file_path, sep='\t')
                    elif output_format == "excel":
                        df.to_excel(file_path)
                    elif output_format == "parquet":
                        df.to_parquet(file_path)
                    elif output_format == "json":
                        df.to_json(file_path)
                    else:
                        raise ValueError(f"Unsupported output format: {output_format}")

            for level, summ_stats_df in values[1].items():
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"summary_statistics_{level}.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    summ_stats_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    summ_stats_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    summ_stats_df.to_excel(file_path)
                elif output_format == "parquet":
                    summ_stats_df.to_parquet(file_path)
                elif output_format == "json":
                    summ_stats_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")

            for level, abundance_metrics_df in values[2].items():
                # Replace ' ' by '_' in taxa names
                abundance_metrics_df.index = abundance_metrics_df.index.str.replace(' ', '_')

                # Save dataframes
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"preMapping_abundanceMetrics_{level}.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    abundance_metrics_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    abundance_metrics_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    abundance_metrics_df.to_excel(file_path)
                elif output_format == "parquet":
                    abundance_metrics_df.to_parquet(file_path)
                elif output_format == "json":
                    abundance_metrics_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")

            for level, abundance_metrics_df in values[3].items():
                # Replace ' ' by '_' in taxa names
                abundance_metrics_df.index = abundance_metrics_df.index.str.replace(' ', '_')

                # Save dataframes
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"mapped_abundanceMetrics_{level}.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    abundance_metrics_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    abundance_metrics_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    abundance_metrics_df.to_excel(file_path)
                elif output_format == "parquet":
                    abundance_metrics_df.to_parquet(file_path)
                elif output_format == "json":
                    abundance_metrics_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")
    
            for level, abundance_metrics_df in values[4].items():
                # Replace ' ' by '_' in taxa names
                abundance_metrics_df.index = abundance_metrics_df.index.str.replace(' ', '_')

                # Save dataframes
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"unmapped_abundanceMetrics_{level}.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    abundance_metrics_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    abundance_metrics_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    abundance_metrics_df.to_excel(file_path)
                elif output_format == "parquet":
                    abundance_metrics_df.to_parquet(file_path)
                elif output_format == "json":
                    abundance_metrics_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")

            for level, beta_diversity_df in values[5].items():
                # Save dataframes
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"preMapping_brayCurtisDissimilarity_{level}.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    beta_diversity_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    beta_diversity_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    beta_diversity_df.to_excel(file_path)
                elif output_format == "parquet":
                    beta_diversity_df.to_parquet(file_path)
                elif output_format == "json":
                    beta_diversity_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")
        
            for level, beta_diversity_df in values[6].items():
                # Save dataframes
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"mapped_brayCurtisDissimilarity_{level}.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    beta_diversity_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    beta_diversity_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    beta_diversity_df.to_excel(file_path)
                elif output_format == "parquet":
                    beta_diversity_df.to_parquet(file_path)
                elif output_format == "json":
                    beta_diversity_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")
        
            for level, beta_diversity_df in values[7].items():
                # Save dataframes
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"unmapped_brayCurtisDissimilarity_{level}.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    beta_diversity_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    beta_diversity_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    beta_diversity_df.to_excel(file_path)
                elif output_format == "parquet":
                    beta_diversity_df.to_parquet(file_path)
                elif output_format == "json":
                    beta_diversity_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")

        elif "_stratified_metrics" in group_name:
            for level, metrics_dataframes in values[0].items():
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                for metric_name, df in metrics_dataframes.items():
                    file_name = f"{metric_name}_stratified.{output_format}"
                    file_path = os.path.join(level_output_path, file_name)

                    if output_format == "csv":
                        df.to_csv(file_path)
                    elif output_format == "txt" or output_format == "tsv":
                        df.to_csv(file_path, sep='\t')
                    elif output_format == "excel":
                        df.to_excel(file_path)
                    elif output_format == "parquet":
                        df.to_parquet(file_path)
                    elif output_format == "json":
                        df.to_json(file_path)
                    else:
                        raise ValueError(f"Unsupported output format: {output_format}")

            for level, summ_stats_df in values[1].items():
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"summary_statistics_{level}_stratified.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    summ_stats_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    summ_stats_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    summ_stats_df.to_excel(file_path)
                elif output_format == "parquet":
                    summ_stats_df.to_parquet(file_path)
                elif output_format == "json":
                    summ_stats_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")

            for level, abundance_metrics_df in values[2].items():
                # Replace ' ' by '_' in taxa names
                abundance_metrics_df.index = abundance_metrics_df.index.str.replace(' ', '_')

                # Save dataframes
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"preMapping_abundanceMetrics_{level}_stratified.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    abundance_metrics_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    abundance_metrics_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    abundance_metrics_df.to_excel(file_path)
                elif output_format == "parquet":
                    abundance_metrics_df.to_parquet(file_path)
                elif output_format == "json":
                    abundance_metrics_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")

            for level, abundance_metrics_df in values[3].items():
                # Replace ' ' by '_' in taxa names
                abundance_metrics_df.index = abundance_metrics_df.index.str.replace(' ', '_')

                # Save dataframes
                level_output_path = os.path.join(group_output_path, level)
                os.makedirs(level_output_path, exist_ok=True)

                file_name = f"mapped_abundanceMetrics_{level}_stratified.{output_format}"
                file_path = os.path.join(level_output_path, file_name)

                if output_format == "csv":
                    abundance_metrics_df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    abundance_metrics_df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    abundance_metrics_df.to_excel(file_path)
                elif output_format == "parquet":
                    abundance_metrics_df.to_parquet(file_path)
                elif output_format == "json":
                    abundance_metrics_df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")
    
        else:
            for level, df in values.items():
                # Replace ' ' by '_' in taxa names
                df.index = df.index.str.replace(' ', '_')

                # Save dataframes
                file_name = f"{group_name}_{level.lower()}.{output_format}"
                file_path = os.path.join(group_output_path, file_name)

                if output_format == "csv":
                    df.to_csv(file_path)
                elif output_format == "txt" or output_format == "tsv":
                    df.to_csv(file_path, sep='\t')
                elif output_format == "excel":
                    df.to_excel(file_path)
                elif output_format == "parquet":
                    df.to_parquet(file_path)
                elif output_format == "json":
                    df.to_json(file_path)
                else:
                    raise ValueError(f"Unsupported output format: {output_format}")
            


def load_input_and_preprocess(input_file1, input_file2=None, taxaSplit=';', input_is_df=False):
    """
    Loads in input data, and preprocesses it (merging dataframes, handling NaN,
    targeting duplicate entries, subsetting for taxa with species information &
    checking whether input is in absolute read counts or already normalized to relative abundances).

    Args:
        input_file1 (chars/string):     Path to input data (can be tab seperated, comma seperated, etc.),
                                        containing either only read counts or additionally already a index column
                                        with taxa names.
        input_file1 (chars/string):     Path to additional input data (can be tab seperated, comma seperated, etc.),
                                        containing taxa names associated with input_file1.
                                        Defaults to None, only required if input_file1 does not contain taxa names.
        taxaSplit (string):             Seperator by which input taxonomic levels are seperated.
        input_is_df (boolean):          Boolean whether input file is a path to a file or (in case o the streamlit
                                        app.py) a pandas dataframe.

    Returns:
        uniqueSpecies_dataframe (pd.dataframe):     A dataframe containing the preprocessed data from input_file1
                                                    (and potentially input_file2).
        dfvalues_are_rel_abundances (boolean):      A boolean indicating whether the input data is already normalized
                                                    to relative abundances or not. Relevant for MARS subsequent steps.
    """
    logger = logging.getLogger('main.operations')
    logger.info("Loading input data & preprocessing it.")

    # Read in input file(s) & merge dataframes in case abundances plus taxa names are stored seperately
    merged_dataframe = merge_files(input_file1, input_file2, input_is_df=input_is_df)

    # Replace NaN by 0
    merged_dataframe_woNaN = merged_dataframe.fillna(0)

    # Sum read counts (or relative abundances) for same taxa in case there are duplicate entries
    uniqueTaxa_dataframe = merged_dataframe_woNaN.groupby(merged_dataframe_woNaN.index.name).sum()
    print(uniqueTaxa_dataframe)

    # Add ["k__", "p__", "c__", "o__", "f__", "g__", "s__"] to taxa names to indicate taxonomic levels, if not present already
    # & remove any leading whitespaces in front of taxa names per taxonomic level
    uniqueTaxa_dataframe.index = uniqueTaxa_dataframe.index.map(lambda x: standardize_prefixes(x, taxaSplit=taxaSplit))
    
    # Remove rows/taxa from dataframe which do not contain species level information (& turned None in previous step)
    uniqueSpecies_dataframe = uniqueTaxa_dataframe[uniqueTaxa_dataframe.index.notnull()]

    if uniqueSpecies_dataframe.empty:
        raise ValueError("Either no species present in input data or taxonomic level indication in naming does not follow MARS convention. If so, please adapt the taxa naming according to the template or guidance from the tutorials.")
    
    # Replace potential "_" between Genus & epithet in species name by " "
    uniqueSpecies_dataframe.index = uniqueSpecies_dataframe.index.str.replace(r'(?<!_)_(?!_)', ' ', regex=True)

    # Check whether the dataframe contains absolute read counts or relative abundances and 
    # set boolean accordingly for subsequent steps
    dfvalues_are_rel_abundances = check_df_absolute_or_relative_counts(uniqueSpecies_dataframe)

    return uniqueSpecies_dataframe, dfvalues_are_rel_abundances


def standardize_prefixes(df_index, taxaSplit=';'):
    """
    Add ["k__", "p__", "c__", "o__", "f__", "g__", "s__"] to taxa names to 
    indicate taxonomic levels, if "s__" prefix is not present already. If "s__"
    is present, reorder the taxonomic levels which are present, into standard order
    from highest taxonomic level (kingdom) to lowest (species) to have a standard format
    for subsequent MARS steps & replaces "d__" by "k__".
    Any leading whitespaces in front of taxa names are also removed.
    """
    logger = logging.getLogger('main.operations')
    
    prefixes = ["k__", "p__", "c__", "o__", "f__", "g__", "s__"]
    parts = df_index.split(taxaSplit)

    # Replace "d__" with "k__" in the input
    df_index = df_index.replace("d__", "k__")
    parts = [part.replace("d__", "k__") for part in parts]

    if "s__" not in df_index and len(parts) == 7:
        # Remove leading whitespace from each part
        parts = [part.lstrip() for part in parts]
        
        # Combine prefixes with cleaned parts
        prefixed_parts = []
        for prefix, part in zip(prefixes, parts):
            if part:  # Only add prefix if part is not empty
                prefixed_parts.append(f"{prefix}{part}")
            else:
                prefixed_parts.append("")  # Keep empty parts as empty strings
        modified_taxaNames = taxaSplit.join(prefixed_parts)

        return modified_taxaNames
    elif "s__" in df_index:
        # If 's__' is already present, remove leading whitespace from each part
        cleaned_parts = [part.lstrip() for part in parts]

        # Create a dictionary to store parts by their prefix
        sorted_parts = {prefix: "" for prefix in prefixes}
        
        # Assign parts to their corresponding prefixes
        for part in cleaned_parts:
            for prefix in sorted_parts.keys():
                if part.startswith(prefix):
                    sorted_parts[prefix] = part
                    break
        
        modified_taxaNames = taxaSplit.join(sorted_parts[prefix] for prefix in prefixes if sorted_parts[prefix])
        
        return modified_taxaNames
    else:
        return None


def check_df_absolute_or_relative_counts(df):
    """
    Checks whether dataframe read counts are absolute read counts or were 
    already normalized to relative abundances. Sets boolean, which is used
    in subsequent steps to avoid re-normalization & unwanted metrics calculation.

    Arg:
        df (pd.DataFrame): The input DataFrame with taxonomic groups in the index.
    
    Returns:
        dfvalues_are_rel_abundances (boolean): Boolean stating if dataframe values are normalized relative abundances or not.
    """
    logger = logging.getLogger('main.operations')
    
    read_counts = df.sum()
    dfvalues_are_rel_abundances = True
    if (read_counts > 1).any():
        dfvalues_are_rel_abundances = False
    
    logger.info(f"Evaluating whether input data was already normalized to relative abundances. Already normalized: {dfvalues_are_rel_abundances}.")

    return dfvalues_are_rel_abundances


def remove_clades_from_taxaNames(preprocessed_dataframe, taxaSplit=';'):
    """ 
    Removes clade extensions from taxonomic names at any taxonomic level (if present)
    and sums counts of all clades of each taxa together, because most taxa with clade extensions
    will not be present in AGORA2/APOLLO and would therefore be set to absent, but their
    according taxa without the clade extension could be present.

    Args:
        preprocessed_dataframe (pd.DataFrame):  The input DataFrame with taxonomic groups in the index.
        taxaSplit (string):                     Seperator by which input taxonomic levels are seperated.
    
    Returns:
        grouped_df (pd.DataFrame): The input DataFrame with taxonomic groups in the index, without clade seperation anymore.
    """
    logger = logging.getLogger('main.operations')
    logger.info("Removing clade extensions from taxa names.")
    
    preprocessed_dataframe = preprocessed_dataframe.reset_index()
    taxa = preprocessed_dataframe['Taxon']
    
    # Abbreviations for all taxonomic levels
    taxUnitAbbreviations = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    # Specifies the column for the current taxonomic level from the split name (increases by 1 with each for-loop iteration)
    columnTracker = 0

    for taxLevel in taxUnitAbbreviations:
        # Filter taxa that contain taxLevel & split the strings by the specified taxaSplit
        taxa_wTaxLevel = taxa[taxa.str.contains(taxLevel)]
        if len(taxa_wTaxLevel) > 0:
            taxaSep = taxa_wTaxLevel.str.split(taxaSplit, expand=True)
        
            # Extract taxa names from the taxLevel column
            taxName = taxaSep[columnTracker].astype(str)
            
            # Remove clades from taxa names using regex
            original_taxName = taxName.copy()
            taxName = taxName.str.replace(r'(?i)\s+(clade\s+)?[A-Z]$', '', regex=True)
            taxName = taxName.str.replace(r'(?i)\s+(clade\s+)?[A-Z]\s', '', regex=True)
        
            # Find the taxa names that were changed
            changed_taxa = original_taxName[original_taxName != taxName]
            
            # Log all taxa for which clades were found & removed, if any were found
            if not changed_taxa.empty:
                logger_output = ", ".join(changed_taxa.values)
                logger.info(f"List of {taxLevel} taxa for which clade extension was found & removed: {logger_output}.")
            else:
                logger.info(f"For {taxLevel} taxa no clade extension was found & removed.")
                        
            # Update full taxa name (inlcuding all taxLevels) with cleaned taxa
            if taxLevel == 's__':
                taxaUpdate = pd.concat([taxaSep.iloc[:, :columnTracker], taxName], axis=1)
            else:
                taxaUpdate = pd.concat([taxaSep.iloc[:, :columnTracker], taxName, taxaSep.iloc[:, columnTracker+1:]], axis=1)
            taxaUpdate = taxaUpdate.apply(lambda x: taxaSplit.join(x.astype(str)), axis=1)
        
            # Replace original taxa with updated taxa names & update the merged_df
            taxa[taxa.str.contains(taxLevel)] = taxaUpdate.values

            # Increase column to match corresponding taxonomic level for next for-loop iteration
            columnTracker += 1
    
    # Replace the clade-containing taxa names by the updated taxa names without clades
    preprocessed_dataframe['Taxon'] = taxa

    # Group by 'Taxon' and sum, remove 'GroupCount' column if it exists &
    # rename columns to remove 'sum_' prefix
    grouped_df = preprocessed_dataframe.groupby('Taxon').sum().reset_index()
    if 'GroupCount' in grouped_df.columns:
        grouped_df.drop(columns='GroupCount', inplace=True)
    grouped_df.columns = [col.replace('sum_', '') for col in grouped_df.columns]

    grouped_df = grouped_df.set_index('Taxon')

    return grouped_df


def concatenate_genus_and_species_names(preprocessed_dataframe, taxaSplit=';'):
    """ 
    Concatenates genus name with species epithet to form full species names
    and replaces the species epithet in the dataframe index with the full name
    for subsequent steps.

    Args:
        preprocessed_dataframe (pd.DataFrame):  The input DataFrame with taxonomic groups in the index.
        taxaSplit (string):                     Seperator by which input taxonomic levels are seperated.
    
    Returns:
        grouped_df (pd.DataFrame): The input DataFrame with taxonomic groups in the index, with complete species names.
    """
    logger = logging.getLogger('main.operations')
    logger.info("Concatenating genus name with species epithet.")
    
    preprocessed_dataframe = preprocessed_dataframe.reset_index()
    taxa = preprocessed_dataframe['Taxon']

    # Filter taxa that contain species names & split the strings by the specified taxaSplit
    taxa_wSpecies = taxa[taxa.str.contains('s__')]
    if len(taxa_wSpecies) > 0:
        taxaSep = taxa_wSpecies.str.split(taxaSplit, expand=True)
        
        # Get number of columns in taxaSep to identify number of taxonomic level present 
        # (used to find position of genus & species in taxaSep)
        numTaxLevels = taxaSep.shape[1]
    
        # Extract genus name & species epithet
        genusName = taxaSep[numTaxLevels-2].astype(str)
        speciesName = taxaSep[numTaxLevels-1].astype(str)

        genusName = genusName.str.replace('g__', '')
        speciesName = speciesName.str.replace('s__', '')

        # Merge genus name & species epithet to form updated species name
        new_speciesName = 's__' + genusName + ' ' + speciesName

        # Update full taxa name (inlcuding all taxLevels) with cleaned taxa
        taxaUpdate = pd.concat([taxaSep.iloc[:, :numTaxLevels-1], new_speciesName], axis=1)
        taxaUpdate = taxaUpdate.apply(lambda x: taxaSplit.join(x.astype(str)), axis=1)
    
        # Replace original taxa with updated taxa names & update the merged_df
        taxa[taxa.str.contains('s__')] = taxaUpdate.values
    
    # Replace the clade-containing taxa names by the updated taxa names without clades
    preprocessed_dataframe['Taxon'] = taxa

    # Group by 'Taxon' and sum, remove 'GroupCount' column if it exists &
    # rename columns to remove 'sum_' prefix
    grouped_df = preprocessed_dataframe.groupby('Taxon').sum().reset_index()
    if 'GroupCount' in grouped_df.columns:
        grouped_df.drop(columns='GroupCount', inplace=True)
    grouped_df.columns = [col.replace('sum_', '') for col in grouped_df.columns]

    grouped_df = grouped_df.set_index('Taxon')

    return grouped_df


def rename_taxa(preprocessed_dataframe):
    """
    Rename taxa in the preprocessed dataframe by applying alterations and homosynonyms.

    Args:
        preprocessed_dataframe (pandas dataframe):The preprocessed input DataFrame with taxonomic groups in the index.

    Returns:
        renamed_dataframe: The input DataFrame with taxonomic groups in the index, which in case they had initially
                            a different naming convention than AGORA2/APOLLO, are renamed.
    """
    logger = logging.getLogger('main.operations')
    logger.info("Renaming taxa for database compatibility.")

    renaming_json_path = pkg_resources.files('persephone').joinpath('data/sp/renaming.json')

    # Read the dictionaries from the JSON file
    with open(renaming_json_path, 'r') as f:
        loaded_dicts = json.load(f)

    # Access the dictionaries
    alterations, specific_alterations, homosynonyms = loaded_dicts

    renamed_df = preprocessed_dataframe.copy()

    # Apply alterations
    for pattern in alterations:
        renamed_df.index = renamed_df.index.str.replace(pattern, '', regex=True)

    # Apply specific alterations
    for pattern, replacement in specific_alterations.items():
        renamed_df.index = renamed_df.index.str.replace(pattern, replacement, regex=True)

    # Apply homosynonyms
    for pattern, replacement in homosynonyms.items():
        renamed_df.index = renamed_df.index.str.replace(pattern, replacement, regex=True)
    
    # Identify & log replaced entries
    replaced_entries = preprocessed_dataframe.index[preprocessed_dataframe.index != renamed_df.index]
    replacement_pairs = []
    for original_taxa_name, replacement in zip(replaced_entries, renamed_df.index[preprocessed_dataframe.index != renamed_df.index]):
        replacement_pairs.append(f"Original taxa name:{original_taxa_name} - Replacement: {replacement}")
    
    if replacement_pairs:
        logger_output = ', '.join(replacement_pairs)
        logger.info(f"Original taxa name(s) with their replacement(s): {logger_output}.")
    else:
        logger.info(f"No taxa namings were replaced.")
    
    # Group by index and sum the rows with the same name (as some taxa might be non-unique after renaming)
    renamed_df = renamed_df.groupby(renamed_df.index.name).sum()

    return renamed_df


def filter_samples_low_read_counts(dataframe, sample_read_counts_cutoff=1):
    """
    Filters the renamed dataframe (in case it contains absolute read counts) 
    for samples which contain less total read counts than a specified cutoff, with min.=1
    & exclude them from the output dataframes. This ensures there are no samples with read counts = 0, which
    could lead to downstream errors in the pipeline, as well as that the sequencing depth
    is high enough in each sample that saturation for new species detection is reached.

    Args:
        dataframe (pandas dataframe):       The input DataFrame with taxonomic groups in the index.
        sample_read_counts_cutoff (int):    A cutoff for minimal read counts in a sample to be included in downstream analysis.
                                            Defaults to 1, and is min = 1.

    Returns:
        filtered_dataframe (pandas dataframe): The input DataFrame with taxonomic groups in the index,
                                            without samples whose total read count is below the threshold.
    """
    logger = logging.getLogger('main.operations')
    logger.info('Filtering the merged dataframe for samples with read counts below the cutoff.')
    
    # If sample_read_counts_cutoff is lower than 1, set it to 1 to ensure that there are no samples
    # with 0 read counts, which would cause MgPipe to crash
    if sample_read_counts_cutoff < 1:
        logger.warning('The sample_read_counts_cutoff was tried to set below 1, and therefore replaced by 1, as this is the required minimum of read counts in a sample for Persephone to run properly.')
        sample_read_counts_cutoff = 1
    
    read_counts = dataframe.sum()
    
    # Subset the dataframe only including those samples with read count higher than cutoff
    samples_equal_or_higher_than_cutoff = read_counts[read_counts >= sample_read_counts_cutoff].index
    if len(samples_equal_or_higher_than_cutoff) > 0:
        filtered_dataframe = dataframe[samples_equal_or_higher_than_cutoff]

        # Log which samples have a total read count below the cutoff and are therefore removed
        samples_lower_than_cutoff = read_counts[read_counts < sample_read_counts_cutoff].index
        if len(samples_lower_than_cutoff) > 0:
            logger_output = ', '.join(samples_lower_than_cutoff)
            logger.info(f"Following samples had a total read count below the cutoff & were removed: {logger_output}.")
        else:
            logger.info("No samples were below the read counts cutoff & removed.")
    else:
        raise too_low_read_counts_error()
    
    return filtered_dataframe


class too_low_read_counts_error(Exception):
    """
    Exception raised in function 'filter_samples_low_read_counts' when all 
    samples have read counts below the threshold and subsequent analysis can not be run.
    """
    def __init__(self, message="All samples have read counts below the sample_read_counts_cutoff (which defaults to 1, with min=1). Therefore, no subsequent analysis can be run."):
        self.message = message
        super().__init__(self.message)


def check_presence_in_modelDatabase(renamed_dataframe, whichModelDatabase="full_db", userDatabase_path="", taxaSplit=';'):
    """
    Check if entries from the input DataFrame are in the model-database DataFrame under the same taxonomic level column.
    Split the input DataFrame into two DataFrames: present and absent. 

    Args:
        renamed_dataframe (pandas dataframe):   The input DataFrame to be 
                                                checked against the model-database (which is
                                                AGORA2 or APOLLO, combination of both or a user-defined one).
        whichModelDatabase (string):            A string defining if AGORA2, APOLLO, a 
                                                combination of both or a user-defined database should be used as model
                                                database to check presence in. 
                                                Allowed Input (case-insensitive): "AGORA2", "APOLLO", "full_db", "user_db".
                                                Default: "full_db".
        userDatabase_path (string):             A string containing the full path to the user-defined database,
                                                which should be in .csv, .txt, .parquet or .xlsx format and
                                                have column names = taxonomic levels.

    Returns:
        pandas dataframes: Two dataframes containing present and absent taxa.
    """
    logger = logging.getLogger('main.operations')
    logger.info('Checking presence of taxa in model database.')

    # Check, if user wants to use the AGORA2 or APOLLO database, the combination of both or
    # a user-defined database & load the according one
    dataPath = pkg_resources.files('persephone').joinpath('data/sp/')
    modelDatabase_df = pd.read_parquet(dataPath.joinpath(whichModelDatabase + '.parquet'))
    
    if whichModelDatabase.lower() == "agora2":
        # Read in model-database as dataframe
        updatedModelDatabase_df = modelDatabase_df[modelDatabase_df['Resource'] == 'AGORA2'].drop('Resource', axis=1)
    elif whichModelDatabase.lower() == "apollo":
        # Read in model-database as dataframe
        updatedModelDatabase_df = modelDatabase_df[modelDatabase_df['Resource'] == 'APOLLO'].drop('Resource', axis=1)
    elif whichModelDatabase.lower() == "user_db":
        # Read in model-database as dataframe
        updatedModelDatabase_df = read_file_as_dataframe(userDatabase_path, header=0)               
    else:
        # Read in model-database as dataframe
        updatedModelDatabase_df = modelDatabase_df.drop('Resource', axis=1)

    # Split the indeces (taxa namings) by the taxaSplit, grab the species names & remove "s__" for model database comparison
    species = renamed_dataframe.index.str.split(taxaSplit).str[-1].str.replace("s__", "", regex=False)

    # Find entries present in the model database
    present_mask = species.isin(updatedModelDatabase_df['Species'])
    present_df = renamed_dataframe.loc[present_mask]

    # Find entries absent in the model database
    absent_mask = ~present_mask
    absent_df = renamed_dataframe.loc[absent_mask]

    return present_df, absent_df


def split_taxonomic_groups(df, flagLoneSpecies=False, taxaSplit=';'):
    """
    Split the taxonomic groups in the index of the input DataFrame and create separate DataFrames for each taxonomic level.

    Args:
        merged_df (pd.DataFrame): The input DataFrame with taxonomic groups in the index, after model database mapping & normalization.

    Returns:
        dict: A dictionary with keys as taxonomic levels and values as the corresponding DataFrames.
    """
    logger = logging.getLogger('main.operations')
    logger.info("Splitting taxonomic levels in seperate dataframes, one dataframe per level.")

    # Define the mapping of prefixes to full taxonomic level names
    prefix_to_level = {
        'k': 'Kingdom',
        'p': 'Phylum',
        'c': 'Class',
        'o': 'Order',
        'f': 'Family',
        'g': 'Genus',
        's': 'Species'
    }
    
    # Split the index column into separate columns for each taxonomic level
    taxonomic_split_df = df.index.str.split(taxaSplit, expand=True)

    # Convert taxonomic_split_df to DataFrame
    taxonomic_split_df = taxonomic_split_df.to_frame()
    
    # Determine which levels are present in the data
    present_levels = []
    for col in taxonomic_split_df.columns:
        # Get the first character of the prefix (without "__")
        prefix = taxonomic_split_df[col].str[:1].iloc[0]
        if prefix in prefix_to_level:
            present_levels.append(prefix_to_level[prefix])
            # Remove the full prefix (with "__") from the data
            taxonomic_split_df[col] = taxonomic_split_df[col].str[3:]
    
    # Assign the determined levels as column names
    taxonomic_split_df.columns = present_levels
    
    # Replace NaN by empty strings
    taxonomic_split_df = taxonomic_split_df.fillna('')
    
    # Reset the index
    taxonomic_split_df = taxonomic_split_df.reset_index(drop=True)

    # Concatenate the taxonomic_split_df and the abundance data from taxonomic_levels_df
    taxonomic_levels_df = pd.concat([taxonomic_split_df, df.reset_index(drop=True)], axis=1)

    # Initialize a dictionary to store DataFrames for each taxonomic level
    taxonomic_dfs = {}

    # Iterate through the present taxonomic levels and create a DataFrame for each level
    for level in present_levels:
        level_df = taxonomic_levels_df[[level] + list(taxonomic_levels_df.columns[len(present_levels):])]
        level_df = level_df.rename(columns={level: 'Taxon'})

        # Set the 'Taxon' column as the index and remove rows with an empty string in the 'Taxon' column
        level_df = level_df.set_index('Taxon')
        level_df = level_df.loc[level_df.index != '']

        # Group by index and sum the rows with the same name (as especially taxa from higher taxonomic level might occour in multiple rows after seperation)
        level_df = level_df.groupby(level_df.index.name).sum()

        # Add the DataFrame to the dictionary
        taxonomic_dfs[level] = level_df

    return taxonomic_dfs


def calculate_metrics(dataframes_normalized, dataframes, group=None):
    """
    Calculate & then summarize alpha & beta diversity metrics, read counts & Firmicutes to Bacteroidetes ratio (for pyhlum)
    to compare pre-mapping & post-mapping status & evaluate the mapping coverage.

    Args:
        dataframes (dict): A dictionary with keys as taxonomic levels and values as the corresponding DataFrames.
        group (dict):      Stratified dictionary with keys as taxonomic levels and values as the corresponding DataFrames.

    Returns:
        dicts: Dictionaries with keys as taxonomic levels and values as the calculated metrics (metrics, abundance metrics, coverage summary statistics).
    """
    logger = logging.getLogger('main.operations')
    
    metrics = {}
    beta_diversity_metrics = {}
    abundance_metrics = {}
    summ_stats = {}

    for level in dataframes_normalized.keys():
        df = dataframes[level]
        df_normalized = dataframes_normalized[level]
        if group is not None:
            df_complete = dataframes[level]
            df_complete_normalized = dataframes_normalized[level]
            df = df_complete[group]
            df_normalized = df_complete_normalized[group]

        # Calculate read counts for dataframe
        read_counts_df = df.sum()
        mean_read_counts_df = np.mean(read_counts_df)
        std_read_counts_df = np.std(read_counts_df)
        
        # Get total number of taxa
        num_taxa = df_normalized.shape[0]

        # Estimate total number of named taxa, by checking for presence of "-"
        # or multiple uppercase letters in a row in taxa name
        taxa_name_contains_dash = df_normalized.index.str.contains('-')
        taxa_name_contains_uppercase = df_normalized.index.str.contains(r'[A-Z]{2,}')
        est_taxa_wo_standard_name = df_normalized[taxa_name_contains_dash | taxa_name_contains_uppercase]
        est_num_taxa_w_standard_name = len(df_normalized.index) - len(est_taxa_wo_standard_name)

        # Calculate non-zero entries per column to get species richness per sample
        species_richness = (df_normalized != 0).sum()
        mean_species_richness = np.mean(species_richness)
        std_species_richness = np.std(species_richness)

        # Calculate alpha diversity using the pielou's evenness
        shannon_index = -1 * (df_normalized * df_normalized.apply(np.log)).sum()
        pielous_evenness = shannon_index / species_richness.apply(np.log)
        mean_pielous_evenness = np.mean(pielous_evenness)
        std_pielous_evenness = np.std(pielous_evenness)
        
        level_metrics = {
            'read_counts': read_counts_df,
            'pielous_evenness': pielous_evenness,
        }

        # Calculate beta diversity using bray-curtis dissimilarity
        beta_diversity_metrics[level] = calc_bray_curtis_matrix(df_normalized)
        
        # Calculate mean + SD, min & max relative abundance of all taxa & in 
        # how many samples a taxa is present
        level_abundance_metrics = pd.DataFrame({
            'mean': df_normalized.mean(axis=1),
            'SD': df_normalized.std(axis=1),
            'minimum': df_normalized.min(axis=1),
            'maximum': df_normalized.max(axis=1),
            'non_zero_count': (df_normalized != 0).sum(axis=1)
        })
        
        abundance_metrics[level] = level_abundance_metrics

        if level == 'Phylum':
            # Calculate Firmicutes to Bacteroidetes ratio
            firmicutes = df_normalized.loc['Firmicutes'] if 'Firmicutes' in df_normalized.index else 0
            bacteroidetes = df_normalized.loc['Bacteroidetes'] if 'Bacteroidetes' in df_normalized.index else 0

            try:
                fb_ratio = firmicutes / bacteroidetes
            except ZeroDivisionError:
                fb_ratio = 0

            level_metrics.update({
                'firmicutes_bacteroidetes_ratio': fb_ratio,
            })

        # Add the metrics to the main dictionary
        metrics[level] = pd.DataFrame.from_dict(level_metrics)
        summ_stats[level] = [num_taxa, est_num_taxa_w_standard_name, \
                             mean_species_richness, std_species_richness, \
                             mean_pielous_evenness, std_pielous_evenness, \
                             mean_read_counts_df, std_read_counts_df]

    return metrics, abundance_metrics, beta_diversity_metrics, summ_stats

def calc_bray_curtis_matrix(rel_abundances):
    """
    Calculates Bray-Curtis dissimilarity matrix for all samples in a taxa abundances dataframe.

    Args:
        rel_abundances (pandas dataframe): A dataframe containing (relative) 
                                            abundances of taxa per samples, with
                                            columns = samples, rows = taxa.

    Returns:
        Pandas dataframe with Bray-Curtis dissimilarity matrix of pairwise sample comparisons.
    """
    logger = logging.getLogger('main.operations')
    
    samples = rel_abundances.columns
    num_samples = len(samples)
    
    # Convert DataFrame to numpy array for faster calculations
    abundance_array = rel_abundances.values
    
    # Pre-calculate sum of relative abundances for each sample
    sample_sums = abundance_array.sum(axis=0)
    
    # Pre-allocate the result matrix
    bray_curtis_matrix = np.zeros((num_samples, num_samples))
    
    # Calculate bray-curtis dissimilarity between all samples in pairwise manner & store in result matrix
    for i in range(num_samples):
        for j in range(i+1, num_samples):
            numerator = np.sum(np.abs(abundance_array[:, i] - abundance_array[:, j]))
            denominator = sample_sums[i] + sample_sums[j]
            dissimilarity = numerator / denominator
            bray_curtis_matrix[i, j] = dissimilarity
            bray_curtis_matrix[j, i] = dissimilarity
    
    return pd.DataFrame(bray_curtis_matrix, index=samples, columns=samples)

def process_microbial_abundances(input_file1, input_file2, output_path=None, cutoff=0.000001, output_format="csv", stratification_file=None, flagLoneSpecies=False, taxaSplit=";", removeCladeExtensionsFromTaxa=True, whichModelDatabase="full_db", userDatabase_path="", sample_read_counts_cutoff=1):
    # Initialize logger to generate a MARS log file    
    logger = setup_logger('main', os.path.join(output_path, 'MARS.log'))
    logger_taxa_below_cutoff = setup_logger('taxa_below_cutoff', os.path.join(output_path, 'MARS_taxaBelowCutoff.log'))
    
    logger.info(f'INPUT VARIABLES - taxaSplit: {taxaSplit}, flagLoneSpecies: {flagLoneSpecies}, sample_read_counts_cutoff: {sample_read_counts_cutoff}, cutoff: {cutoff}, removeCladeExtensionsFromTaxa: {removeCladeExtensionsFromTaxa}, whichModelDatabase: {whichModelDatabase}, userDatabase_path (if whichModelDatabase is set to "user_db"): {userDatabase_path}.')
    
    # Run MARS
    # Step 1: Check input data & preprocess
    [preprocessed_dataframe, dfvalues_are_rel_abundances] = load_input_and_preprocess(input_file1, input_file2, taxaSplit=taxaSplit, input_is_df=False)

    # Optional Step: Remove potential clade extensions (e.g. "clade A"; " A") from taxa namings if set true
    if removeCladeExtensionsFromTaxa == True:
        preprocessed_dataframe = remove_clades_from_taxaNames(preprocessed_dataframe, taxaSplit=taxaSplit)
    
    # Optional Step: Concatenate genus name with species epithet if both are present, otherwise leave species column unchanged
    if flagLoneSpecies:
        preprocessed_dataframe = concatenate_genus_and_species_names(preprocessed_dataframe, taxaSplit=taxaSplit)

    # Step 2: Rename taxa according to resources/renaming.json to share same nomenclature as the model-database
    renamed_dataframe = rename_taxa(preprocessed_dataframe)

    # Step 3: Filter out samples with too few total reads from subsequent analysis in case input 
    # is in absolute read counts (not relative abundance)
    if dfvalues_are_rel_abundances == False:
        renamed_dataframe = filter_samples_low_read_counts(renamed_dataframe, sample_read_counts_cutoff=sample_read_counts_cutoff)
    
    # Step 4: Normalize & apply taxa relative abundance cutoff to renamed_dataframe
    logger.info('Normalizing & applying the species relative abundance cutoff to the pre-mapping dataframe.')
    [dataframe_afterCutoff, normalized_dataframe_afterCutoff] = normalize_dataframe(renamed_dataframe, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=cutoff)

    # Step 5: Check for presence of input taxa in a specified model database (AGORA2, APOLLO, 
    # combination of both or user-defined)
    present_dataframe, absent_dataframe = check_presence_in_modelDatabase(dataframe_afterCutoff, whichModelDatabase=whichModelDatabase, userDatabase_path=userDatabase_path, taxaSplit=taxaSplit)
    if present_dataframe.empty:
        logger.critical('No species from the input data were found & could be mapped to the reconstruction database. Led to termination of code execution with ValueError.')
        raise ValueError("No species from the input data were found & could be mapped to the reconstruction database.")

    # Step 6.1: Normalize present_dataframe & absent_dataframe (with cutoff = 0,
    # because cutoff was already applied on dataframe_afterCutoff)
    logger.info('Normalizing the post-mapping dataframes (present & absent).')
    [present_dataframe_afterCutoff, normalized_present_dataframe] = normalize_dataframe(present_dataframe, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=0, dataframe_to_normalize_to=dataframe_afterCutoff)
    [absent_dataframe_afterCutoff, normalized_absent_dataframe] = normalize_dataframe(absent_dataframe, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=0, dataframe_to_normalize_to=dataframe_afterCutoff)
    # Step 6.2: Additionally normalize present_dataframe to its own total read count to be valid input for modelling 
    # (with relative abundances of present species per sample summing to 1)
    [_, normalized_present_dataframe_adjForModelling] = normalize_dataframe(present_dataframe, dfvalues_are_rel_abundances=dfvalues_are_rel_abundances, cutoff=0)

    # Step 7.1: Seperate normalized dataframes by taxonomic levels 
    # (one df per taxonomic level, for both normalized & not-normalized dataframes)
    preMapped_dataframes = split_taxonomic_groups(dataframe_afterCutoff, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    present_dataframes = split_taxonomic_groups(present_dataframe_afterCutoff, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    absent_dataframes = split_taxonomic_groups(absent_dataframe_afterCutoff, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    preMapped_dataframes_normalized = split_taxonomic_groups(normalized_dataframe_afterCutoff, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    present_dataframes_normalized = split_taxonomic_groups(normalized_present_dataframe, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    absent_dataframes_normalized = split_taxonomic_groups(normalized_absent_dataframe, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)
    present_dataframes_adjForModelling = split_taxonomic_groups(normalized_present_dataframe_adjForModelling, flagLoneSpecies=flagLoneSpecies, taxaSplit=taxaSplit)

    # Step 7.2: Add "pan" prefix to species names (in index) of present_dataframes_adjForModelling 
    # (because pan-species reconstructions will be used in MgPipe)
    present_dataframes_adjForModelling['Species'].index = 'pan' + present_dataframes_adjForModelling['Species'].index

    # Step 8.1: Calculate metrics on mapping coverage & microbiome composition
    logger.info('Calculating metrices for pre-mapping dataframes.')
    pre_mapping_metrics, pre_mapping_abundance_metrics, pre_mapping_beta_diversity, pre_mapping_summ_stats = calculate_metrics(preMapped_dataframes_normalized, preMapped_dataframes)
    logger.info('Calculating metrices for post-mapping present dataframes.')
    present_post_mapping_metrics, present_post_mapping_abundance_metrics, present_post_mapping_beta_diversity, present_post_mapping_summ_stats = calculate_metrics(present_dataframes_adjForModelling, present_dataframes)
    logger.info('Calculating metrices for post-mapping absent dataframes.')
    absent_post_mapping_metrics, absent_post_mapping_abundance_metrics, absent_post_mapping_beta_diversity, absent_post_mapping_summ_stats = calculate_metrics(absent_dataframes_normalized, absent_dataframes)
    
    # Step 8.2: Combine pre- and postMapping information of metrices, where needed
    logger.info('Combining metrics dataframes.')
    combined_metrics = combine_metrics(pre_mapping_metrics, present_post_mapping_metrics, df_type="metrics", dfvalues_are_rel_abundances=dfvalues_are_rel_abundances)
    combined_summ_stats = combine_metrics(pre_mapping_summ_stats, present_post_mapping_summ_stats, df_type="summ_stats", dfvalues_are_rel_abundances=dfvalues_are_rel_abundances)
    
    # Optional Step: If stratification groups are provided, stratify dataframe on these 
    # groups & calculate metrices on them too
    stratification_groups = {}
    stratification_groupnames = []
    if stratification_file is not None:
        stratification = pd.read_csv(stratification_file)
        for group in stratification["group"].unique():
            # Select the columns for this group
            group_columns = list(stratification[stratification["group"] == group]["samples"])
            pre_group_metrics, pre_group_abundance_metrics, pre_group_summ_stats = calculate_metrics(preMapped_dataframes_normalized, preMapped_dataframes, group=group_columns)
            post_group_metrics, post_group_abundance_metrics, post_group_summ_stats = calculate_metrics(present_dataframes_normalized, present_dataframes, group=group_columns)

            combined_group_metrics = combine_metrics(pre_group_metrics, post_group_metrics, df_type="metrics")
            combined_group_summ_stats = combine_metrics(pre_group_summ_stats, post_group_summ_stats, df_type="summ_stats")
            group_name = f"{group.lower()}_stratified_metrics"
            stratification_groupnames.append(group_name)
            stratification_groups[group_name] = [combined_group_metrics, combined_group_summ_stats, \
                                                pre_group_abundance_metrics, post_group_abundance_metrics]
        logger.info(f'Stratifying taxa dataframes using following groups: {stratification_groupnames}.')

    # Step 9: Store all result dataframes in a structure & save, if output-path is provided
    dataframe_groups = {'normalized_preMapped': preMapped_dataframes_normalized, 
                        'renormalized_mapped_forModelling': present_dataframes_adjForModelling,
                        'normalized_mapped': present_dataframes_normalized,
                        'normalized_unmapped': absent_dataframes_normalized,
                        'metrics': [combined_metrics, combined_summ_stats, \
                                    pre_mapping_abundance_metrics, \
                                   present_post_mapping_abundance_metrics, \
                                   absent_post_mapping_abundance_metrics, \
                                   pre_mapping_beta_diversity, \
                                   present_post_mapping_beta_diversity, \
                                   absent_post_mapping_beta_diversity]}
    
    dataframe_groups.update(stratification_groups)

    if output_path is not None:
        logger.info(f'Saving output to {output_path}.')
        save_dataframes(dataframe_groups, output_path, output_format)

        renamed_dataframe.to_csv(os.path.join(output_path, 'preprocessedInput_afterRenaming.csv'), sep=',')


    return dataframe_groups

def runMARS(counts: pd.DataFrame, taxa: pd.DataFrame, database: str = 'AGORA2-APOLLO', output_format: str = 'csv', taxa_delimiter: str = ';', 
            flag_lone_species: bool = False, cutoff: float = 0.000001, sample_read_counts_cutoff: int = 1, remove_clades: bool = True,
            output_dir: str = os.path.join('.','output','MARS')) -> (pd.DataFrame, pd.DataFrame):
    # Generate output folders
    os.makedirs(output_dir, exist_ok = True)
    
    # Rename database (if required)
    if database == 'AGORA2-APOLLO':
        database = 'full_db'
    
    # Build temporal file
    tmp_df = pd.concat([taxa, counts], axis = 1)
    tmp_df.to_csv('tmpAbundance.txt', sep = '\t', index = False)
    
    # Run MARS
    process_microbial_abundances('tmpAbundance.txt', None, output_path = output_dir, cutoff=cutoff, output_format=output_format, stratification_file=None,
                                 flagLoneSpecies=flag_lone_species, taxaSplit=taxa_delimiter, removeCladeExtensionsFromTaxa=remove_clades, whichModelDatabase = database, userDatabase_path="", sample_read_counts_cutoff=sample_read_counts_cutoff)
    
    # Remove temporal files
    os.remove('tmpAbundance.txt')
    
    # Load MARS result
    out_mars = pd.read_csv(os.path.join(output_dir, 'renormalized_mapped_forModelling','renormalized_mapped_forModelling_species.' + output_format))
    
    # Adapt MARS output to QIIME2
    countsNew = out_mars
    taxaNew = pd.DataFrame({'Taxon': countsNew.loc[:,'Taxon']})
    countsNew.drop('Taxon', inplace = True, axis = 1)
    
    countsNew.index = [f'OTU_{i+1}' for i in range(len(countsNew))]
    countsNew.index.name = '#OTU ID'
    taxaNew.index = [f'OTU_{i+1}' for i in range(len(countsNew))]
    taxaNew.index.name = 'Feature ID'
    
    # Output variable
    return (countsNew, taxaNew)
