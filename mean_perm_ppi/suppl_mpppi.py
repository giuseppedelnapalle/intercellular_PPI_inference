#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created in Jun 2020
@title: supplementary functions for mean permutation PPI inference algorithm
@description: 
    required for infer_PPI function
    version 1.7
@author: nikola
"""

"""
TO BE MODIFIED.
"""

import numpy as np
import pandas as pd
import itertools as it
import random


# SECTION_1 
# function from CSN module (user-defined)
def merge_GEP(*args,    # GEPs, one for each cell_type
              **kwargs    # keyword argument, cell_type = cell_type
              ) -> pd.DataFrame:
    """
    merge GEPs & append cell_type to column name
    """
    data = pd.concat(args, axis=1)
    
    # rename columns
    lst = list(kwargs.values())[0]
    k = len(lst)
    n = int(data.shape[1]/k)
    
    col_name = list()
    for i in range(k):
        for j in range(n):
            fld1 = lst[i]
            fld2 = data.columns[int(n*i+j)]
            col_name.append(fld1 + "_" + fld2)
            
    data.columns = col_name
    return data


# SECTION_2 
# functions for prefilters (method_2)

# 2.1 filter function 1
def filter_genes_by_value(GEPs: pd.DataFrame):
    """
    filter genes if all values are 0 or NaN
    """
    row_sum = GEPs.sum(axis = 'columns')
    GEPs_filtered = GEPs.loc[~(row_sum == 0)]
    GEPs_filtered = GEPs_filtered.dropna(axis=0, how='all')
    
    return GEPs_filtered


# 2.2 filter function 2
def filter_interactions_by_genes(interactions: pd.DataFrame, 
                                 GEPs: pd.DataFrame, 
                                 proteinComplex: pd.DataFrame,
                                 separators: list = ['_', ' ', ':'],
                                 suffixes: tuple = ('ligand', 'receptor'),
                                 interaction_name: str = 'interaction_name') -> pd.DataFrame:
    """
    remove interaction if one or two components not in gene list
    at least one gene out of many required to keep the component for proteinComplex
    """
    # obtain genes for each component
    isComplex = is_protein_complex(interactions, proteinComplex, suffixes=suffixes)
    # isSingleGene = is_single_gene(interactions, separators = separators, suffixes = suffixes)
    pd_upd = get_genes_for_component(isComplex, 
                                     interactions, proteinComplex, 
                                     suffixes = suffixes, interaction_name = interaction_name)
            
    # remove interaction if one of components missing
    filt = (pd_upd[suffixes[0]].str.len() == 0) | (pd_upd[suffixes[1]].str.len() == 0)
    pd_upd = pd_upd.loc[~filt]
            
    # remove interaction if all genes of one of components missing
    filt2 = (is_component_in_GEPs(pd_upd[suffixes[0]], GEPs)) & (is_component_in_GEPs(pd_upd[suffixes[1]], GEPs))
    pd_upd = pd_upd.loc[filt2]
    
    return pd_upd


# functions for filter_interactions_by_genes
# 2.2.1
def is_protein_complex(interactions: pd.DataFrame, 
                       proteinComplex: pd.DataFrame,
                       suffixes: tuple = ('ligand', 'receptor')) -> np.array:
    """
    inspect if component is protein complex
    component present in proteinComplex considered protein complex
    """
    for i in range(2):
        components = interactions[suffixes[i]]  # pd.Series
        col = components.isin(proteinComplex.index)
        
        if i == 0:
            isComplex = np.array(col)
        else:
            isComplex = np.vstack([isComplex, np.array(col)])
            
    return isComplex


# 2.2.2
def get_genes_for_component(isComplex: np.array, 
                            interactions: pd.DataFrame, 
                            proteinComplex: pd.DataFrame,
                            suffixes: tuple = ('ligand', 'receptor'),
                            interaction_name: str = 'interaction_name')-> pd.DataFrame:
    """
    find genes for each component 
    updated interactions with same gene_symbol as GEPs returned
    """
    for i in range(2):
        components = interactions[suffixes[i]]  # pd.Series
        genes = ['missing' for j in range(interactions.shape[0])]
        
        # genes of protein complex
        m = get_genes_multiple(isComplex[i], components, proteinComplex)
        # update genes
        genes = update_genes(genes, isComplex[i], components, m)
        
        if i == 0:
            dict_genes = {suffixes[i]: genes}
        else:
            dict_genes[suffixes[i]] = genes
            
    pd_genes = pd.DataFrame(dict_genes)
    pd_genes.index = interactions[interaction_name]
    
    return pd_genes


# 2.2.2.1
def update_genes(genes: list, 
                 isComplex_arr: np.array,  # NOT isComplex[i] (single array of isComplex)
                 components: pd.Series, 
                 genes_mult: list) -> list:
    """
    update genes for the component according to if the component composed of single gene
    """
    lst_compnt = list(components)
    
    for i, item in enumerate(genes):
        if isComplex_arr[i] == False:
            genes[i] = list(lst_compnt[i].split())
        else:
            genes[i] = genes_mult[i]
    
    return genes


# 2.2.2.2
def get_genes_multiple(isComplex_arr: np.array,  # isComplex[i]
                       components: pd.Series,  # interactions[suffixes[i]]
                       proteinComplex: pd.DataFrame) -> list:
    
    """
    find genes for proteinComplex for each pd.Series of components
    """
    genes = ['missing' for i in range(len(isComplex_arr))]
    
    for i in range(len(components)):
        if isComplex_arr[i] == True:
            gene_list = proteinComplex.loc[proteinComplex.index == components[i]].values.tolist()
            # gene_list = proteinComplex.loc[proteinComplex.index.isin(components[i])].values.tolist()  # equivalent
            gene_list = [item for elem in gene_list for item in elem]
            gene_list = [x for x in gene_list if str(x) != 'nan']
            genes[i] = gene_list
    
    return genes


# 2.2.3
def is_component_in_GEPs(components: pd.Series, 
                         GEPs: pd.DataFrame) -> pd.Series:
    """
    determine if genes of component in GEPs
    at least one gene of component required
    """
    gene_list = GEPs.index
    IsInGEPs = components.apply(lambda x: is_in_geneList(x, gene_list))
    
    return IsInGEPs
    

# 2.2.3.1
def is_in_geneList(x: list, 
                   gene_list: pd.Index) -> bool:
    """
    check if genes in specified gene list
    """
    IsInGenelist = False
    for i in range(len(x)):
        IsInGenelist = IsInGenelist or (x[i] in gene_list)
    
    return IsInGenelist


# 2.3 filter function 3
def filter_genes_by_interactions(GEPs: pd.DataFrame, 
                                 interactions: pd.DataFrame) -> pd.DataFrame:
    """
    filter genes in GEPs by interactions
    """
    gene_list = interactions.values.tolist()
    
    genes_flat = [item for elem in gene_list for item in elem]
    genes_flat = list(set([item for elem in genes_flat for item in elem]))
    
    filt = [i for i in range(GEPs.shape[0]) if GEPs.index[i] in genes_flat]
    GEPs_filtered = GEPs.iloc[filt]
    
    return GEPs_filtered


# SECTION_3 
# functions for infer_PPI (method_3)

def get_cellType_combinations(cell_type: list) -> list:
    """
    calculate and sort combinations including itself
    """
    res = sorted(it.product(cell_type, repeat = 2))
    return res
  
      
def get_samplePair_index(GEPs: pd.DataFrame,
                         cell_type: list,
                         sample_index: int,  # single number
                         cell_type_comb: tuple) -> list:
    """
    obtain index for sample pair
    """
    num_sample = int(GEPs.shape[1] / len(cell_type))
    try:
        t = 1 if sample_index <= num_sample else 0
        test = 1 / t == 1
    except ZeroDivisionError:
        test = False
        print('Sample_index must be smaller than number of samples in GEPs.')
    
    if test:
        idx = []
        for i in range(2):
            rank = [j for j in range(len(cell_type)) if cell_type[j] == cell_type_comb[i]]
            rank = rank[0]
            x = sample_index + num_sample * rank
            idx.append(x)
        return idx
    else:
        return None


def _compute_mean_statistic(GEPs: pd.DataFrame, 
                           interaction_index: pd.Series,  # row of df returned by get_interaction_index() func, e.g. df.iloc[1]
                           sample_index: list) -> float:  # optimised
    """
    calculate the mean statistic for each interaction, 
    sample, cell_type_pair combination
    """
    # mean_statistics for both components
    mean_comp = [
        GEPs.iloc[interaction_index[i], sample_index[i]].mean() if 
        
        # count number of missing genes
        [j for j in range(len(interaction_index[i])) if interaction_index[i][j]==[]] == [] else 
        
        # add 0 before computing mean
        pd.concat([GEPs.iloc[list(filter(None, interaction_index[i])), sample_index[i]],
                   
                  # number of 0 = number of missing genes
                  pd.Series([0 for k in range(len(
                      [j for j in range(len(interaction_index[i])) if interaction_index[i][j]==[]]
                      ))])]).mean()
        
        for i in range(2)
        ]  # list comprehension

    return pd.Series(mean_comp).mean() if [i for i in range(len(mean_comp)) if mean_comp[i]==[]] == [] else np.nan  # conditional expression


def _mean_statistic_itparam(GEPs: pd.DataFrame, 
                           interaction_index: pd.DataFrame,  # df returned by get_interaction_index() func
                           sample_index: pd.Series,  # array returned by build_samplePairIndex_array() func
                           cell_type_comb: list,
                           param: list  # list of 3 indices for interaction_index, sample_index, cell_type_comb
                           ) -> float:  # optimised
    """
    calculate the mean statistic for each interaction, sample, cell_type_pair combination
    modified for iterable parameters
    """
    iidx, sidx, cidx = param[0], param[1], param[2]
    
    # mean_statistics for both components
    mean_comp = [
        GEPs.iloc[interaction_index.iloc[iidx][i], sample_index.iloc[sidx*len(cell_type_comb)+cidx][i]].mean() if 
        
        # count number of missing genes
        [j for j in range(len(interaction_index.iloc[iidx][i])) if interaction_index.iloc[iidx][i][j]==[]] == [] else 
        
        # add 0 before computing mean
        pd.concat([GEPs.iloc[list(filter(None, interaction_index.iloc[iidx][i])), sample_index.iloc[sidx*len(cell_type_comb)+cidx][i]],
                  # number of 0 = number of missing genes
                  pd.Series([0 for k in range(len(
                      [j for j in range(len(interaction_index.iloc[iidx][i])) if interaction_index.iloc[iidx][i][j]==[]]
                      ))])]).mean()
        
        for i in range(2)
        ]  # list comprehension

    return pd.Series(mean_comp).mean() if [i for i in range(len(mean_comp)) if mean_comp[i]==[]] == [] else np.nan


def get_interaction_index(GEPs: pd.DataFrame,  # filtered GEPs
                          interactions: pd.DataFrame  # filtered interactions
                          ) -> pd.DataFrame:
    """
    obtain gene index from GEPs for each component of interaction
    each element is a list in returned pd.dataframe
    if gene missing in GEPs, value=[]
    """
    num_inter = interactions.shape[0]
    num_gene = GEPs.shape[0]
    
    # loop over 2 components
    for i in range(2):
        inter_idx = [[] for i in range(num_inter)]
        
        # loop over num_inter interactions
        for j in range(num_inter):
            
            # loop over numGinC genes of specified component
            for k in range(len(interactions.iloc[j, i])):
                x = [u for u in range(num_gene) if GEPs.index[u] == interactions.iloc[j, i][k]]
                
                if x != []:
                    x = x[0]
                
                inter_idx[j].append(x)
                    
        if i == 0:
            dict_index = {interactions.columns[i]: inter_idx}
        else:
            dict_index[interactions.columns[i]] = inter_idx
            
    pd_index = pd.DataFrame(dict_index)
    pd_index.index = interactions.index
    
    return pd_index


def _mean_stat_shuffled_meta_itparam(GEPs: pd.DataFrame, 
                                     interaction_index: pd.DataFrame,  # df returned by get_interaction_index function
                                     param: list  # list containing 1. index for interaction_index & 2. iterations
			                         ) -> float:
    """
    calculate the mean statistic for shuffled data
    modified for iterable parameters
    """
    index = param[0]
    
    # obtain shuffled sample index, which is
    # equivalent to randomly choosing 2 items from sequence without replacement
    shuf_idx = random.sample(range(GEPs.shape[1]), 2)

    # mean_statistics for both components
    mean_comp = [
        GEPs.iloc[interaction_index.iloc[index][i], shuf_idx[i]].mean() if 
        
        # count number of missing genes
        [j for j in range(len(interaction_index.iloc[index][i])) if interaction_index.iloc[index][i][j]==[]] == [] else 
        
        # add 0 before computing mean
        pd.concat([GEPs.iloc[list(filter(None, interaction_index.iloc[index][i])), shuf_idx[i]],
                  
                  # number of 0 = number of missing genes
                  pd.Series([0 for k in range(len(
                      [j for j in range(len(interaction_index.iloc[index][i])) if interaction_index.iloc[index][i][j]==[]]
                      ))])]).mean()
        
        for i in range(2)
        ]  # list comprehension

    return pd.Series(mean_comp).mean() if [i for i in range(len(mean_comp)) if mean_comp[i]==[]] == [] else np.nan


def get_sample_id(GEPs: pd.DataFrame,
                  cell_type: list,
                  separator: str = '_') -> list:
    """
    obtain sample ID
    """
    num_samp = int(GEPs.shape[1] / len(cell_type))
    df = pd.DataFrame(list(map(lambda x: x.split(separator), GEPs.columns[:num_samp])))
    
    return list(df.iloc[:,1])


def build_result_df(GEPs: pd.DataFrame, 
                    interactions: pd.DataFrame,
                    sample_id: list,
                    cell_type_comb: list,
                    separator: str = '_') -> pd.DataFrame:
    """
    builds an empty data frame for results
    rows: interactions
    columns: sample*cell_type_comb
    """
    columns = []
    for sample in sample_id:
        for comb in cell_type_comb:
            columns.append('{}{}{}'.format(sample, separator, comb))
    
    result = pd.DataFrame(index = interactions.index, columns = columns, dtype = 'float32')
    return result


def build_samplePairIndex_array(GEPs: pd.DataFrame,
                                cell_type: list,
                                sample_id: list,
                                cell_type_comb: list,
                                separator: str = '_') -> pd.Series:
    """
    builds an empty array for sample_pair_index
    sample*cell_type_comb
    """
    idx = ['{}{}{}'.format(sample, separator, comb) for 
           sample in sample_id for comb in cell_type_comb]
    
    arr = [get_samplePair_index(GEPs, cell_type, i, cell_type_comb[j]) for 
           i in range(len(sample_id)) for j in range(len(cell_type_comb))]
        
    return pd.Series(arr, index=idx)


def _compute_pvalue_itparam(stat_obs: np.array,  # statistic of observed data
                            stat_shuf: np.array,  # statistic of shuffled data
                            iterations: int,
                            param: list):  # list containing index for row & col of stat_obs
    """
    compute p-value based on observed statistics & null distribution 
    derived from shuffled data
    """
    # choose observed mean_statistic by parameter list
    row, col = param[0], param[1]
    x = stat_obs[row, col]
    
    # p-value calculated based on proportion of the statistics which are as large as or 
    # larger than the actual statistics
    pvalue = 1.0 if str(x) == 'nan' else (stat_shuf.iloc[row] >= x).sum() / iterations
    
    # pvalue = 1.0 if str(x) == 'nan' else round(len([i for i in range(iterations) if stat_shuf[row, i] >= x]) / 
    #           iterations, 4)
    
    return pvalue


def update_pvalue(df_result: pd.DataFrame,
                  pvalue: np.array):
    """
    fill p-value to df_result
    """
    for i in range(df_result.shape[0]):
        for j in range(df_result.shape[1]):
            df_result.iloc[i, j] = pvalue[i, j]
            
    return df_result


def get_PPI(df_pvalue: pd.DataFrame,
            alpha : float = .05) -> pd.DataFrame:
    """
    transform df_prob to data frame of PPI, '1' means PPI exists, 
    '0' means PPI does not exist, at significance level of alpha.
    """
    df = df_pvalue.apply(lambda arr: is_smaller(arr, alpha), axis=0)
    return df


def is_smaller(arr: pd.Series, 
              alpha: float) -> pd.Series:
    """
    determine if item in array is <= q,
    if yes, assign '1', else assign '0'.
    """
    res = arr.apply(lambda x: x < alpha).multiply(1)
    return res


def sig_mean_statistic(mean_stat: np.array,
                       df_pvalue: pd.DataFrame,
                       alpha: float = .05) -> pd.DataFrame:
    """
    compute significant mean_statistic
    If p.value < alpha, the value will be the mean. Alternatively, the value is set to 0.
    """
    sig_stat = mean_stat.copy()
    
    for i in range(sig_stat.shape[0]):
        for j in range(sig_stat.shape[1]):
            sig_stat[i, j] = sig_stat[i, j] if df_pvalue.iloc[i, j] < alpha else 0
    
    sig_stat = pd.DataFrame(sig_stat)
    sig_stat.index = df_pvalue.index
    sig_stat.columns = df_pvalue.columns
    
    return sig_stat
