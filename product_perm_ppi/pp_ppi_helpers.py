#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created in Jun 2020
@title: supplementary functions for lr_product PPI inference algorithm
@description: 
    required for preprocess & infer_PPI methods
    version 1.2 (same as v1.1)
@author: nikola
"""

import numpy as np
import pandas as pd
import itertools as it
import random
import re
from numba import jit  #, prange


# SECTION_1 
# function from CSN module (user-defined)
def merge_geps(*args,    # geps, one for each cell_type
              **kwargs    # keyword argument, cell_type = cell_type
              ) -> pd.DataFrame:
    """
    merge geps & append cell_type to column name
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
# function for __init__ method
def strip_prefix(x: str, 
                 char: str) -> str:
    """
    remove prefix in sample_ID
    """
    res = re.sub(char, "", x)
    return res


# SECTION_3
# functions for preprocess (method_2)

# 2.1 filter function 1
def filter_genes_by_value(geps: pd.DataFrame):
    """
    filter genes if all values are 0 or 1 or NaN across columns
    
    The "1" values in the expression matrix txt files are genes with insufficient
    evidence of expression (these genes are either not expressed or have inadequate statistical power to be imputed).
    
    The NA values are genes that have inadequate statistical power to be imputed.
    """
    row_sum = geps.sum(axis=1)  # or axis='columns'
    geps_filt = geps.loc[~(row_sum == 0)]
    
    filt = ~geps_filt.apply(is_all_one, axis=1)
    geps_filt = geps_filt.loc[filt]
    
    geps_filt = geps_filt.dropna(axis=0, how='all')
    
    return geps_filt


# functions for filter_genes_by_value
# 2.1.1
def is_all_one(arr: pd.Series) -> pd.Series:
    """
    determine if all values are 1 in array
    """
    y = arr.apply(lambda x: x == 1)
    res = len(y) == y.sum()
    
    return res


# 2.2 filter function 2
def filter_interactions_by_genes(interactions: pd.DataFrame, 
                                 geps: pd.DataFrame) -> pd.DataFrame:
    """
    remove interaction if one or two components not in gene list
    """
    df = pd.DataFrame(index=interactions.index, columns=interactions.columns, dtype=bool)
    
    for i in range(interactions.shape[0]):
        for j in range(2):
            df.iloc[i, j] = True if interactions.iloc[i, j] in geps.index else False
    
    agg = df.multiply(1).sum(axis=1)  # row sum
    filt = agg == 2
    inter_filt = interactions.loc[filt]
    
    return inter_filt
    

# 2.3 filter function 3
def filter_genes_by_interactions(geps: pd.DataFrame, 
                                 interactions: pd.DataFrame) -> pd.DataFrame:
    """
    filter genes in geps by interactions
    """
    gene_list = interactions.values.tolist()
    genes_flat = list(set([item for elem in gene_list for item in elem]))
    
    filt = [i for i in range(geps.shape[0]) if geps.index[i] in genes_flat]
    geps_filt = geps.iloc[filt]
    
    return geps_filt


# 2.4 transform values of 1 in geps
@jit(nopython=True)
def convert_one(arr: np.array,
                value_thresh: float = .01) -> np.array:
    """
    transform values of 1 in geps
    """
    # convert value of 1 to small number
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            if arr[i, j] == 1:
                arr[i, j] = value_thresh
        
    return arr


# SECTION_4
# functions for compute_lr_product, compute_pvalue, significant_lr_product (method 3_5)

def get_interaction_index(geps: pd.DataFrame,  # filtered geps
                          interactions: pd.DataFrame  # filtered interactions
                          ) -> pd.DataFrame:
    """
    obtain index for interaction from geps
    """
    num_inter = interactions.shape[0]
    num_g = geps.shape[0]
    df = pd.DataFrame(index=interactions.index, columns=interactions.columns)
    
    for i in range(num_inter):
        for j in range(2):
            df.iloc[i, j] = [k for k in range(num_g) if geps.index[k] == interactions.iloc[i, j]][0]
    
    return df


def get_cell_type_combinations(cell_type: list) -> list:
    """
    calculate and sort combinations including itself
    """
    ct_comob = sorted(it.product(cell_type, repeat = 2))
    return ct_comob


def get_sample_pair_index(cell_type: list,
                          exc_cell_type_ctrl: str,
                          sample_index: int,  # single number
                          cell_type_comb: tuple,
                          size_cell_type: list,
                          test_index_mapping: pd.DataFrame) -> list:
    """
    obtain index for sample pair
    """    
    idx = []
    for i in range(2):
        rank = [j for j in range(len(cell_type)) if cell_type[j] == cell_type_comb[i]][0]
        
        if cell_type_comb[i] == exc_cell_type_ctrl:  # if no sample of exc_cell_type removed, x equivalent in 'if' & 'else' section
            x = test_index_mapping.loc[sample_index, 0] + sum(size_cell_type[:rank])  # sum(size_cell_type[:0])==0
            idx.append(x)
        else:
            x = sample_index + sum(size_cell_type[:rank])
            idx.append(x)
    
    return idx


def get_cell_type_size(cell_type: list,
                       sample_type: pd.DataFrame,
                       test_group: str = 'tumour',
                       exc_cell_type_ctrl: str = 'no_cell_type_excluded',) -> list:
    """
    obtain sample size of each cell type in filtered geps
    """
    num_samp = sample_type.shape[0]
    
    # test if argument exc_cell_type_ctrl given
    test = exc_cell_type_ctrl in cell_type
    
    if test:
        num_test = sample_type.loc[sample_type.iloc[:, 0] == test_group].shape[0]
        r_exc_ct = [i for i in range(len(cell_type)) if cell_type[i] == exc_cell_type_ctrl][0]
    
        size_ct = []
        for i in range(len(cell_type)):
            x = num_samp if i != r_exc_ct else num_test
            size_ct.append(x)
    else:
        size_ct = []
        for i in range(len(cell_type)):
            size_ct.append(num_samp)
            
    return size_ct


def sample_pair_index_array(sample_id: pd.Series,
                            cell_type: list,
                            sample_type: pd.DataFrame,
                            test_group: str = 'tumour',
                            exc_cell_type_ctrl: str = 'no_cell_type_excluded',
                            separator: str = '_') -> pd.DataFrame:
    """
    build an array for sample_pair_index
    sample*cell_type_comb
    """
    # cell type combinations of case group
    cell_type_comb = get_cell_type_combinations(cell_type)
    
    # cell type combinations of control group
    ct_ctrl = cell_type.copy()
    if exc_cell_type_ctrl in cell_type:
        ct_ctrl.remove(exc_cell_type_ctrl)
        
    cell_type_comb_ctrl = get_cell_type_combinations(ct_ctrl)
    
    # define parameters for get_sample_pair_index
    # param_1 size_cell_type
    size_ct = get_cell_type_size(cell_type, sample_type, test_group, exc_cell_type_ctrl)
    
    # param_2 test_index_mapping
    num_test = sample_type.loc[sample_type.iloc[:, 0] == test_group].shape[0]
    test_idx_map = pd.DataFrame([i for i in range(num_test)])
    test_idx_map.index = sample_type.loc[sample_type.iloc[:, 0] == test_group].index
    
    # search for sample_pair_index
    idx = []
    arr = []
    for i in range(len(sample_id)):
        if sample_type.iloc[i, 0] == test_group:
            for comb in cell_type_comb:
                x = '{}{}{}'.format(sample_id[i], separator, comb)
                idx.append(x)
                
                y = get_sample_pair_index(cell_type, exc_cell_type_ctrl, i, comb, size_ct, test_idx_map)
                arr.append(y)
        else:
            for comb in cell_type_comb_ctrl:
                x = '{}{}{}'.format(sample_id[i], separator, comb)
                idx.append(x)
                
                y = get_sample_pair_index(cell_type, exc_cell_type_ctrl, i, comb, size_ct, test_idx_map)
                arr.append(y)  
    
    return pd.DataFrame(arr, index=idx)


# note Pandas not supported by Numba
# @jit(nopython=True, parallel=True)  # only one advanced index is allowed. failed to implement numba execution on the function.
# main function
def compute_lr_product(arr_geps: np.array,
                       arr_interaction_idx: np.array,
                       arr_sample_pair_idx: np.array) -> np.array:
    """
    calculate the lr_product statistic for each interaction, sample, cell_type_pair combination
    lr_product defined as sqrt(l*r), 
    where l, r denote expresson of ligand and receptor
    """
    n_row = arr_interaction_idx.shape[0]
    n_col = arr_sample_pair_idx.shape[0]
    
    l_idx = arr_sample_pair_idx[:, 0]
    r_idx = arr_sample_pair_idx[:, 1]
    
    arr_res = np.zeros([n_row, n_col], dtype=np.float64)
    
    for i in range(n_row):
        # lig = arr_geps[arr_interaction_idx[i, 0], l_idx]  # equivalent
        # rec = arr_geps[arr_interaction_idx[i, 1], r_idx]
        # arr_res[i] = np.sqrt(lig*rec)
        arr_res[i] = np.sqrt(arr_geps[arr_interaction_idx[i, 0], l_idx] * arr_geps[arr_interaction_idx[i, 1], r_idx])
    
    return arr_res


def sample_pair_index_sampling(index_list: list,
                               iterations: int = 1000,
                               seed_num: int = 123) -> np.array:
    """
    sample sample_pair_index from input list
    """
    random.seed(seed_num)
    sample = random.choices(index_list, k=iterations*2)  # sampling with replacement
    arr = np.array([sample[:iterations], sample[iterations:iterations*2]]).reshape([iterations, 2])
    
    return arr


def lr_product_shuffled(arr_geps: np.array,
                        arr_interaction_idx: np.array,
                        cell_type: list,
                        sample_type: pd.DataFrame,
                        null_dist_ctrl: bool = True,  # valid when exc_cell_type_ctrl not given
                        test_group: str = 'tumour',
                        exc_cell_type_ctrl: str = 'no_cell_type_excluded',
                        exact: bool = False,  # exact or approximate method
                        iterations: int = 1000,
                        seed_num: int = 123) -> np.array:
    """
    compute lr_product statistic on shuffled data of control group.
    if exc_cell_type_ctrl provided, null distribution estimated from control samples, with exact or approximate method;
    else whether null distribution estimated from control samples can be specified with null_dist_ctrl argument,
    along with selection of exact or approximate method.
    """
    print('If argument exc_cell_type_ctrl provided, null distribution of lr_product estimated from control group, \
          overriding null_dist_ctrl argument.')
    
    # check if arguments consistent
    if null_dist_ctrl:
        # test if control group present
        test = sample_type.loc[sample_type.iloc[:,0] != test_group].shape[0] == 0
        if test:
            print('No control label found in sample_type. null_dist_ctrl set to False.')
            null_dist_ctrl = False
    
    # common parameters
    # param_1 size_cell_type
    size_ct = get_cell_type_size(cell_type, sample_type, test_group, exc_cell_type_ctrl)
    # param_2 num of row of result array
    n_row = arr_interaction_idx.shape[0]
    
    # test if argument exc_cell_type_ctrl given
    test2 = exc_cell_type_ctrl in cell_type
    
    if test2:  # defined cell_type removed in control samples
        # rank of cell_type excluded for control sample
        r_exc_ct = [i for i in range(len(cell_type)) if cell_type[i] == exc_cell_type_ctrl][0]
        
        # indices of control samples for geps
        ctrl_idx0 = sample_type.loc[~(sample_type.iloc[:, 0] == test_group)].index
        ctrl_idx0 = np.array([i for i in ctrl_idx0])  # converted to np.array to use broadcasting
        
        ct_i = [i for i in range(len(cell_type))]
        ct_i.remove(r_exc_ct)
        
        ctrl_idx = []
        for i in ct_i:
            rank = [j for j in range(len(cell_type)) if cell_type[j] == cell_type[i]][0]
            x = [ctrl_idx0 + sum(size_ct[:rank])]
            ctrl_idx.extend(x)  # list of np.array
        ctrl_idx = [item for elem in ctrl_idx for item in elem]  # expand np.array to flat list
        
        # compute lr_product
        if exact:
            if len(ctrl_idx) > 500:
                print('WARNING: number of samples for exact test over 500. Calculation may take a long time.')
            
            arr_sample_pair_idx = np.array(sorted(it.product(ctrl_idx, repeat = 2)))
            l_idx = arr_sample_pair_idx[:, 0]
            r_idx = arr_sample_pair_idx[:, 1]
            n_col = arr_sample_pair_idx.shape[0]
            arr_res = np.zeros([n_row, n_col], dtype=np.float64)
            
            for i in range(n_row):
                arr_res[i] = np.sqrt(arr_geps[arr_interaction_idx[i, 0], l_idx] * arr_geps[arr_interaction_idx[i, 1], r_idx])
        else:
            # sample sample_pair_index from all permutations
            arr_sp_idx_samp = sample_pair_index_sampling(ctrl_idx, iterations, seed_num)
            l_idx = arr_sp_idx_samp[:, 0]
            r_idx = arr_sp_idx_samp[:, 1]
            n_col = arr_sp_idx_samp.shape[0]
            arr_res = np.zeros([n_row, n_col], dtype=np.float64)
            
            for i in range(n_row):
                arr_res[i] = np.sqrt(arr_geps[arr_interaction_idx[i, 0], l_idx] * arr_geps[arr_interaction_idx[i, 1], r_idx])
            
    else:  # test2==False, no cell_type removed in any samples
        if null_dist_ctrl:
            # indices of control samples for geps
            ctrl_idx0 = sample_type.loc[~(sample_type.iloc[:, 0] == test_group)].index
            ctrl_idx0 = np.array([i for i in ctrl_idx0])  # converted to np.array to use broadcasting
            
            ct_i = [i for i in range(len(cell_type))]
            
            ctrl_idx = []
            for i in ct_i:
                rank = [j for j in range(len(cell_type)) if cell_type[j] == cell_type[i]][0]
                x = [ctrl_idx0 + sum(size_ct[:rank])]
                ctrl_idx.extend(x)  # list of np.array
            ctrl_idx = [item for elem in ctrl_idx for item in elem]  # expand np.array to flat list
            
            if exact:
                if len(ctrl_idx) > 500:
                    print('WARNING: number of samples for exact test over 500. Calculation may take a long time.')
                
                arr_sample_pair_idx = np.array(sorted(it.product(ctrl_idx, repeat = 2)))
                l_idx = arr_sample_pair_idx[:, 0]
                r_idx = arr_sample_pair_idx[:, 1]
                n_col = arr_sample_pair_idx.shape[0]
                arr_res = np.zeros([n_row, n_col], dtype=np.float64)
                
                for i in range(n_row):
                    arr_res[i] = np.sqrt(arr_geps[arr_interaction_idx[i, 0], l_idx] * arr_geps[arr_interaction_idx[i, 1], r_idx])
            else:  # exact==False
                # sample sample_pair_index from all permutations
                arr_sp_idx_samp = sample_pair_index_sampling(ctrl_idx, iterations, seed_num)
                l_idx = arr_sp_idx_samp[:, 0]
                r_idx = arr_sp_idx_samp[:, 1]
                n_col = arr_sp_idx_samp.shape[0]
                arr_res = np.zeros([n_row, n_col], dtype=np.float64)
                
                for i in range(n_row):
                    arr_res[i] = np.sqrt(arr_geps[arr_interaction_idx[i, 0], l_idx] * arr_geps[arr_interaction_idx[i, 1], r_idx])
        
        else:  # null_dist_ctrl==False
            arr_sp_idx_samp = sample_pair_index_sampling([i for i in range(arr_geps.shape[1])], iterations, seed_num)
            l_idx = arr_sp_idx_samp[:, 0]
            r_idx = arr_sp_idx_samp[:, 1]
            n_col = arr_sp_idx_samp.shape[0]
            arr_res = np.zeros([n_row, n_col], dtype=np.float64)
            
            for i in range(n_row):
                arr_res[i] = np.sqrt(arr_geps[arr_interaction_idx[i, 0], l_idx] * arr_geps[arr_interaction_idx[i, 1], r_idx])
        
    return arr_res


def _compute_pvalue_itparam(arr_stat: np.array, 
                            arr_stat_shuf: np.array,
                            iterations: int,
                            param: list  # list containing index for row & col of arr_lr_product
                            ) -> float:
    """
    compute p-values for lr_product based on estimated null distribution with multiprocessing framework
    p-value calculated based on proportion of the statistics which are as large as or larger than the actual statistics
    """
    row, col = param[0], param[1]
    
    x = arr_stat[row, col]
    pvalue = (arr_stat_shuf[row] >= x).sum() / iterations
    
    return pvalue
 

def _compute_pvalue_itparam_2(arr_stat: np.array, 
                              arr_stat_shuf: np.array,
                              iterations: int,
                              param: list  # list containing index for row & col of arr_lr_product
                              ) -> float:
    """
    compute p-values for lr_product for arr_stat with NaN value
    """
    row, col = param[0], param[1]
    
    x = arr_stat[row, col]
    pvalue = 1.0 if str(x) == 'nan' else (arr_stat_shuf[row] >= x).sum() / iterations
    
    return pvalue
    

def _sig_lr_product_itparam(arr_stat: np.array,
                            arr_pvalue: np.array,
                            alpha: float,
                            value_thresh: float,
                            row_index: int) -> np.array:
    """
    compute significant lr_product with multiprocessing framework
    If p.value < alpha, the value will be the lr_product. Alternatively, the value is set to value_thresh.
    """
    sig_stat = arr_stat[row_index]
    pvalue = arr_pvalue[row_index]
    
    filt = pvalue >= alpha
    sig_stat[filt] = np.full([filt.sum()], value_thresh)
    
    return sig_stat


def get_ppi(df_pvalue: pd.DataFrame,
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
