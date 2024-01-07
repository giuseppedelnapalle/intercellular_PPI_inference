#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created in Jun 2020
@title: infering intercellular PPI based on permutation test
@reference: CellPhoneDB https://doi.org/10.1038/s41596-020-0292-x
@description:
    algorithm for infering intercellular PPI [optimised version].
    mean of ligand and receptor expression in each cell type as statistic,
    permutation test used to test significicance of enrichment of PPI pairs.
    required functions imported from suppl_mpppi module
    version 1.7 (compatible with suppl_mpppi_v1.7)
    $input 
        1. cell type expression profiles (pd_dataframe; derived form CIBERSORTx)
        2. protein pairs (pd.dataframe)
        3. protein complex (pd.dataframe)
        4. cell types
    $output 
        1. data frame of probabilty of interactions between two proteins 
           from specified cell type pairs (pd.dataframe)
        2. data frame qualitative value of interactions between two proteins 
           from specified cell type pairs (pd.dataframe)
@author: nikola
"""

"""
TO BE MODIFIED.
"""


import numpy as np
import pandas as pd
import itertools as it
from functools import partial
from multiprocessing.pool import Pool
import suppl_mpppi as sp


# define class for the algorithm
# the meanPermPPI object is created from the merged gene expression profile (GEP) matrix derived from CIBERSORTx.
# GEP object returned by merge_GEP of CSN module (pd.dataframe).
# genes in rows and samples in columns (#sample=n). first n columns for cell_type_1, second n columns for cell_type_2, etc.
# the class provides methods for intercellular PPI inference.
class MeanPermPPI:
    
    # method 1
    def __init__(self, 
                 GEPs: pd.DataFrame, 
                 cell_type: list, 
                 interactions: pd.DataFrame, 
                 proteinComplex: pd.DataFrame,
                 sample_type: pd.DataFrame,
                 exc_cellType_ctrl: str):
        """
        create a new MeanPermPPI instance
        remove defined cell_type, e.g. 'malignant', in control samples
        """          
        # test if number of samples in mixture data & sample_type equal
        test = GEPs.shape[1]/len(cell_type) == sample_type.shape[0]
        
        if test:
            order_ct = [i for i in range(len(cell_type)) if cell_type[i] == exc_cellType_ctrl]
            order_ct = order_ct[0]
            
            # column index of to_be_removed_columns
            idx_rm = [i for i in range(sample_type.shape[0]) if sample_type.iloc[i,:].values == 'control']
            idx_rm = [order_ct * sample_type.shape[0] + i for i in idx_rm]
            # column names
            filt = GEPs.columns[idx_rm]
            
            # remove defined cell_type in control samples
            GEPs_filt = GEPs.loc[:, ~GEPs.columns.isin(filt)]
            
            # create new instance
            self.GEPs = GEPs_filt
            self.cell_type = cell_type  
            self.numCellTypes = len(cell_type)
            self.interactions = interactions
            self.proteinComplex = proteinComplex
            self.sample_type = sample_type,
            self.exc_cellType_ctrl = exc_cellType_ctrl
            
        else:
            print('Wrong number of rows in sample_type. Please input another file.')
    
    # method 2                
    def prefilters(self,
                   separators: list = ['_', ' ', ':'],
                   suffixes: tuple = ('ligand', 'receptor'),
                   interaction_name: str = 'interaction_name'):
        """
        filter genes & interactions
        remove: 
            - genes if all row values are 0.
            - ineractions without two components in GEPs
            - genes not in interactons
        """
        # filter 1
        self.GEPs = sp.filter_genes_by_value(self.GEPs)
        
        # filter 2
        self.interactions = sp.filter_interactions_by_genes(self.interactions, self.GEPs, self.proteinComplex,
                                                            separators = separators,suffixes = suffixes, interaction_name = interaction_name)
        
        # filter 3
        self.GEPs = sp.filter_genes_by_interactions(self.GEPs, self.interactions)
        return self
    
    # method 3
    def infer_PPI(self, 
                  alpha: float = .05,
                  iterations: int = 2000, 
                  default_threads: int = 8) -> pd.DataFrame:
        """
        compute PPI enrichment probabiltiy of interactions from interacting cell types
        output: 
            1. real_stat data frame of mean_statistis of observations
            2. df_res data frame of probability;
            3. sig_stat data frame of significant mean_statistis of observations
               If p.value < alpha, the value will be the mean. Alternatively, the value is set to 0.
            4. df_PPI data frame with item of 1 or 0.
        """
        # interaction index (pd.dataframe)
        inter_idx = sp.get_interaction_index(self.GEPs, self.interactions)
    
        # sample pair index (pd.Series)
        sample_id = sp.get_sample_id(self.GEPs, self.cell_type) # list
        cellTypeComb = sp.get_cellType_combinations(self.cell_type) # list
        arr_sidx = sp.build_samplePairIndex_array(self.GEPs, self.cell_type, sample_id, cellTypeComb)

        # generate data frame to store results
        df_res = sp.build_result_df(self.GEPs, self.interactions, sample_id, cellTypeComb)
        
        # SECTION_1 compute mean_statistics for observed values
 
        # generate a list of tuples where each tuple is a combination of parameters
        paramlist = list(it.product(range(inter_idx.shape[0]), 
                                    range(len(sample_id)), 
                                    range(len(cellTypeComb))))
        
        # generate partial function which will process a tuple of parameters
        mean_stat_mp = partial(sp._mean_statistic_itparam, 
                                self.GEPs, 
                                inter_idx, 
                                arr_sidx, 
                                cellTypeComb)
        
        print('Calculating mean statistics for observed values. Please be patient.')    
        
        with Pool(processes = default_threads) as pool:  # context manager
            real_stat = pool.map(mean_stat_mp, paramlist)
                
        # create & reshape np.array
        real_stat = np.array(real_stat, dtype='float32').reshape(self.interactions.shape[0],
                                                            len(sample_id)*len(cellTypeComb))
        
        # SECTION_2 compute p-value based on null distribution from shuffled (columns of GEPs) data 
        # estimate null distribution once for each interaction
        
        # compute mean_statistics on shuffled data
        mean_stat_shuffled_mp = partial(sp._mean_stat_shuffled_meta_itparam,
                                            self.GEPs,
                                            inter_idx)
        
        paramlist2 = list(it.product(range(inter_idx.shape[0]), 
                                     range(iterations)))
        
        print('Infering PPI takes a while to complete. Please be patient.')
        
        with Pool(processes = default_threads) as pool:
            stat_shuf = pool.map(mean_stat_shuffled_mp, paramlist2)
     
        stat_shuf = np.array(stat_shuf, dtype='float32').reshape(real_stat.shape[0], 
                                                                  iterations)   
        
        # compute p-value
        compute_pvalue_mp = partial(sp._compute_pvalue_itparam,
                                    real_stat,
                                    stat_shuf,
                                    iterations)
        
        paramlist3 = list(it.product(range(real_stat.shape[0]), 
                                      range(real_stat.shape[1])))
        
        print('Computing p-values.')
        
        with Pool(processes = default_threads) as pool:
            pvalue = pool.map(compute_pvalue_mp, paramlist3)  # list
        
        pvalue = np.array(pvalue, dtype='float32').reshape(real_stat.shape[0], 
                                                            real_stat.shape[1])
            
        # update df_res
        df_res = sp.update_pvalue(df_res, pvalue)  # pd.dataframe of pvalue
        
        # convert p-value to 1 or 0
        df_PPI = sp.get_PPI(df_res, alpha)
        
        # # compute significant mean_statistics
        # sig_stat = sp.sig_mean_statistic(real_stat, df_res, alpha)
        
        # transform real_stat to pd.dataframe
        df_stat = pd.DataFrame(real_stat)
        df_stat.index = df_res.index
        df_stat.columns = df_res.columns
        
        self.df_stat = df_stat
        self.df_res = df_res
        # self.sig_stat = sig_stat
        self.df_PPI = df_PPI
        
        return self
        # return df_stat, df_res, sig_stat, df_PPI
    
