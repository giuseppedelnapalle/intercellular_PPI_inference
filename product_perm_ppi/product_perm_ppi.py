#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created in Jun 2020
@title: infering intercellular PPI based on LR_product score
@reference: Kumar_Cell_Reports_2018
@description:
    algorithm for infering intercellular PPI [optimised version].
    product of ligand and receptor expression in each cell type as statistic,
    permutation test used to test significicance of enrichment of PPI pairs.
    required functions imported from pp_ppi_helpers module
    version 1.0 (compatible with pp_ppi_helpers V1.0)
    $input 
        1. cell type expression profiles (pd_dataframe; derived form CIBERSORTx)
        2. protein pairs (pd.dataframe)
        3. cell types
        4. sample types of samples in gene expression profiles
    $output 
        1. data frame of product of interactions between two proteins 
           from specified cell type pairs (pd.dataframe)
        2. data frame of probabilty associated with product
        3. data frame of significant product of interactions
@author: nikola
"""

import numpy as np
import pandas as pd
import itertools as it
from functools import partial
from multiprocessing.pool import Pool
import pp_ppi_helpers as hp


# define class for the algorithm
# the ProductPermPPI object is created from the merged gene expression profile (GEP) matrix derived from CIBERSORTx.
# GEP object returned by merge_geps function in rpp_ppi_helpers module (pd.dataframe).
# genes in rows and samples in columns (#sample=n). first n columns for cell_type_1, second n columns for cell_type_2, etc.
# the class provides methods for intercellular PPI inference.
class ProductPermPPI:
    
    # method 1
    def __init__(self, 
                 geps: pd.DataFrame, 
                 cell_type: list, 
                 interactions: pd.DataFrame, 
                 sample_type: pd.DataFrame,
                 test_group: str = 'tumour',
                 exc_cell_type_ctrl: str = 'no_cell_type_excluded',  # e.g. 'Malignant'
                 ):
        """
        create a new ProductPermPPI instance
        remove defined cell_type, e.g. 'Malignant', in control samples, if argument exc_cell_type_ctrl provided
        """          
        # sample_ID (pd.series)
        sample_id = pd.Series(geps.columns[:int(geps.shape[1] / len(cell_type))])
        sample_id = sample_id.apply(lambda x: hp.strip_prefix(x, "".join(('^', cell_type[0], '_'))))
        
        # test if number of samples in mixture data & sample_type equal
        test = geps.shape[1]/len(cell_type) == sample_type.shape[0]
        
        if test:
            # test if argument exc_cell_type_ctrl given
            test2 = exc_cell_type_ctrl in cell_type
            
            if test2:
                rank_exc_ct = [i for i in range(len(cell_type)) if cell_type[i] == exc_cell_type_ctrl][0]
            
                # column index of to_be_removed_columns
                idx_rm = [i for i in range(sample_type.shape[0]) if sample_type.iloc[i].values != test_group]
                idx_rm = [rank_exc_ct*sample_type.shape[0] + i for i in idx_rm]
                # column names
                filt = geps.columns[idx_rm]
                
                # remove defined cell_type in control samples
                geps_filt = geps.loc[:, ~geps.columns.isin(filt)]
                
                # create new instance
                self.geps = geps_filt
                
            else:  # test2
                print('Argument exc_cell_type_ctrl not provided, no sample removed in gene expression profiles.')
                
                # create new instance
                self.geps = geps
                
            # append other attributes
            self.sample_id = sample_id
            self.cell_type = cell_type  
            self.interactions = interactions
            self.sample_type = sample_type
            self.test_group = test_group
            self.exc_cell_type_ctrl = exc_cell_type_ctrl
            
            print('Attributes of class ProductPermPPI: geps, sample_id, cell_type, \
                  interactions, sample_type, test_group, exc_cell_type_ctrl.')
            
        else:  # test==False
            print('Wrong number of rows in sample_type. Please input another file.')
            
    # method 2                
    def preprocess(self,
                   value_thresh: float = .01):
        """
        filter genes & interactions
        remove: 
            - genes if all row values are 0.
            - ineractions without two components in geps
            - genes not in interactons
        """
        # filter 1
        self.geps = hp.filter_genes_by_value(self.geps)
        
        # filter 2
        self.interactions = hp.filter_interactions_by_genes(self.interactions, self.geps)
        
        # filter 3
        self.geps = hp.filter_genes_by_interactions(self.geps, self.interactions)
        
        # transform values of 1 and NaN in geps
        arr = np.array(self.geps)
        arr = hp.convert_one(arr, value_thresh)
        arr[np.isnan(arr)] = value_thresh
        
        df = pd.DataFrame(arr)
        df.index = self.geps.index
        df.columns = self.geps.columns
        self.geps = df
        
        print('{}{}'.format('Shape of data frame geps: ', df.shape))
        print('{}{}'.format('Shape of data frame interactions: ', self.interactions.shape))
        
        return self
    
    # method 3
    def log2_transform(self):
        """
        log2 transform expression data
        """
        self.geps = np.log2(self.geps + 1)
        
        return self
    
    # method 4
    def compute_lr_product(self) -> pd.DataFrame:
        """
        compute lr_product_statistics measuring interaction potential from interacting cell types
        """
        # interaction index
        inter_idx = hp.get_interaction_index(self.geps, self.interactions)
        
        # sample pair index (pd.dataframe)
        sp_idx = hp.sample_pair_index_array(self.sample_id, self.cell_type, self.sample_type, 
                                            self.test_group, self.exc_cell_type_ctrl)
        
        # convert data frame to np.array for vectorized operation
        arr_geps = self.geps.to_numpy(dtype=np.float64)
        arr_inter_idx = inter_idx.to_numpy(dtype=np.int32)
        arr_sp_idx = sp_idx.to_numpy(dtype=np.int32)
        
        print('Computing lr_product statistics for observed values.')
        
        arr_stat = hp.compute_lr_product(arr_geps, arr_inter_idx, arr_sp_idx)
        
        # transform np.array to pd.dataframe
        df_stat = pd.DataFrame(arr_stat, index=self.interactions.index, columns=sp_idx.index)
        
        # append attributes
        self.df_stat = df_stat
        
        print('New attributes: df_stat.')
    
        return self
    
    # method_5
    def compute_pvalue(self, 
                       null_dist_ctrl: bool = True,
                       exact: bool = False,
                       iterations: int = 2000, 
                       default_threads: int = 8,
                       seed_num: int = 123):  # pd.dataframe added as attribute
        """
        compute emperical p-value based on null distribution from shuffled (columns of geps) data of control or all samples
        null distribution estimated from either control or all samples
        """
        # num control samples
        num_ctrl = self.sample_type.loc[self.sample_type.iloc[:, 0] != self.test_group].shape[0]
        print('{}{}{}'.format('Number of control samples ', num_ctrl, '.'))
        
        # interaction index
        inter_idx = hp.get_interaction_index(self.geps, self.interactions)
        
        # convert data frame to np.array for vectorized operation implementation
        arr_geps = self.geps.to_numpy(dtype=np.float64)
        arr_inter_idx = inter_idx.to_numpy(dtype=np.int32)
        
        arr_stat_shuf = hp.lr_product_shuffled(arr_geps, arr_inter_idx, self.cell_type, 
                                                self.sample_type, null_dist_ctrl, self.test_group, 
                                                self.exc_cell_type_ctrl, exact, iterations, seed_num)
        
        # compute p-value
        # define required parameter
        arr_stat = self.df_stat.to_numpy(dtype=np.float64)
        
        # check if arr_stat contains any NaN value
        test = np.isnan(arr_stat.sum())
        print('{}{}'.format('Does arr_stat contain any NaN value? ', test))
        
        if test:
            compute_pvalue_mp = partial(hp._compute_pvalue_itparam_2,
                                        arr_stat,
                                        arr_stat_shuf,
                                        iterations)
        else:
            compute_pvalue_mp = partial(hp._compute_pvalue_itparam,  # faster function
                                        arr_stat,
                                        arr_stat_shuf,
                                        iterations)
        
        paramlist = list(it.product(range(arr_stat.shape[0]), 
                                    range(arr_stat.shape[1])))  # list of tuple
        
        print('Computing p-values. Please be patient.')
        
        with Pool(processes = default_threads) as pool:
            pvalue = pool.map(compute_pvalue_mp, paramlist)  # list
    
        arr_p = np.array(pvalue, dtype=np.float64).reshape(arr_stat.shape[0], 
                                                            arr_stat.shape[1])
        
        # transform np.array to pd.dataframe
        df_pvalue = pd.DataFrame(arr_p, index=self.interactions.index, columns=self.df_stat.columns)
        self.df_pvalue = df_pvalue
        
        print('New attributes: df_pvalue.')
        
        return self
    
    # method_6
    def significant_lr_product(self, 
                                alpha: float = .05,
                                iterations: int = 2000, 
                                default_threads: int = 8,
                                value_thresh: float = .01,
                                seed_num: int = 123):  # pd.dataframe added as attribute
        """
        compute significant lr_product_statistic
        If p.value < alpha, the value will be the lr_product. Otherwise, the value is set to value_thresh.
        """
        # define required parameters
        arr_stat = self.df_stat.to_numpy(dtype=np.float64)
        arr_p = self.df_pvalue.to_numpy(dtype=np.float64)
        
        sig_lr_product_mp = partial(hp._sig_lr_product_itparam,
                                    arr_stat, arr_p,
                                    alpha, value_thresh)
        
        print('Computing significant lr_product statistics. Please be patient.')
        
        with Pool(processes = default_threads) as pool:
            sig_stat = pool.map(sig_lr_product_mp, range(arr_stat.shape[0]))  # list of np.array
        
        
        arr_sig_stat = np.array(sig_stat, dtype=np.float64)
        
        # transform np.array to pd.dataframe
        df_sig_stat = pd.DataFrame(arr_sig_stat, index=self.interactions.index, columns=self.df_stat.columns)
        self.df_sig_stat = df_sig_stat
        
        print('New attributes: df_sig_stat.')
        
        return self
