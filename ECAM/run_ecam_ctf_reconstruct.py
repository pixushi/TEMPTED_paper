import os
import csv
import pandas as pd
from pandas import concat
from pandas import DataFrame

import warnings
import numpy as np
np.seterr(divide = 'ignore')
from numpy.linalg import norm

import biom
from biom import Table, load_table

import skbio
from skbio import TreeNode, read, OrdinationResults, DistanceMatrix, TreeNode

from typing import Optional

from gemelli.factorization import TensorFactorization
from gemelli.ctf import ctf_table_processing
from gemelli.preprocessing import (build,
                                   fast_unifrac,
                                   bp_read_phylogeny,
                                   retrieve_t2t_taxonomy,
                                   create_taxonomy_metadata)
from gemelli._defaults import (DEFAULT_COMP, DEFAULT_MSC,
                               DEFAULT_MFC, DEFAULT_BL,
                               DEFAULT_MTD, DEFAULT_MFF,
                               DEFAULT_TENSALS_MAXITER,
                               DEFAULT_FMETA as DEFFM)
from gemelli.factorization import construct_tensor
# hide pandas Future/Deprecation Warning(s) for tutorial
warnings.filterwarnings("ignore", category=DeprecationWarning) 
warnings.simplefilter(action='ignore', category=FutureWarning)


### Functions for reconstruction

def reconstruction_helper(table: biom.Table,
                          sample_metadata: DataFrame,
                          individual_id_column: str,
                          state_column: str,
                          n_components: int = DEFAULT_COMP,
                          min_sample_count: int = DEFAULT_MSC,
                          min_feature_count: int = DEFAULT_MFC,
                          min_feature_frequency: float = DEFAULT_MFF,
                          max_iterations_als: int = DEFAULT_TENSALS_MAXITER,
                          max_iterations_rptm: int = DEFAULT_TENSALS_MAXITER,
                          n_initializations: int = DEFAULT_TENSALS_MAXITER,
                          feature_metadata: DataFrame = DEFFM):

    # check the table for validity and then filter
    state_column = [state_column]
    process_results = ctf_table_processing(table,
                                           sample_metadata,
                                           individual_id_column,
                                           state_column,
                                           min_sample_count,
                                           min_feature_count,
                                           min_feature_frequency,
                                           feature_metadata)
    (table, sample_metadata,
     all_sample_metadata, feature_metadata) = process_results

    # tables to dataframe
    table_df = DataFrame(table.matrix_data.toarray(),
                      table.ids('observation'),
                      table.ids('sample'))
    sample_metadata = sample_metadata.copy().reindex(table.ids())
    # tensor building
    tensor = build()
    tensor.construct(table_df, sample_metadata,
                     individual_id_column,
                     state_column,
                     branch_lengths=None)
    # factorize
    TF = TensorFactorization(
        n_components=n_components,
        max_als_iterations=max_iterations_als,
        max_rtpm_iterations=max_iterations_rptm,
        n_initializations=n_initializations, center=False).fit(tensor.rclr_transformed)
    # flatten for ease of testing
    M_rclr_obs = flatten_tensor(sample_metadata.copy(),
                                individual_id_column, state_column[0],
                                tensor.rclr_transformed.copy(),
                                tensor.subject_order.copy(),
                                tensor.condition_orders[0].copy(),
                                tensor.feature_order.copy())
    T_hat = construct_tensor(TF.eigvals, TF.loadings)
    M_hat_rclr = flatten_tensor(sample_metadata.copy(),
                                individual_id_column, state_column[0],
                                T_hat.copy(),
                                tensor.subject_order.copy(),
                                tensor.condition_orders[0].copy(),
                                tensor.feature_order.copy())
    # match tables
    tensor_table_df = pd.DataFrame(table.matrix_data.toarray(), 
                                   table.ids('observation'), table.ids())
    tensor_table_df = tensor_table_df.loc[M_rclr_obs.index, M_rclr_obs.columns]
    M_hat_rclr = M_hat_rclr.loc[M_rclr_obs.index, M_rclr_obs.columns]
    # ensure masked tables
    mask = np.zeros(tensor_table_df.shape)  # mask of zeros the same shape
    mask[abs(tensor_table_df.values) > 0] = 1  # set masked (missing) values to one
    M_rclr_obs = np.multiply(M_rclr_obs, mask)
    M_hat_rclr = np.multiply(M_hat_rclr, mask)
    
    return M_rclr_obs, M_hat_rclr
    #return tensor.rclr_transformed, T_hat
    

def flatten_tensor(sample_metadata,
                   individual_id_column,
                   state_column,
                   tensor_build,
                   subject_order,
                   condition_order,
                   feature_order):
    # construct matrix to make sanity checking easier
    table_matrix = []
    for i, subject_ in enumerate(subject_order):
        for j, state_ in enumerate(condition_order):
            id_ = sample_metadata[(sample_metadata[individual_id_column] == subject_) & (sample_metadata[state_column] == state_)]
            if id_.shape[0] == 0:
                continue
            id_ = id_.index[0]
            df_ = pd.DataFrame(tensor_build[i, :, j], feature_order, [id_])
            table_matrix.append(df_)
    table_matrix = pd.concat(table_matrix, axis=1).fillna(0)
    return table_matrix



### Simulation setup

rmax = 10
nsim = 100
subject_label = 'studyid'
state_label = 'visit'

count_tab = pd.read_csv("/Users/pixushi/Box/0-ftsvd/TEMPTED_paper/ECAM/data/otu_count_cleaned_q2.csv", index_col=0)

err = np.random.rand(nsim, rmax)*0
err = pd.DataFrame(err)
err.columns = [f'V{i}' for i in range(1, rmax+1)]
err.index = err.index + 1


for ss in range(nsim):
    print(ss)
    # Reading the metadata, setting the first column as the index
    meta_sub = pd.read_csv(f"/Users/pixushi/Box/0-ftsvd/TEMPTED_paper/ECAM/simdata/realsim{ss+1}_ntime9.csv", header=0, index_col=0)
    meta_sub['studyid'] = 'id_' + meta_sub['studyid'].astype('str')
    meta_sub['delivery_ind'] = meta_sub['delivery'] == 'Vaginal'
    meta_sub = meta_sub[['visit', 'studyid']]
    meta_sub['visit'] = meta_sub['visit'].astype('int')
    #meta_sub.index.name = '#SampleID'

    count_sub = count_tab.loc[meta_sub.index]
    # Filtering columns in count_sub where <= 95% of entries are zero
    count_sub = count_sub.loc[:, (count_sub == 0).mean() <= 0.95]

    table_sub = Table(count_sub.T.values, count_sub.columns, count_sub.index)

    for ii in range(rmax):
        T, T_hat = reconstruction_helper(table_sub.copy(),
                                                   meta_sub.reindex(table_sub.ids()),
                                                   subject_label, state_label,
                                                   n_components=ii+1)
        err.iloc[ss,ii] = 1 - (norm(T - T_hat)**2) / (norm(T_hat) ** 2)
    if ss % 10 == 9:
        err.to_csv("/Users/pixushi/Box/0-ftsvd/TEMPTED_paper/ECAM/result/ecam_reconstruction_ctf_nomean_ntime2.csv", index=True)
for ss in range(nsim):
    print(ss)
    # Reading the metadata, setting the first column as the index
    meta_sub = pd.read_csv(f"/Users/pixushi/Box/0-ftsvd/TEMPTED_paper/ECAM/simdata/realsim{ss+1}_ntime9.csv", header=0, index_col=0)
    meta_sub['studyid'] = 'id_' + meta_sub['studyid'].astype('str')
    meta_sub['delivery_ind'] = meta_sub['delivery'] == 'Vaginal'
    meta_sub = meta_sub[['visit', 'studyid']]
    meta_sub['visit'] = meta_sub['visit'].astype('int')

    count_sub = count_tab.loc[meta_sub.index]
    # Filtering columns in count_sub where <= 95% of entries are zero
    count_sub = count_sub.loc[:, (count_sub == 0).mean() <= 0.95]

    table_sub = Table(count_sub.T.values, count_sub.columns, count_sub.index)

    for ii in range(rmax):
        T, T_hat = reconstruction_helper(table_sub.copy(),
                                                   meta_sub.reindex(table_sub.ids()),
                                                   subject_label, state_label,
                                                   n_components=ii+1)
        err.iloc[ss,ii] = 1 - (norm(T - T_hat)**2) / (norm(T_hat) ** 2)
    if ss % 10 == 9:
        err.to_csv("/Users/pixushi/Box/0-ftsvd/TEMPTED_paper/ECAM/result/ecam_reconstruction_ctf_nomean_ntime9.csv", index=True)
print(err)

