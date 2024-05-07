import os
import biom
import csv
import pandas as pd
import qiime2 as q2
import numpy as np
import skbio
from qiime2 import Artifact
from qiime2 import Metadata
from qiime2.plugins.gemelli.actions import ctf

from gemelli.factorization import construct_tensor
from gemelli.preprocessing import build
from numpy.linalg import norm           



def reconstruct_tensor(table, metadata, 
                       subject_label, state_label,
                       state_subject_ordination,
                       subject_loadings, 
                       feature_loadings, 
                       state_loadings,
                       eigenvalues):
    """
    Generate the tensor and the reconstruction
    based on the loadings from CTF.
    """
    # table to dataframe
    table_df = table.transpose()
    # match the table to the end-result from CTF
    table_df = table_df.reindex(feature_loadings.index, axis=0)
    table_df = table_df.reindex(state_subject_ordination.index, axis=1)
    sample_metadata = metadata.copy().reindex(table_df.columns)
    # tensor building
    tensor = build()
    tensor.construct(table_df, sample_metadata,
                     subject_label,
                     [state_label],
                     branch_lengths=None)
    T = tensor.rclr_transformed.copy()
    # obtain mask array where values of tensor are zero
    mask = np.zeros(T.shape)  # mask of zeros the same shape
    mask[abs(T) > 0] = 1  # set masked (missing) values to one
    # build the reconstructed tensor
    U = [subject_loadings.reindex(tensor.subject_order).values, 
         feature_loadings.reindex(tensor.feature_order).values, 
         state_loadings.reindex(tensor.condition_orders[0]).values]
    T_hat = construct_tensor(eigenvalues.values, U)
    T_hat = np.multiply(T_hat, mask)
    return T, T_hat



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
    meta_sub['studyid'] = meta_sub['studyid'].astype(str)
    meta_sub['delivery_ind'] = meta_sub['delivery'] == 'Vaginal'
    meta_sub = meta_sub[['visit', 'studyid']]
    meta_sub['visit'] = meta_sub['visit'].astype('int')
    count_sub = count_tab.loc[meta_sub.index]
    # Filtering columns in count_sub where <= 95% of entries are zero
    count_sub = count_sub.loc[:, (count_sub == 0).mean() <= 0.95]

    # convert to qiime2 format
    meta_sub.index.name = '#SampleID'
    meta_q2 = Metadata(meta_sub)
    count_q2 = Artifact.import_data("FeatureTable[Frequency]",count_sub)

    for ii in range(rmax):
        ctf_results  = ctf(count_q2, meta_q2,
                                   subject_label,
                                   state_label,
                                   min_sample_count = 0,
                                   min_feature_count = 0, 
                                   min_feature_frequency = 0,
                                   max_iterations_rptm = 5 ,
                                   n_initializations = 5 ,
                                   max_iterations_als = 5,
                                   n_components = ii+1)

        # expand the results
        subject_biplot = ctf_results[0].view(skbio.OrdinationResults)
        state_biplot = ctf_results[1].view(skbio.OrdinationResults)
        distance_matrix = ctf_results[2].view(skbio.DistanceMatrix)
        state_subject_ordination = ctf_results[3].view(pd.DataFrame)

        subject_loadings=subject_biplot.samples.copy()
        feature_loadings=subject_biplot.features.copy()
        state_loadings=state_biplot.samples.copy()
        eigenvalues=subject_biplot.eigvals.copy()

        state_loadings.index = state_loadings.index.astype(float).astype(int)
        subject_loadings.index = subject_loadings.index.astype(float).astype(int).astype(str)

        # get reconstruction
        T, T_hat = reconstruct_tensor(count_sub, meta_sub, 
                                    subject_label, state_label,
                                    state_subject_ordination,
                                    subject_loadings, 
                                    feature_loadings, 
                                    state_loadings,
                                    eigenvalues)
        err.iloc[ss,ii] = 1 - np.sum((T - T_hat) ** 2) / np.sum(T ** 2)
    if ss % 10 == 0:
        err.to_csv("/Users/pixushi/Box/0-ftsvd/TEMPTED_paper/ECAM/result/ecam_reconstruction_ctf_nomean.csv", index=True)



