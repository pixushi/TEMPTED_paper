import os
import biom
import csv
import pandas as pd
import qiime2 as q2
import numpy as np
from qiime2 import Artifact
from qiime2 import Metadata
from qiime2.plugins.gemelli.actions import ctf


# obtain meta data and otu table
meta_table = pd.read_csv('data/metadata_cleaned.csv', sep=',', index_col=0)
otu_table = pd.read_csv('data/count_cleaned.csv', sep=',', index_col=0)
meta_table_ctf = meta_table.loc[:,['hostID', 'week','group']]
meta_table_ctf['week_int'] = [np.round(x) for x in meta_table_ctf['week']] 

# convert to qiime2 format
meta_table_ctf.index.name = '#SampleID'
meta_q2 = Metadata(meta_table_ctf)
otu_q2 = Artifact.import_data("FeatureTable[Frequency]", otu_table)

# run ctf using python directly
ctf_results  = ctf(otu_q2,
                   meta_q2,
                   'hostID',
                   'week_int',
                   min_feature_frequency=10,
                   n_components=3)
# save results
for id_, art_ in ctf_results.__dict__.items():
    if id_ != '_fields':
        art_.save(os.path.join('result',id_.replace('_', '-')+'_leukemia'))

