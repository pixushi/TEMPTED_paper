import os
import biom
import csv
import pandas as pd
import qiime2 as q2
import numpy as np
from qiime2 import Artifact
from qiime2 import Metadata
from qiime2.plugins.gemelli.actions import ctf

# reading in the txt data table
npc = int(r.npc)
ss = str(int(r.ss))
pm = str(int(r.pm))
dp = str(int(r.depthkk))
# obtain meta data and otu table
meta_table = r.meta_sub
meta_table['study_day'] = meta_table['study_day'].astype('int')
meta_table = meta_table.drop(['study_group'], axis=1)
otu_table = r.count_sub

# convert to qiime2 format
meta_table.index.name = '#SampleID'
meta_q2 = Metadata(meta_table)
otu_q2 = Artifact.import_data("FeatureTable[Frequency]",otu_table)

# run ctf using python directly
ctf_results  = ctf(otu_q2, meta_q2,
                           'SubjectID',
                           'study_day',
			   min_sample_count = 0, 
			   min_feature_count = 0, 
			   min_feature_frequency = 0,
                           max_iterations_rptm = 5,
                           n_initializations = 5,
                           max_iterations_als = 5,
                           n_components = npc)


# save results
for id_, art_ in ctf_results.__dict__.items():
    if id_ != '_fields' and ('subject_biplot' in id_ or 'distance' in id_) :
        art_.save(os.path.join('simresult_ctf',id_.replace('_', '-')+'_sim'+ss+'_pmiss'+pm+'_depth'+dp+'K'))
