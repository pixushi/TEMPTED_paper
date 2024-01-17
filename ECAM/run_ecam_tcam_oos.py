import os
import csv
import pandas as pd
import numpy as np
from mprod import table2tensor
from mprod.dimensionality_reduction import TCAM

# get the directory where the current script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# set the current working directory to the script directory
os.chdir(script_dir)

# reading in the txt data table
npc = int(r.npc)
fnameout = str(r.fnameout)
# obtain meta data and log 2 fold over baseline transformed otu table
data_train = r.dat_train
data_train = data_train.set_index(['studyid', 'timepoint'])
data_test = r.dat_test
data_test = data_test.set_index(['studyid', 'timepoint'])

# format into data tensor
data_tensor_train, map1_train, map3_train =  table2tensor(data_train)
data_tensor_test, map1_test, map3_test =  table2tensor(data_test)

# run TCAM on training data
tca = TCAM()
tca_trans_train = tca.fit_transform(data_tensor_train)
tca_trans_test = tca.transform(data_tensor_test)

# save results
tca_df_train = pd.DataFrame(tca_trans_train)
tca_df_train.rename(index = dict(map(reversed, map1_train.items())), inplace = True)  
tca_df_train = tca_df_train.iloc[:, :npc]
tca_df_train.to_csv(fnameout+'_train.csv', sep='\t', encoding='utf-8')

tca_df_test = pd.DataFrame(tca_trans_test)
tca_df_test.rename(index = dict(map(reversed, map1_test.items())), inplace = True)  
tca_df_test = tca_df_test.iloc[:, :npc]
tca_df_test.to_csv(fnameout+'_test.csv', sep='\t', encoding='utf-8')
