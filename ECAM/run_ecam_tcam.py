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
# obtain meta data and clr transformed otu table
data_table = r.dat_sub
data_table = data_table.set_index(['studyid', 'timepoint'])
# format into data tensor
data_tensor, map1, map3 =  table2tensor(data_table)

# run TCAM
tca = TCAM()
tca_trans = tca.fit_transform(data_tensor)

# save results
tca_df = pd.DataFrame(tca_trans)
tca_df.rename(index = dict(map(reversed, map1.items())), inplace = True)  
tca_df = tca_df.iloc[:, :npc]
tca_df.to_csv(fnameout, sep='\t', encoding='utf-8')
