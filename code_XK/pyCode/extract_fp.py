import os
import pandas as pd
import pybel
from libs.paths import path_data,path_out_data
from libs.helper import SMILE2fp


data_test = pd.read_csv(os.path.join(path_data,"round_1_template.csv"))
fps = data_test.iloc[:,0].apply(SMILE2fp)

# convert a list of integers to a binary indcator vector
fps = pd.Series(fps)
fps_binary = pd.get_dummies(fps.apply(pd.Series).stack()).sum(level = 0)
fps_binary.columns =  [int(c) for c in fps_binary.columns]  

# convert each binary indcator vector to a list of integers to double check
cols = fps_binary.columns
fps_back = fps_binary.apply(lambda x: x > 0, raw = True).apply(lambda x: list(cols[x.values]), axis = 1)
num = 10 # check the first few to save time
equal_or_not = [fp_back != fp for fp_back, fp in zip(fps_back[:num], fps[:num])] 
assert(sum(equal_or_not) == 0)

# insert the SMILE column and save the finger prints
fps_binary.insert(loc=0, column="SMILE",value = data_test.iloc[:,0])
fps_binary.to_csv("../../data/fingerprint.csv",index = False)
