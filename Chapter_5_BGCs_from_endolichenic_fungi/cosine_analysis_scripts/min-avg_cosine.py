#!/usr/bin/envpython
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from os import path
from scipy.spatial.distance import pdist, squareform
import csv
import os, argparse, sys

# define input arguments
csv_filename = sys.argv[1]
identifier=sys.argv[2]

# define output names
root_name = os.path.splitext(csv_filename)[0]
stat_filename= root_name + "-mincosine.csv"
stat_filename2= root_name + "-avgcosine.csv"

# parse data
bgc_features=pd.read_csv(csv_filename,index_col=0,header=None)

# convert to cosine vectors
X=pd.DataFrame(squareform(pdist(bgc_features,metric='cosine')), columns=bgc_features.index, index=bgc_features.index)

# parse metadata
index=pd.read_csv(identifier,header=None)[0].tolist()
col=[i for i in bgc_features.index.tolist() if i not in index]

# define the min cosine value for each BGC and write to file
X1=X.loc[index]
X2=X1.loc[:, col]
minValues0 = X2.min(axis=1)
minValues0.to_csv(stat_filename)
