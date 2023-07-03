#!/usr/bin/envpython
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from os import path
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import AgglomerativeClustering

bgc_features0=pd.read_csv('fungal_hmm.csv',index_col=0,header=None)

X=pd.DataFrame(squareform(pdist(bgc_features0,metric='cosine')), columns=bgc_features0.index.unique(), index=bgc_features0.index.unique())
X.to_csv('fungal_matrix.csv')

gcf_cluster = AgglomerativeClustering(n_clusters=None, linkage="average", affinity="precomputed", distance_threshold=0.2)
gcc_cluster = AgglomerativeClustering(n_clusters=None, linkage="average", affinity="precomputed", distance_threshold=0.8)
gcf = ['GCF_'+str(i) for i in gcf_cluster.fit_predict(X)]
gcc = ['GCC_'+str(i) for i in gcc_cluster.fit_predict(X)]

pd.DataFrame ([list(a) for a in zip (X.index.tolist(),gcf,gcc)], columns = ['BGC', 'GCF','GCC']).to_csv('fungal_GCFC.csv')