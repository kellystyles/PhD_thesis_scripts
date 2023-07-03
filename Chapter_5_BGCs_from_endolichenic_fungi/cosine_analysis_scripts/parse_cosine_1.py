#!/usr/bin/envpython
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import sqlite3
import csv
import os, argparse, sys

db_filename = sys.argv[1]
root_name = os.path.splitext(db_filename)[0]
stat_filename= root_name + "_hmm.csv"
# stat_filename2= root_name + "-dedupe_GCFC.csv"

con = sqlite3.connect(db_filename) # change to 'sqlite:///your_filename.db'
cur = con.cursor()
cur.execute("DROP TABLE IF EXISTS hmm;")
cur.execute("CREATE TABLE IF NOT EXISTS hmm (id,accession,name,db_id,model_length);")

with open('hmm.csv','r') as fin: # `with` statement available in 2.5+
    # csv.DictReader uses first line in file for column headings by default
    dr = csv.DictReader(fin) # comma is default delimiter
    to_db = [(int(i['id']),i['accession'], i['name'],str(i['db_id']),int(i['model_length'])) for i in dr]
cur.executemany("INSERT INTO hmm (id, accession, name, db_id, model_length) VALUES (?, ?, ?, ?, ?);", to_db)

conn = sqlite3.connect(db_filename, isolation_level=None,
                       detect_types=sqlite3.PARSE_COLNAMES)
db_df = pd.read_sql_query("SELECT * FROM bgc_features", conn)

con.commit()
con.close()

with sqlite3.connect(db_filename) as con:
    cur = con.cursor()
    bgc_ids, bgc_names = list(zip(*cur.execute("select id, orig_filename from bgc").fetchall()))
    hmm_ids, hmm_names = list(zip(*cur.execute("select id, name from hmm").fetchall()))
    bgc_features = pd.DataFrame(
        np.zeros((len(bgc_ids), len(hmm_ids)), dtype=np.uint8),
        index=bgc_ids,
        columns=hmm_ids
    )

    for bgc_id in bgc_ids:
        for hmm_id, value in cur.execute(("select hmm_id, value"
                                          " from bgc_features,hmm"
                                          " where bgc_features.bgc_id=?"
                                          " and hmm.id=bgc_features.hmm_id"
                                          " and hmm.db_id='1'"), (bgc_id,)).fetchall():
            bgc_features.at[bgc_id, hmm_id] = value

    bgc_features.index = bgc_names
    bgc_features.columns = hmm_ids

bgc_features0=bgc_features.loc[~(bgc_features==0).all(axis=1)]
bgc_features0.to_csv(stat_filename)
