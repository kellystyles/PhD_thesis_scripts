#!/usr/bin/env/python3
"""
The aim of this script is to prepare 'n' databses containing 100 fungal genomes 
from two pools; one pool of IDT BGC negative genomes and the other of IDT BGC 
positive genomes.

Usage: python3 prepare_synthetic_dbs.py number_dbs_to_prepare

Notes:
 * negative genomes confirmed negative by BLASTn for core IDT biosynthetic genes
   (B,C,G,M)
 * take random selection of positive genomes and add negative genomes up to 100 
   for 'n' iterations; producing 'n' databases positive_dir negative_dir
 * this is completed in a blinded way; with the number of positives in each 
   database written to a file
"""
import sys
import os
import random
import glob
import shutil
import pandas as pd

negative = "negative"
positive = "positive"

def main():
    # randomly chooses positive genomes and copies them to the database folder, then fills the remaining spaces with negative genomes
    # change the range value to determine the number of databases prepared
    master_df = pd.DataFrame()

    # parse genome information
    df = pd.read_csv("genome_info.tsv", sep="\t", header=None)
    pos = df[df[4] == "positive"]
    neg = df[df[4] == "negative"]
    print(len(pos), len(neg))
    for i in range(int(sys.argv[1])):
        n = i + 1
        print(f"Preparing database {n} of {sys.argv[1]}")
        path = os.path.join("databases", str(n))
        os.makedirs(path, exist_ok=True)
        pos_rows = pos.sample(random.randint(0, 50))
        neg_rows = neg.sample(100 - len(pos_rows))
        df2 = pd.concat([pos_rows, neg_rows], axis=0)
        df2.columns = ["Accession", "Assembly", "Taxid", "Species", "IDTs"]
        df2 = df2.insert(loc=0, column='ID', value=n)
        pos_accs = pos_rows[0].to_list()
        neg_accs = neg_rows[0].to_list()

        try:
          for g in pos_accs:
              print(g)
              old_path = glob.glob(os.path.join("positive", g + '*.gbff'))[0]
              genome = old_path.split("/")[-1]
              print(genome)
              new_path = os.path.join("databases", str(n), genome)
              shutil.copy(old_path, new_path)
          for g in neg_accs:
              print(g)
              old_path = glob.glob(os.path.join("negative", g + '*.gbff'))[0]
              genome = old_path.split("/")[-1]
              print(genome)
              new_path = os.path.join("databases", str(n), genome)
              shutil.copy(old_path, new_path)
        except FileNotFoundError:
            print(f"{genome} not found :/")
            shutil.rmtree(os.path.join("databases", str(n)))
            break
        master_df = pd.concat([master_df, df2], axis=0)
    
    master_df.to_csv("databases.csv")

if __name__ == '__main__':
    main()