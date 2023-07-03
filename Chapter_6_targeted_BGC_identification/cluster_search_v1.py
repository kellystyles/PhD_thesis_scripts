#!/usr/bin/env python
"""
Created on Wed Feb 17 13:20:39 2021

@author: styleske
version 1

Usage: python3 --input 'gbk_dir' --directory 'pHMM_dir'

Please read prior to usage:
 -  check that there is no file called 'trusted_cutoffs.txt' as it will attempt 
    to parse the results from that file and the format changed between versions 
    2 and 3 of cluster_search.
"""

import os
from os import path
import subprocess
import sys, argparse
import ast
import fnmatch
from platform import python_version
import pandas as pd
import Bio
from Bio import SeqIO
from collections import Counter
from datetime import datetime

def get_args():
    """Parse command line arguments"""

    try:
        parser = argparse.ArgumentParser(
            description="Finds and extracts clusters of biosynthetic genes in an annotated genome using given HMM models")
        parser.add_argument('-i', '--input', action='store',
                            help='genbank file(s)')
        parser.add_argument('-d', '--directory', action='store',
                            help='directory with pHMMs and multi-FASTA file with proteins used to generate the pHMMs. Filenames must share the same basename i.e., idtB.hmm and idtB.fasta')

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            exit(1)
    except:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")

    return parser.parse_args()

def check_versions():
    
    print("Installed Software:\n")
    
    # Python modules
    #if python_version()
    print("Python: {}".format(python_version()))
    print("BioPython: {}".format(Bio.__version__))
    print("pandas: {}".format(pd.__version__))
    
    # Bash modules
    cmd="blastp --version"
    print(cmd)
    subprocess.run([cmd], stdout=subprocess.PIPE)
    
def get_cds(gbk):
    """
    Input a GBK file and returns a multi-FASTA file with all protein sequences.
    """
    from Bio import SeqIO
    
    if path.exists("cds") == False:
        os.mkdir("cds")

    input_handle  = open(gbk)
    out_filename = "cds/" + gbk.rsplit(".",1)[0].rsplit("/",1)[1] + ".fa"
    output_handle = open(out_filename, "w")

    for seq_record in SeqIO.parse(input_handle, "genbank") :
        for seq_feature in seq_record.features :
            if seq_feature.type=="CDS":
                try:
                    assert len(seq_feature.qualifiers['translation']) == 1
                    output_handle.write(">%s \n%s\n" % (
                           seq_feature.qualifiers['protein_id'][0],     # change protein_id to appropriate name value, depending on GBK source
                           seq_feature.qualifiers['translation'][0]))
                except KeyError:
                    continue
    
    output_handle.close()
    input_handle.close()
    
    return out_filename

def make_blastdb(infile):
    """
    Input a multi-FASTA protein sequence file and make blastdb
    """
    db_name = infile.rsplit(".",1)[0].split("/")[1]
    if path.exists("blast_db") == False:
        os.mkdir("blast_db")
    
    if os.path.isfile('blast_db/' + db_name + ".pin") == False:
        cmd = "makeblastdb -in " + infile + " -input_type fasta -dbtype prot -title " + db_name + " -out blast_db/" + db_name
        print(cmd + "\n")
        subprocess.run(cmd, shell=True)
    else:
        print("The BLASTDB for " + db_name + " has already been prepared; re-using...\n")

def hmmsearch(hmmprofile,seq_file,out):
    """
    Function that runs hhmsearch against sequence database using a given HMM profile
    input: profile HMM in HMM format;
           protein sequence file in multi-FASTA format
    output: table with scores, one line per homologous target found
    """
    
    cmd = "hmmsearch " + hmmprofile + " " + seq_file + " > " + out
    print("hmmsearch command:",cmd)
    subprocess.run(cmd, shell=True)
    
    return out

def trusted_cutoffs(directory):
    """
    Calculates trusted cutoff thresholds for each pHMM.
    Input is a directory containing the pHMMs and a multi-Fasta file of the protein sequences used to generate the pHMMs.
    Returns dictionary with trusted cutoffs for each pHMM.
    -------
    """
    
    TC = {}
    files = os.listdir(directory)
    output_file = "trusted_cutoffs.txt"

    if os.path.isfile(output_file) == False or os.stat(output_file).st_size == 0:
        for file in files:
            if file.endswith('.hmm'):
                filename = file.split(".")[0]
                hmmsearch("hmms/" + file, "hmms/" + filename + ".fasta", "hmms/" + filename+".txt")
                results = pd.DataFrame()
                clean_coll = []
                
                # read in input file and obtain hit lines
                with open("hmms/"+filename+".txt", 'r') as read_obj:
                    for line in read_obj:
                        cleaned = line.strip()
                        if cleaned.startswith(tuple('0123456789')) and ('!' not in cleaned) and ('?' not in cleaned) and ('*' not in cleaned) and ('PP' not in cleaned):
                            clean_coll.append(cleaned.split(maxsplit=9)[:9])
                read_obj.close()
                
                # add data to pandas dataframe
                results = pd.DataFrame(clean_coll)
                results.columns = ['seq_E-value', 'seq_score', 'seq_bias', 'dom_E-value', 'dom_score', 'dom_bias', 'exp', 'N', 'Accession']   
                
                seq_min = float(results.seq_score.min())
                dom_min = float(results.dom_score.min())
                
                # calculate TCs
                if len(results) >= 5:
                    c = 0.2
                else:
                    c = 0.5
                
                seq_TC = seq_min - (c * seq_min)
                dom_TC = dom_min - (c * dom_min)
                
                # write TCs to file
                TC[filename]={'Name':filename,'seq_TC':round(seq_TC,2),'dom_TC':round(dom_TC,2)}
                with open("trusted_cutoffs.txt", "w") as tc:
                    tc.write(str(TC))

    else:
        print("Trusted cutoffs already calculated; reusing...")
        TC = {}
        with open("trusted_cutoffs.txt", "r") as f:
            f_in = f.read()
            TC = ast.literal_eval(f_in)

    if len(TC) != len(fnmatch.filter(os.listdir(directory), '*.hmm')):
        sys.stderr.write("Not all trusted cutoffs could be calculated for your pHMMs. Please check your inputs\n")
        exit(1)

    return TC
    
def search_output(search_results, seq_TC, dom_TC):
    """
    Function that filters out hits below the trusted cutoff bitscores then outputs results to a dictionary
    input: HMMSEARCH output file, and trusted cutoff bitscores (TCs) for sequence (seq_TC) and domain (dom_TC)
    output: nested dictionary with hit information, bitscores, and e-values for sequence and domain hits
    """
    parsed = {}
    results = pd.DataFrame()
    clean_coll = []

    #read in input file and obtain hit lines
    with open(search_results, 'r') as read_obj:
        for line in read_obj:
            cleaned = line.strip()
            if cleaned.startswith(tuple('0123456789')) and ('!' not in cleaned) and ('?' not in cleaned) and ('*' not in cleaned) and ('PP' not in cleaned):
                clean_coll.append(cleaned.split(maxsplit=9)[:9])
            elif cleaned == "[No hits detected that satisfy reporting thresholds]":
                hmm = str(search_results).split("_")[0]
                raise Exception("No hits detected that satisfy reporting thresholds were detected for " + hmm)

    results = pd.DataFrame(clean_coll)
    results.columns = ['seq_E-value', 'seq_score', 'seq_bias', 'dom_E-value', 'dom_score', 'dom_bias', 'exp', 'N', 'Accession']

    # Filters out hits that have lower sequence and domain bitscores than the TCs
    for row in results.iterrows():
        df = results[results['seq_score'].astype(float) >= seq_TC]
        df2 = df[df['dom_score'].astype(float) >= dom_TC]
    
    prot_query = str(search_results).split("/")[1].split("_",1)[0]

    # makes dictionary from each row
    n = 1
    if len(df2) != 0:
        for row in df2.iterrows():
            parsed["hit_" + str(n)] = {
                'seq_E-value': df2['seq_E-value'][n-1],
                'seq_score': df2['seq_score'][n-1], 
                'seq_bias': df2['seq_bias'][n-1], 
                'dom_E-value': df2['dom_E-value'][n-1],  
                'dom_score': df2['dom_score'][n-1], 
                'dom_bias': df2['dom_bias'][n-1],
                'exp': df2['exp'][n-1], 
                'Accession': df2['Accession'][n-1], 
                'hmm_model': str(prot_query),
                'Hit_number': n
                }
            sys.stderr.write("\n")
            sys.stderr.write("HMM hit no. {}\n".format(n,parsed["hit_" + str(n)]['Hit_number']))
            sys.stderr.write("\tTop hit: {}\n".format(parsed["hit_" + str(n)]['Accession']))
            sys.stderr.write("\t%seq_E-value: {}\n".format(parsed["hit_" + str(n)]['seq_E-value']))
            sys.stderr.write("\tseq_score: {}\n".format(parsed["hit_" + str(n)]['seq_score']))
            sys.stderr.write("\tdom_E-value: {}\n".format(parsed["hit_" + str(n)]['dom_E-value']))
            sys.stderr.write("\tdom_score: {}\n".format(parsed["hit_" + str(n)]['dom_score']))
            n += 1
    else:
        print("No hits for " + prot_query)

    return parsed

def parse_gbk(infile, query):
    """
    Function that retrieves information from a gbk file
    input:  annotated genome genbank file in GENBANK format;
            query, a string with protein accession
    output: returns organism, contig, and sequence information
    """
    values = {}
    recs = [rec for rec in SeqIO.parse(infile, "genbank")]
    for rec in recs:
        for feat in rec.features:
            if feat.type == "CDS":
                try:       
                    if feat.qualifiers['locus_tag'][0] == query:
                        print("Query is locus_tag")
                        nquery = query
                    elif feat.qualifiers['protein_id'][0] == query:
                        nquery = feat.qualifiers['locus_tag'][0]
                        print("Query is protein_id")
                        print("New query =",nquery)
                except KeyError:
                    continue

    for rec in recs:
        feats_gene = [feat for feat in rec.features if feat.type == "gene"]
        feats_cds = [feat for feat in rec.features if feat.type == "CDS"]
        try:    
            for feat in feats_cds:
                locus_tag = feat.qualifiers['locus_tag'][0]            # get source (gene accession)
                protein_id = feat.qualifiers['protein_id'][0]
                if locus_tag == nquery:
                    organism = rec.annotations['organism']        # gets organism
                    contig = rec.id                                     # contig is the rec.id
                    seq = feat.qualifiers['translation']
                    loc_cds = feat.location
                    print("Hit found:",locus_tag,"on contig",contig)
                    values['Organism'] = organism
                    values['Contig'] = contig
                    values['Sequence'] = seq
                    values['location_cds'] = loc_cds
                    values['Locus_tag'] = locus_tag
                    values['Protein_id'] = protein_id             
        except KeyError:
            continue
        try:
            for feat in feats_gene:            
                    locus_tag = feat.qualifiers['locus_tag'][0]
                    if locus_tag == nquery:
                        loc_mrna = feat.location
                        values['location_gene'] = loc_mrna
        except KeyError:
            continue
    
    return values     

def update_parsed(infile, parsed):
    """
    Function that updates 'parsed' with contig/scaffold and sequence information
    input: annotated genbank file in GENBANK format;
           nested dictionary - the output of the search_output() function
    output: updates 'parsed' with contig and sequence information
    """
    n = 1
    for i in range(len(parsed)): # loops through nested dictionary 'parsed'
        query = parsed['hit_' + str(n)]['Accession']   # retrieves the accession for each hit in parsed    
        match = parse_gbk(infile, query)                    # finds the matching contig and sequence info for this query
        parsed['hit_' + str(n)].update(match)
        n += 1

    return parsed

def blastp(parsed, db):
    """
    Function that loops through nested dictionary and BLASTp each sequence against the nr database
    input: parsed = output nested dictionary of the search_output() function (dictionary);
           image = Singularity image with NCBI-BLAST tools installed in SIF format (string);
           db = 'nr' NCBI database location (string);
    output: updates parsed with a further nested dictionary with top blast hit (dictionary)
    """
    n = 1
    for i in range(len(parsed)):                            # loops through nested dictionary 'parsed'
        hit_num = str(parsed['hit_' + str(n)]['Hit_number'])
        prot_query = parsed['hit_' + str(n)]['hmm_model']
        acc_query = parsed['hit_' + str(n)]['Accession']   # retrieves the accession for each hit in parsed
        seq = parsed['hit_' + str(n)]['Sequence']
        if path.exists("blast_output") == False:
            os.mkdir("blast_output")
        filename = "blast_output/" + acc_query + ".fasta"
        output = "blast_output/" + prot_query + "_" + hit_num + "_" + acc_query + "_blastp.out"

        if os.path.isfile(filename) == False:
            sys.stderr.write("infile: {}\n".format(filename))
            sys.stderr.write("outfile: {}\n".format(output))
            with open(filename, 'w') as f:
                f.write(">" + acc_query + '\n')
                f.write(str(seq[0]))
        else:
            print("Using existing fasta file for record " + acc_query)    
        
        if os.path.isfile(output) == False:
            cmd1 = "blastp" + " -query " + filename + \
            " -db " + db + " -outfmt" +  " \"6" + " std\"" + \
            " > " + output
            sys.stderr.write("Running BLASTp to find the top hit for {}...\n".format(acc_query))
            try:
                print("BLASTp command:",cmd1)
                subprocess.run(cmd1, shell=True) # doesn't work with Jupyter Notebook
            except subprocess.CalledProcessError as e:
                print(e.output)
        else:
            print("Reusing BLASTp results for " + acc_query)

        print("opening blast results")
        with open(output, "r") as o:
            first_line = o.readline().rstrip().split('\t')
            blastp_hit = first_line[1]
            blastp_perid = first_line[2]
            blastp_eval = first_line[10]
            blastp_bits = first_line[11]
            blastp = {"\tTop_hit":blastp_hit,"%_identity":blastp_perid,"E-value":blastp_eval,"Bitscore":blastp_bits}
            sys.stderr.write("\n")
            sys.stderr.write("BLASTp hit no. {}\n".format(n))            
            sys.stderr.write("\tTop hit: {}\n".format(blastp_hit))
            sys.stderr.write("\t% identity: {}\n".format(blastp_perid))
            sys.stderr.write("\tE-value: {}\n".format(blastp_eval))
            sys.stderr.write("\tBitscore: {}\n".format(blastp_bits))
            parsed['hit_' + str(n)]['blastp'] = blastp
        
        n += 1

def count_contigs(listed_dicts):
    """
    Function that retrieves all contig hits from multiple dictionaries and counts occurences of each contig
    """

    new_vals = []
    for i in listed_dicts:
        for x in range(len(i)):
            x += 1
            new_vals.append(i['hit_' + str(x)]['Contig'])        

    counted = Counter(new_vals).most_common() # counts each element in the list

    return counted

def combine_hits(contig, list_hits, i):    
    """
    Combines hits that match the same scaffold into a single dictionary
    """
    x = {}

    for index in range(len(list_hits)):
        for key in list_hits[index]:
            if list_hits[index][key]['Contig'] == contig:
                hit_name = list_hits[index][key]['Name'] 
                x[hit_name] = list_hits[index][key]
    
    return x

def split_hits(listed_dicts):
    """
    Splits each hit into individual dictionaries and appends them to a list; 
    counts the number of contigs;
    applies the combine_hits() funtion to hits for each contig;
    returns a dictionary with all hits sorted by each contig
    """

    list_hits = [] 
    counted = count_contigs(listed_dicts)    
    
    for dictionary in listed_dicts:
        n = 1
        if len(dictionary) <= len(counted) and len(dictionary) != 0:
            new_name = dictionary['hit_' + str(n)]['hmm_model']
            n += 1
            for i in dictionary:
                dict_new = {}
                iter_name = new_name + "_" + str(i)
                dict_new[iter_name] = dictionary[i]
                dict_new[iter_name]['Name'] = iter_name
                list_hits.append(dict_new)
    nn = 1
    final = {}
    for obj in counted:
        name = obj[0]
        final[name] = combine_hits(obj[0], list_hits, nn)
        print("Cluster_" + str(nn) + " on " + name + " has " + str(len(final[name])) + " hits!")
        nn += 1
    
    return final

def list_loci(ins):
    """
    Returns all locus_tags in dictionary as a list (as a nested list for each cluster), 
    if there is more than one gene in the cluster
    """

    bb = []
    for i in ins.keys():
        lst = []
        print("i = ",i)
        for j in ins[i].keys():
            print("j =",j)
            print(ins[i][j]["Locus_tag"])
            lst.append(ins[i][j]["Locus_tag"])        

        bb.append(lst) 
        # uncomment the '#' to and add '#' to the if statement lines below to include single hits
        #if len(ins[i]) > 1:
         #   bb.append(lst)
        #else:
         #   print("\tCluster_{} on {} was removed as it did not have enough hits for a cluster".format(n, i))
          #  i += 1
    return bb 

def get_range(infile, gs, n):
    """
    Takes locus_tags for two genes, gets range (the distance between these two genes), 
    and extracts 15 kb on either side of this range, outputting as a Genbank file
    Input: infile = target Genbank file;
           g1,g2 = genes to use for targetting range;
           n = iteratable integer for labelling purposes
    This code is adapted from "https://www.biostars.org/p/340270/" by user Joe, with some modifications
    """
    print(gs)

    print('Fetching cluster_{} - {}\n'.format(n,gs))
    
    recs = SeqIO.parse(infile, 'genbank')   
    
    for rec in recs:
        loci = [feat for feat in rec.features if feat.type == "CDS"]
        organism = rec.annotations['organism']
        try:
            start = min([int(l.location.start) for l in loci if l.qualifiers['locus_tag'][0] in gs])
            end = max([int(l.location.end) for l in loci if l.qualifiers['locus_tag'][0] in gs])
            print("Located records in {}!\n".format(rec.id))
        except ValueError:
            continue
        
        try:
            (start and end)
            print("Start position: {} bp \nEnd position: {} bp\n".format(start,end))
            left_edge = start - 10000
            right_edge = end + 10000
            subrecord = []
            subrecord = rec[left_edge:right_edge]
            organism_ = organism.replace(" ","_")
            cluster = "cluster_" + str(n) + "_" + infile.rsplit("/")[1].rsplit(".")[0]
            cluster = cluster.rsplit(".",1)[0]
            SeqIO.write(subrecord, "clusters/" + cluster + ".gb", "genbank")           
        except NameError:
            print('Didn\'t get any indices even though the genes seemed to match. Couldn\'t slice.\n')
            pass

def dict_print(d):
    """
    Prints out all elements of a nested dictionary - uses a recursive method
    """
    for k, v in d.items():
        if isinstance(v, dict):
            print(k)
            dict_print(v)
        else:
            print("{0} : {1}".format(k, v))

def update_cluster(cluster, dict):
    """
    Adds predicted protein name to cluster Genabnk files as 'product'
    """
    parsed = SeqIO.parse(cluster, 'genbank')
    cluster_name = cluster.split("/")[1].rsplit(".",1)[0]
    hit_num = cluster_name.split("_",2)[1]

    for record in parsed:    
        for feat in record.features:
            if feat.type == "CDS" or feat.type == "CDS":
                for x in dict:
                    for y in dict[x]:
                        if dict[x][y]['Locus_tag'] == feat.qualifiers['locus_tag'][0]:
                            feat.qualifiers['product'] = dict[x][y]['hmm_model']

    SeqIO.write(record, cluster, "genbank")

def main():
    sys.stderr.write("~~~ Cluster Searcher ~~~\n")
    sys.stderr.write("------------------------\n")
    sys.stderr.write(str(datetime.now()).split(".")[0])
    sys.stderr.write("\n")
    
    args = get_args()
    input = args.input
    hmmdir = args.directory
    # bdb = args.blastdb
    # af = args.assembly

    sys.stderr.write("--- Parameters ---\n")
    sys.stderr.write("GBK file: {}\n".format(input))
    sys.stderr.write("HMM directory: {}\n".format(hmmdir))
    sys.stderr.write("\n")

    # Calculates or retrieves trusted cut off values for each pHMM
    sys.stderr.write("--- Calculating trusted cutoff bitscores for each pHMM ---\n")
    TC = trusted_cutoffs(hmmdir)
    #for i in TC:
    #    sys.stderr.write("{} TCs:\n".format(TC[i]['Name']))
    #    sys.stderr.write("Seq_TC: {}\n".format(TC[i]['seq_TC']))
    #    sys.stderr.write("Seq_TC: {}\n".format(TC[i]['dom_TC']))   
    sys.stderr.write("\n")

    # Loops over every pHMM file and every CDS file
    phmms = os.listdir(hmmdir)
    gbk = input
    list_of_dicts = []
    
    if path.exists("hmm_output") == False:
        os.mkdir("hmm_output")
    for phmm in phmms:
        if phmm.endswith('.hmm'):
            # A try loop so if one pHMM doesn't have any hits it will continue to search with other pHMMs
            try:
                print(phmm,gbk)   
                hmm_name = phmm.rsplit(".",1)[0]
                cd_name = gbk.rsplit(".",1)[0].rsplit("/",1)[1]
                print(hmm_name,cd_name)
                
                # Makes multi-fasta file of protein seqs from Genbank file and makes BLAST db from it; exits if there are no CDS seqs
                af = get_cds(gbk)
                sys.stderr.write("--- Searching for hits in {} with {} ---\n".format(af,phmm))
                if os.stat(af).st_size == 0:
                    print(af + " has no CDS records; exiting...")
                    pass
                else:
                    print("multi-FASTA CDS file:",af)
                    make_blastdb(af)

                out = phmm.rsplit(".",1)[0] + "_" + cd_name + ".txt"
                
                # Checks if the pHMM has already been searched against the infile, otherwise searches with pHMM against infile
                if os.path.isfile("hmm_output/" + out) == False:
                    search_out = hmmsearch("hmms/" + phmm, af, "hmm_output/" + out)
                else:
                    search_out = "hmm_output/" + out
                    print("Reusing hmmsearch results from " + search_out)
                
                # Filters positive pHMM hits with trusted cutoff values via a pandas dataframe
                file_dict = hmm_name + "_hits"
                file_dict = search_output(search_out,TC[hmm_name]['seq_TC'],TC[hmm_name]['dom_TC']) # filtering HMMsearch hits by trusted cutoffs
                
                # Retrieves information for each hit from input genbank file
                update_parsed(gbk,file_dict)                   # adding information from GBK file to dictionary
                
                # BLASTp the hit against the BLAST database ('nr' by default)
                blast_out = "blast_db/" + cd_name
                blastp(file_dict,blast_out)                         # BLASTp the blast database with the top hit accession from HMMsearch
                list_of_dicts.append(file_dict)                     # adding all info to main dictionary
                sys.stderr.write("\n")

            except Exception:

                pass
    sys.stderr.write("--- Combining results ---\n")
    
    # Prints out flattened dictionary to see all the info
    for i in list_of_dicts:
        dict_print(i)
    
    # Splits all the hits into separate dictionaries based on which contig the hit is on
    hits = split_hits(list_of_dicts)            # splitting each hit into individual nested dictionaries for each cluster within a main dictionary
    sys.stderr.write("\n")

    # Makes a list with just the accessions of each hit from 'hits' and removes clusters with one gene (could turn off for single-copy orphans)
    sys.stderr.write("Removing single gene clusters...\n")
    accessions = list_loci(hits)                # making a list with just the accessions
    sys.stderr.write("\n")
    sys.stderr.write("--- Extracting clusters ---\n")
    if path.exists("clusters") == False:
        os.mkdir("clusters")    

    # Extracts each cluster from genbank infile  
    cluster_num = 1
    for i in range(len(accessions)):
        #print(i)
        #for infile in gbks:
           # get_range(input + "/" + infile,accessions[i],nn)          # extracting clusters
          #  nn += 1
        get_range(gbk,accessions[i],cluster_num)
        cluster_num += 1

    # Adds product names to each hit in each cluster genbank file
    clusters = os.listdir("clusters")
    for cluster in clusters:
        update_cluster("clusters/" + cluster, hits)

    # Reports the number of clusters extracted
    #print("LEN_ACCESSIONS = ",len(accessions))
    if len(accessions) == 1:
        sys.stderr.write("{} cluster was found and is available in the \'clusters\' directory\n".format(len(accessions)))
    elif len(accessions) == 0:
        sys.stderr.write("No clusters found")
    else:
        sys.stderr.write("{} clusters were found and are available in the \'clusters\' directory\n".format(len(accessions)))
    sys.stderr.write("\n")

    sys.stderr.write("\n")
    sys.stderr.write("Script finished!\n")
    
if __name__ == '__main__':
    main()

# Written by Kelly Styles
# Version 1