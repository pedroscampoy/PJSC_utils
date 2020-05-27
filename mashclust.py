#!/usr/bin/env python

# Standard library imports
import os
import sys
import re
import logging
import subprocess

# Third party imports
import argparse
import datetime
import pandas as pd
import numpy as np 
from Bio import Entrez
from Bio import SeqIO

logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
FUNCTION: Reduces redundancy in multifasta files using kmer mash distance,
takes the longest sequence per cluster as representative

INSTITUTION:CNM-ISCIII
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
d^v^b
VERSION=0.1
CREATED: 27 May 2020
REVISION: 

TODO:

================================================================
END_OF_HEADER
================================================================
"""

END_FORMATTING = '\033[0m'
WHITE_BG = '\033[0;30;47m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
MAGENTA = '\033[35m'
BLUE =  '\033[34m'
CYAN = '\033[36m'
YELLOW = '\033[93m'
DIM = '\033[2m'

def check_create_dir(path):
    #exists = os.path.isfile(path)
    #exists = os.path.isdir(path)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

def mash_dist(input_file, output_dir, threads=10):
    # mash dist -i database.filtered_0.9_term.fasta database.filtered_0.9_term.fasta > database.filtered_0.9_term.mash.distances.tab
    #'reference_ID', 'query_ID', 'distance', 'p_value', 'shared_hashes'

    input_file = os.path.abspath(input_file)
    file_prefix = input_file.split('/')[-1]
    prefix = ('.').join(file_prefix.split('.')[0:-1])

    output_filename = prefix + ".mash.distances.tab"
    
    output_file = os.path.join(output_dir, output_filename)

    cmd = ["mash", "dist", "-i", "-p", str(threads), input_file, input_file]

    prog = cmd[0]
    param = cmd[1:]

    try:
        with open(output_file, "w+") as outfile:
            #calculate mash distance and save it in output file
            command = subprocess.run(cmd,
            stdout=outfile, stderr=subprocess.PIPE, universal_newlines=True)
        if command.returncode == 0:
            logger.info(GREEN + "Program %s successfully executed" % prog + END_FORMATTING)
        else:
            logger.info (RED + BOLD + "Command %s FAILED\n" % prog + END_FORMATTING
                + BOLD + "WITH PARAMETERS: " + END_FORMATTING + " ".join(param) + "\n"
                + BOLD + "EXIT-CODE: %d\n" % command.returncode +
                "ERROR:\n" + END_FORMATTING + command.stderr)
    except OSError as e:
        sys.exit(RED + BOLD + "failed to execute program '%s': %s" % (prog, str(e)) + END_FORMATTING)

def find_mash_dist_file(folder):
    folder = os.path.abspath(folder)
    files = os.listdir(folder)
    mash_file = [file for file in files if 'mash.distances.tab' in file][0]
    mash_file_path = os.path.join(folder, mash_file)

    return mash_file_path

def mash_dist_to_pairwise(distance_file, distance_type='hash_distance'):
    df = pd.read_csv(distance_file, sep='\t', names=['reference_ID', 'query_ID', 'distance', 'p_value', 'shared_hashes'])
    df[['hash_1', 'hash_2']] = df['shared_hashes'].str.split('/', expand=True)
    df.hash_1 = df.hash_1.astype(float)
    df.hash_2 = df.hash_2.astype(float)
    df['hash_distance'] = 1 - (df.hash_1 / df.hash_2)
    dfpair = df[['reference_ID', 'query_ID', distance_type]]
    
    return dfpair

def pairwise_to_cluster(pw,threshold = 0.5):
    groups = {}
    sorted_df = pw[(pw[pw.columns[0]] != pw[pw.columns[1]]) & (pw[pw.columns[2]] <= threshold)].sort_values(by=[pw.columns[2]])
    
    def rename_dict_clusters(cluster_dict):
        reordered_dict = {}
        for i, k in enumerate(list(cluster_dict)):
            reordered_dict[i] = cluster_dict[k]
        return reordered_dict
    
    def regroup_clusters(list_keys, groups_dict, both_samples_list):
        #sum previous clusters
        list_keys.sort()
        new_cluster = sum([groups_dict[key] for key in list_keys], [])
        #add new cluster
        cluster_asign = list(set(new_cluster + both_samples_list))
        #Remove duped cluster
        first_cluster = list_keys[0]
        groups_dict[first_cluster] = cluster_asign
        rest_cluster = list_keys[1:]
        for key in rest_cluster:
            del groups_dict[key]
        groups_dict = rename_dict_clusters(groups_dict)
        return groups_dict
        
    for index, row in sorted_df.iterrows():
        group_number = len(groups)
        cluster_name = 'cluster_' + str(group_number)

        sample_1 = str(row[0])
        sample_2 = str(row[1])
        both_samples_list = row[0:2].tolist()
                
        if group_number == 0:
            groups[group_number] = both_samples_list
        
        all_samples_dict = sum(groups.values(), [])
                
        if sample_1 in all_samples_dict or sample_2 in all_samples_dict:
            #extract cluster which have the new samples
            key_with_sample = {key for (key,value) in groups.items() if (sample_1 in value or sample_2 in value)}
            
            cluster_with_sample = list(key_with_sample)
            cluster_with_sample_name = cluster_with_sample[0]
            number_of_shared_clusters = len(key_with_sample)
            if number_of_shared_clusters > 1:
                groups = regroup_clusters(cluster_with_sample, groups, both_samples_list)
            else:
                groups[cluster_with_sample_name] = list(set(groups[cluster_with_sample_name] + both_samples_list))
        else:
            groups[group_number] = both_samples_list
            
    for index, row in pw[(pw[pw.columns[0]] != pw[pw.columns[1]]) & (pw[pw.columns[2]] > threshold)].iterrows():
        sample_1 = str(row[0])
        sample_2 = str(row[1])
        all_samples_dict = sum(groups.values(), [])
        if sample_1 not in all_samples_dict:
            group_number = len(groups)
            groups[group_number] = [sample_1]
        
        if sample_2 not in all_samples_dict:
            group_number = len(groups)
            groups[group_number] = [sample_2]
            
    cluster_df = pd.DataFrame(groups.values(),index=list(groups))
    
    cluster_df_return = cluster_df.stack().droplevel(1).reset_index().rename(columns={'index': 'group', 0: 'id'})
            
    return cluster_df_return
    

def main():

    def get_arguments():

        parser = argparse.ArgumentParser(prog = 'common_mash_reference.py', description= 'Search for all mash files and find the representative reference')
        
        parser.add_argument('-i', '--input', dest="input_file", metavar="input_directory", type=str, required=True, help='REQUIRED.Input FASTA file')
        parser.add_argument('-o', '--output', type=str, required=False, default=False, help='Output directory to extract clusteres FASTA')
        parser.add_argument('-d', '--distance', required=False, help='Threshold distance to cluster sequences[0-1] 0(identical) 1(unrelated) (default 0.5)')

        arguments = parser.parse_args()

        return arguments

    args = get_arguments()

    input_file = os.path.abspath(args.input_file)

    if args.output == False:
        output_dir = ('/').join(input_file.split('/')[0:-1])
    else:
        output_dir = os.path.abspath(args.output)
    
    print(output_dir)

    check_create_dir(output_dir)


    #LOGGING
    #Create log file with date and time
    right_now = str(datetime.date.today())
    right_now_full = "_".join(right_now.split(" "))

    log_filename = 'mashclust' + "_" + right_now_full + ".log"
    log_full_path = os.path.join(output_dir, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    #stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    #####################START PIPELINE################

    logger.info(args)

    #CALCULATE MASH DISTANCE
    logger.info('Obtaining mash distance')
    mash_dist(input_file, output_dir, threads=10)
    logger.info('Obtaining cluster from distance')
    mash_file = find_mash_dist_file('/home/pjsola/TMP/mashclust_test/')
    pairwise_distance = mash_dist_to_pairwise(mash_file)
    cluster_list = pairwise_to_cluster(pairwise_distance, threshold=args.distance)
    #EXTRACT FILES IN THE DESIRED FOLDER



if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise