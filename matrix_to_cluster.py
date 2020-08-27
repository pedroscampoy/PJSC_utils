#!/usr/bin/env python

# Standard library imports
import os
import sys
import re
import logging

# Third party imports
import argparse
import datetime
import pandas as pd
import numpy as np 


logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
FUNCTION: Asign a group number according to the distance between samples supplied in csv lower left matrix format

INSTITUTION:CNM-ISCIII
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
d^v^b
VERSION=0.1
CREATED: 21 August 2020
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

def calculate_distance_stat(dataframe, list_sample=False):
    if list_sample != False:
        dataframe = dataframe.loc[list_sample,list_sample]
    stacked_df = dataframe.stack()
    #np.nanmean(dfdist.loc[cluster_test,cluster_test].values)
    #np.nanmin(dfdist.loc[cluster_test,cluster_test].values)
    #np.nanmax(dfdist.loc[cluster_test,cluster_test].values)
    mean_distance = stacked_df.mean(skipna = True)
    min_distance = stacked_df.min(skipna = True)
    max_distance = stacked_df.max(skipna = True)
    return ("This cluster has %s samples, with a mean distance on %.2f, range [%.0f - %.0f]" % (len(dataframe.columns), mean_distance, min_distance, max_distance))

def pairwise_to_cluster(pw,threshold = 20):
    groups = {}
    columns = pw.columns.tolist()
    sorted_df = pw[(pw[columns[0]] != pw[columns[1]]) & (pw[columns[2]] <= threshold)].sort_values(by=[columns[2]])
    
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
        
    for _, row in sorted_df.iterrows():
        group_number = len(groups)
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
            
    for _, row in pw[(pw[pw.columns[0]] != pw[pw.columns[1]]) & (pw[pw.columns[2]] > threshold)].iterrows():
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

def seqsphere_matrix(df):
    #df = matrix.iloc[:,:-1]
    index_name = df.index.name
    new_index = [index_name] + df.index.tolist()
    df.columns.name = None
    df.columns = new_index
    df1 = pd.DataFrame([[np.nan] * len(df.columns)], columns=new_index, index=[index_name])
    df2 = df1.append(df)
    df2 = df2.astype('float')
    return df2

def extraxt_decimal(file):
    with open(file, 'r') as f:
        content = f.read()
        content_list = content.split('\n')
        #half_line_no = int(len(content_list) / 2)
        second_line = content_list[2].split(',')[1:]
        if any(['.' in i for i in second_line]):
            return '.'
        else:
            return ','

def calculate_N(row):
    return len(row.samples)

def calculate_mean_distance(row, df):
    if row.N > 1:
        list_sample = row.samples
        dataframe = df.loc[list_sample,list_sample]
        stacked_df = dataframe.stack()
        mean_distance = stacked_df.mean(skipna = True)
        min_distance = stacked_df.min(skipna = True)
        max_distance = stacked_df.max(skipna = True)
        return round(mean_distance, 2), min_distance, max_distance
    else:
        return 'NaN'

    

def main():

    def get_arguments():

        parser = argparse.ArgumentParser(prog = 'common_mash_reference.py', description= 'Search for all mash files and find the representative reference')
        
        parser.add_argument('-i', '--input', dest="input_file", metavar="input_directory", type=str, required=True, help='REQUIRED.Input FASTA file')
        parser.add_argument('-o', '--output', type=str, required=False, default=False, help='Output directory to extract clusteres FASTA')
        parser.add_argument('-d', '--distance', type=float, required=False, default=15, help='Threshold distance to cluster sequences (default 15)')

        arguments = parser.parse_args()

        return arguments

    args = get_arguments()

    input_file = os.path.abspath(args.input_file)

    if args.output == False:
        output_dir = ('/').join(input_file.split('/')[0:-1])
    else:
        output_dir = os.path.abspath(args.output)
    
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
    logger.info('Reading Matrix')
    decimal = extraxt_decimal(input_file)
    dfdist = pd.read_csv(input_file, index_col=0, sep=",", decimal=decimal)
    if len(dfdist.index) != len(dfdist.columns) :
        logger.info('SeqSphere matrix detected')
        dfdist = seqsphere_matrix(dfdist)
    logger.info('Making pairwise')
    pairwise = dfdist.stack().reset_index(name='distance').rename(columns={'level_0': 'sample_1', 'level_1': 'sample_2'})
    clusters = pairwise_to_cluster(pairwise,threshold = args.distance)
    cluster_summary = clusters.groupby('group')['id'].apply(list).reset_index(name='samples')
    cluster_summary['N'] = cluster_summary.apply(calculate_N, axis=1)
    cluster_summary = cluster_summary.sort_values(by=['N'], ascending=False)
    logger.info('Reseting group number by length')
    sorted_index = cluster_summary.index.to_list()
    sorted_index.sort()
    sorted_index = [x + 1 for x in sorted_index]
    cluster_summary['group'] = sorted_index
    cluster_summary = cluster_summary.sort_values(by=['N'], ascending=False)

    cluster_summary[['mean', 'min', 'max']] = cluster_summary.apply(lambda x: calculate_mean_distance(x, dfdist), axis=1, result_type="expand")

    final_cluster = cluster_summary[["group", "samples"]].explode("samples").reset_index(drop=True)
    final_cluster = final_cluster.sort_values(by=['group'], ascending=True)

    overal_stats = calculate_distance_stat(dfdist)

    final_cluster_file = os.path.join(output_dir, "group_table_" + str(args.distance) + ".tsv")
    cluster_summary_file = os.path.join(output_dir, "group_summary_" + str(args.distance) + ".tsv")
    overall_summary_file = os.path.join(output_dir, "group_summary_stats_" + str(args.distance) + ".txt")

    cluster_summary.to_csv(cluster_summary_file, sep='\t', index=False)
    final_cluster.to_csv(final_cluster_file, sep='\t', index=False)
    with open (overall_summary_file, 'w+') as f:
        f.write(overal_stats)

    logger.info('DONE')

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise