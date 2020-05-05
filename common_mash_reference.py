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
from ftplib import FTP

#pd.set_option('display.max_colwidth', None)

logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
FUNCTION: search for all mash files (ends with .screen.tab), find and count the representative reference.
    Also download the assembly or all files from the ftp server in a specific folder

INSTITUTION:CNM-ISCIII
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
d^v^b
VERSION=0.1
CREATED: 24 February 2020
REVISION: 

TODO:
    independent scripts: common reference | download reference
================================================================
END_OF_HEADER
================================================================
"""

def check_create_dir(path):
    #exists = os.path.isfile(path)
    #exists = os.path.isdir(path)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)

def return_best_complete_match(screen_file):
    df = pd.read_csv(screen_file, sep='\t', names=['identity', 'shared-hashes', 'median-multiplicity', 'p-value', 'query-ID', 'query-comment'])
    df_complete = df[df['query-comment'].str.contains('complete genome')]
    df_complete = df_complete[~df_complete['query-comment'].str.contains(' phage')]
    df_complete = df_complete[~df_complete['query-comment'].str.contains(' Phage')]
    df_complete = df_complete[~df_complete['query-comment'].str.contains('shotgun')]
    df_complete = df_complete.sort_values(by=['identity'], ascending=False)
    df_complete.reset_index(inplace=True)
    logger.debug(df_complete.head())
    return df_complete.iloc[0][['query-ID','query-comment']].tolist()

def return_best_match(screen_file):
    df = pd.read_csv(screen_file, sep='\t', names=['identity', 'shared-hashes', 'median-multiplicity', 'p-value', 'query-ID', 'query-comment'])
    #df_complete = df[df['query-comment'].str.contains('complete genome')]
    df_complete = df.sort_values(by=['identity'], ascending=False)
    df_complete.reset_index(inplace=True)
    logger.debug(df_complete.head())
    return df_complete.iloc[0][['query-ID','query-comment']].tolist()

def find_common_reference(folder):
    '''
    Dependencies:   -input folder
                    -return_best_complete_match(filename)
                    -return_best_match(filename)
    '''
    #create empty df
    counter_record_complete = pd.DataFrame(columns=['query-ID','query-comment'])
    counter_record_all = pd.DataFrame(columns=['query-ID','query-comment'])
    #Create output tab files
    output_complete_tab = os.path.join(folder, 'counter_complete_mash.tab')
    output_all_tab = os.path.join(folder, 'counter_all_mash.tab')
    #Find common reference from mash result
    for root, _, files in os.walk(folder):
        for name in files:
            if name.endswith("screen.tab"):
                filename = os.path.join(root, name)
                try:
                    best_match_complete = return_best_complete_match(filename)
                    counter_record_complete.loc[len(counter_record_complete)] = best_match_complete
                except:
                    logger.debug("Format in " + name + "incorrect")
                try:
                    best_match = return_best_match(filename)
                    counter_record_all.loc[len(counter_record_all)] = best_match
                except:
                    logger.debug("Format in " + name + "incorrect")
    #pd.rename_axis and reset_index turn count_values() into a dataframe
    #https://stackoverflow.com/questions/47136436/python-pandas-convert-value-counts-output-to-dataframe
    #df = value_counts.rename_axis('unique_values').to_frame('counts')
    #counter_comment = counter_record['query-comment'].value_counts().rename_axis('unique_values').reset_index(name='counts')
    #counter_id = counter_record['query-ID'].value_counts().rename_axis('unique_values').reset_index(name='counts')
    counter_df_complete = counter_record_complete.groupby(['query-comment', 'query-ID']).size().reset_index(name='counts')\
    .sort_values(by=['counts'], ascending=False).reset_index(drop = True)
    
    counter_df = counter_record_all.groupby(['query-comment', 'query-ID']).size().reset_index(name='counts')\
    .sort_values(by=['counts'], ascending=False).reset_index(drop = True)
    
    gfc_complete = ('_').join(counter_df_complete.iloc[0]['query-ID'].split('_')[0:2])
    gfc_description = re.sub(r'[ ]?\[.{1,9}\][ ]?','',counter_df_complete.iloc[0]['query-comment'])
        
    counter_df_complete.to_csv(output_complete_tab, sep='\t', index=False)
    counter_df.to_csv(output_all_tab, sep='\t', index=False)

    return gfc_complete, gfc_description

def gcf_to_ftp_path(gcf_value):
    try:
        dfrefseq = pd.read_csv('ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt', skiprows=0, sep='\t', header=1)
    except:
        logger.info('There was a problem obtaining assembly_summary.txt\n \
        Check: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt')
        sys.exit(1)
    ftp_path = dfrefseq['ftp_path'][dfrefseq['# assembly_accession'] == gcf_value]
    ftp_path_no_domain = ftp_path.values[0].split('.gov')[-1]
    
    return ftp_path_no_domain

def download_gcf (ftp_address, output_dir, all_data=False ):
    
    output_dir = os.path.abspath(output_dir)
    
    if not os.path.exists(output_dir):
        logger.info("path " + output_dir + " doesn't exist and will be created")
        try:
            os.mkdir(output_dir)
        except:
            logger.info("Folder " + output_dir + "can't be created")
            sys.exit(1)
    
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    
    # Get All Files
    ftp.cwd(ftp_address)
    files = ftp.nlst()
    ftp_folder = ftp_address.split('/')[-1]
    assembly_file = ftp_folder + '_genomic.fna.gz'
    
    # logger.info out the files
    if all_data == False:
        for file in files:
            if file == assembly_file:
                local_path = os.path.join(output_dir, file)
                with open(local_path, 'wb') as f:
                    logger.info("Downloading.." + file)
                    #ftp.retrbinary("RETR " + file ,open(output_dir + file, 'wb').write)
                    ftp.retrbinary('RETR ' + file, f.write)
    else:
        for file in files:
            local_path = os.path.join(output_dir, file)
            with open(local_path, 'wb') as f:
                logger.info("Downloading.." + file)
                #ftp.retrbinary("RETR " + file ,open(output_dir + file, 'wb').write)
                ftp.retrbinary('RETR ' + file, f.write)

    ftp.close()


def main():

    def get_arguments():

        parser = argparse.ArgumentParser(prog = 'common_mash_reference.py', description= 'Search for all mash files and find the representative reference')
        
        parser.add_argument('-i', '--input', dest="input_dir", metavar="input_directory", type=str, required=True, help='REQUIRED.Input directory containing all screen.tab files')
        parser.add_argument('-o', '--output', type=str, required=True, help='REQUIRED. Output directory to extract all results')
        parser.add_argument('-a', '--all', required=False, action='store_true', help='Download all files instead of just the fasta')

        arguments = parser.parse_args()

        return arguments

    args = get_arguments()
    

    output_dir = os.path.abspath(args.output)
    input_dir = os.path.abspath(args.input_dir)
    
    check_create_dir(output_dir)

    #LOGGING
    #Create log file with date and time
    right_now = str(datetime.date.today())
    right_now_full = "_".join(right_now.split(" "))

    log_filename = 'mash_reference' + "_" + right_now_full + ".log"
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

    #RETRIEVE COMMON REFERENCE
    logger.info('Obtaining common reference')
    common_reference, common_description = find_common_reference(input_dir)
    logger.info('Common complete reference is: ' + common_reference + ': ' + common_description)
    #USE TO EXTRACT FTP PATH
    ftp_path = gcf_to_ftp_path(common_reference)
    #EXTRACT FILES IN THE DESIRED FOLDER
    download_gcf(ftp_path, output_dir, args.all)



if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise