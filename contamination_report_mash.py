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

pd.set_option('display.max_colwidth', None)

logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
FUNCTION: search for all mash files (ends with .winner.tab), find and count the representative reference,
    and look for contamination, several species in winner mode
INSTITUTION:CNM-ISCIII
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
d^v^b
VERSION=0.1
CREATED: 23 March 2020
REVISION: 

TODO:

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

def return_best_matches(screen_file):
    df = pd.read_csv(screen_file, sep='\t', names=['identity', 'shared-hashes', 'median-multiplicity', 'p-value', 'query-ID', 'query-comment'])
    if df.shape[1] == 6:
        df = df[~df['query-comment'].str.contains(' phage')]
        df = df[~df['query-comment'].str.contains(' Phage')]
        df = df[~df['query-comment'].str.contains('Bacteriophage')]
        df = df.sort_values(by=['identity'], ascending=False)
        df.reset_index(inplace=True)
        return df[df.identity > 0.9]
    else:
        logger.info("Wrong mash file:")
        logger.info(screen_file)

def extract_species(df_line):
    split_query = df_line['query-comment'].split(' ')
    if re.match(r'^\[.{4,10}\]', df_line['query-comment']):
        return ' '.join(split_query[3:5])
    else:
        return ' '.join(split_query[1:3])

def report_contamination(df_mash):
    '''
    Dependencies:   -df_mash
                    -extract_species(df_line)
    '''
    try:
        main_species = extract_species(df_mash.iloc[0])
    except:
        main_species = "NO REPRESENTATIVE ESPECIES"
    contamination_report = "Main species: " + main_species + "\n"
    if df_mash.shape[0] > 1:
        for _, data_row in df_mash.iloc[1:].iterrows():
            hashes = data_row['shared-hashes'] #.split('/')[0]
            contamined_species = extract_species(data_row)
            if contamined_species != main_species:
                contamination_report = contamination_report + 'Contamination: ' + str(hashes) + " " + contamined_species + "\n"
    return contamination_report

def find_contamination_mash(folder, output_file):
    '''
    Dependencies:   -input folder
                    -def return_best_matches(screen_file)
                    -report_contamination(df_mash)
    '''
    with open(output_file, 'w+') as outf:
        #Find common reference from mash result
        for root, _, files in os.walk(folder):
            for name in files:
                if name.endswith(".winner.tab"):
                    filename = os.path.join(root, name)
                    sample = name.split('.')[0]
                    logger.info(sample)
                    best_matches = return_best_matches(filename)
                    outf.write(sample + "\n" + report_contamination(best_matches))
def main():

    def get_arguments():

        parser = argparse.ArgumentParser(prog = 'common_mash_reference.py', description= 'Search for all mash files and find the representative reference')
        
        parser.add_argument('-i', '--input', dest="input_dir", metavar="input_directory", type=str, required=True, help='REQUIRED.Input directory containing all screen.tab files')
        parser.add_argument('-o', '--output', type=str, required=False, help='REQUIRED. Output directory to extract all results')

        arguments = parser.parse_args()

        return arguments

    args = get_arguments()
    

    #output_dir = os.path.abspath(args.output)
    input_dir = os.path.abspath(args.input_dir)
    
    #check_create_dir(output_dir)

    #LOGGING
    #Create log file with date and time
    right_now = str(datetime.date.today())
    right_now_full = "_".join(right_now.split(" "))

    log_filename = 'mash_contamination' + "_" + right_now_full + ".log"
    log_full_path = os.path.join(input_dir, log_filename)

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
    if args.output == None:
        output_file = os.path.join(input_dir, 'contamination_report.txt')
    else:
        output_file = args.output

    logger.info('Finding contaminants')
    find_contamination_mash(input_dir, output_file)
    logger.info("DONE: file can be found in " + output_file)


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise