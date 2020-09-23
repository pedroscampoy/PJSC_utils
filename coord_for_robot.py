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

logger = logging.getLogger()

"""
=============================================================
HEADER
=============================================================
FUNCTION: outputs coordinates and sample name from a specific excel template,

INSTITUTION:IISGM
AUTHOR: Pedro J. Sola (pedroscampoy@gmail.com)
d^v^b
VERSION=0.1
CREATED: 21 Sep 2020
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

def extractPlate(df, rowPos1, rowPos2):
    plate = df.iloc[rowPos1:rowPos2,0:13].set_index(['Unnamed: 0'])
    del plate.index.name
    return plate

def extractCoordinatesAndSAmples(plateDf):
    def transformCoordinates(coordinate):
        index = plateDf.index[coordinate[0]]
        return index + str(coordinate[1]+1)

    def obtainSampleName(coordinate):
        return str(int(plateDf.iloc[coordinate[0],coordinate[1]]))

    coordinateList = np.argwhere(plateDf.notnull().values).tolist()
    finalCoordinates = [transformCoordinates(x) for x in coordinateList]

    finalCoordinates = ["'%s'" % (x) for x in finalCoordinates]

    sampleNames = [obtainSampleName(x) for x in coordinateList]

    return finalCoordinates, sampleNames

def estractfilledCoordinates(excellFile):
    df = pd.read_excel(excellFile)
    plate1 = extractPlate(df, 0, 8)
    plate2 = extractPlate(df, 10, 18)
    plate3 = extractPlate(df, 20, 28)
    plate4 = extractPlate(df, 30, 38)
    coor1, samples1 = extractCoordinatesAndSAmples(plate1)
    coor2, samples2 = extractCoordinatesAndSAmples(plate2)
    coor3, samples3 = extractCoordinatesAndSAmples(plate3)
    coor4, samples4 = extractCoordinatesAndSAmples(plate4)

    logger.info("POSITIVE_POS_1=[%s]" % ",".join(coor1))
    logger.info("POSITIVE_POS_2=[%s]" % ",".join(coor2))
    logger.info("POSITIVE_POS_3=[%s]" % ",".join(coor3))
    logger.info("POSITIVE_POS_4=[%s]" % ",".join(coor4))
    logger.info("")
    logger.info("SAMPLES_1 %s" % ",".join(samples1))
    logger.info("SAMPLES_2 %s" % ",".join(samples2))
    logger.info("SAMPLES_3 %s" % ",".join(samples3))
    logger.info("SAMPLES_4 %s" % ",".join(samples4))

    """
    logger.info('POSITIVE_POS_1=',coor1)
    logger.info('POSITIVE_POS_2=',coor2)
    logger.info('POSITIVE_POS_3=',coor3)
    logger.info('POSITIVE_POS_4=',coor4)
    logger.info('')
    logger.info('SAMPLES_1', samples1)
    logger.info('SAMPLES_2', samples2)
    logger.info('SAMPLES_3', samples3)
    logger.info('SAMPLES_4', samples4)
    """


def main():

    def get_arguments():

        parser = argparse.ArgumentParser(prog = 'common_mash_reference.py', description= 'Search for all mash files and find the representative reference')
        
        parser.add_argument('-i', '--input', dest="input_file", metavar="input_directory", type=str, required=True, help='REQUIRED.Input Fexcel teplate')
        parser.add_argument('-o', '--output', type=str, required=False, default=False, help='Output directory to extract log')

        arguments = parser.parse_args()

        return arguments

    args = get_arguments()

    input_file = os.path.abspath(args.input_file)

    if args.output == False:
        output_dir = ('/').join(input_file.split('/')[0:-1])
    else:
        output_dir = os.path.abspath(args.output)

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
    #CALCULATE COORD
    logger.info(WHITE_BG + 'Obtaining coordinates' + END_FORMATTING)
    estractfilledCoordinates(input_file)
    logger.info(WHITE_BG + 'DONE' + END_FORMATTING)

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logger.exception(e)
        raise