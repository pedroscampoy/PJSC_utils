#!/usr/bin/env python


import os
import pandas as pd
import numpy as np
import logging
import datetime

import re
import sys
import concurrent.futures
import multiprocessing
num_processes = multiprocessing.cpu_count() - 2


def check_create_dir(path):
    # exists = os.path.isfile(path)
    # exists = os.path.isdir(path)
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)


def recalibrate_ddbb_vcf_intermediate(snp_matrix_ddbb_file, bam_folder, min_cov_low_freq=10):

    df_matrix = pd.read_csv(snp_matrix_ddbb_file, sep="\t")

    sample_list_matrix = df_matrix.columns[3:]
    n_samples = len(sample_list_matrix)

    # calculate the chunk size as an integer
    chunk_size = int(df_matrix.shape[0]/num_processes)
    chunks = [df_matrix.iloc[df_matrix.index[i:i + chunk_size]]
              for i in range(0, df_matrix.shape[0], chunk_size)]

    # Iterate over non unanimous positions

    def review_with_vcf(df):
        return_df = df.reset_index(drop=True)
        for index, data_row in df[df.N < n_samples].iloc[:, 3:].iterrows():
            # Extract its position
            whole_position = df.loc[index, "Position"]
            row_reference = whole_position.split('|')[0]
            row_position = int(whole_position.split('|')[2])
            row_alt_snp = whole_position.split('|')[3]

            # Use enumerate to retrieve column index (column ondex + 3)
            # find positions with frequency >80% in mpileup execution
            # Returns ! for coverage 0
            new_presence_row = [recheck_variant_rawvcf_intermediate(
                row_reference, row_position, row_alt_snp, df.columns[n + 3], x, bam_folder, min_cov_low_freq=10) for n, x in enumerate(data_row)]

            new_presence_row_str = [str(x) for x in new_presence_row]

            logger.info(whole_position + " " +
                        (',').join(new_presence_row_str))

            return_df.loc[index, 3:] = new_presence_row

        return return_df

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_processes) as executor:
        futures = []
        for chunk in chunks:
            future = executor.submit(review_with_vcf, chunk)
            futures.append(future)
        for future in concurrent.futures.as_completed(futures):
            logger.info(future.result())

    # pool = multiprocessing.Pool(processes=num_processes)
    # final_df = pd.concat(pool.map(review_with_vcf, chunks))
    # pool.close()
    # pool.join()

    def estract_sample_count(row):
        count_list = [i not in ['!', 0, '0'] for i in row[3:]]
        samples = np.array(df_matrix.columns[3:])
        # samples[np.array(count_list)] filter array with True False array
        return (sum(count_list), (',').join(samples[np.array(count_list)]))

    df_matrix[['N', 'Samples']] = df_matrix.apply(
        estract_sample_count, axis=1, result_type='expand')

    return df_matrix


def recheck_variant_rawvcf_intermediate(reference_id, position, alt_snp, sample, previous_binary, variant_folder, min_cov_low_freq=10):
    """
    CU458896.1	3068036	.	G	A	262.784	.	AB=0.8;ABP=14.7363;AC=1;AF=0.5;AN=2;AO=12;CIGAR=1X;DP=15;DPB=15;DPRA=0;EPP=9.52472;EPPR=3.73412;GTI=0;LEN=1;MEANALT=1;MQM=60;MQMR=60;NS=1;NUMALT=1;ODDS=0.121453;PAIRED=0.333333;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=410;QR=108;RO=3;RPL=7;RPP=3.73412;RPPR=9.52472;RPR=5;RUN=1;SAF=4;SAP=5.9056;SAR=8;SRF=1;SRP=3.73412;SRR=2;TYPE=snp	GT:DP:AD:RO:QR:AO:QA:GL	0/1:15:3,12:3:108:12:410:-32.709,0,-5.55972
    """
    if previous_binary != 0 and previous_binary != '0':
        # logger.info('NON0: {} in sample {} is not 0: {}'.format(position, sample, previous_binary))
        return previous_binary
    else:
        # previous_binary = int(previous_binary)
        position = int(position)

        # Identify correct vcf
        variant_sample_folder = os.path.join(variant_folder, sample)
        for root, _, files in os.walk(variant_sample_folder):
            for name in files:
                if name == "snps.raw.vcf":
                    filename = os.path.join(root, name)

        # Open file and retrieve output

        position_present = False

        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    headers = line.split("\t")
                elif str(position) in line and not line.startswith('#'):
                    line_split = line.split("\t")
                    if line_split[1] == str(position):
                        logger.info(
                            'recalibrating position: {} in sample: {}'.format(position, sample))
                        position_present = True
                        vcf_reference = line_split[0]
                        vcf_position = line_split[1]
                        # vcf_ref_base = line_split[3]
                        vcf_alt_base = line_split[4]
                        params = line_split[-2].split(":")
                        value_params = line_split[-1].split(":")
                        depth_index = params.index('DP')
                        vcf_depth = int(value_params[depth_index])
                        # ref_depth_indef = params.index('RO')
                        alt_depth_indef = params.index('AO')
                        try:
                            vcf_alt_depth = int(value_params[alt_depth_indef])
                        except:
                            vcf_alt_depth = int(
                                value_params[alt_depth_indef].split(',')[-1])
                            vcf_alt_base = vcf_alt_base.split(',')[-1]
                        vcf_alt_freq = vcf_alt_depth/vcf_depth

        if position_present == False:
            logger.debug('Position: {} not present in {}'.format(
                position, filename))
            return 0
        else:
            logger.debug('RECALIBRATE: CHROM: {} POS: {},  ALT: {}, DP: {}, FREQ: {}\nORI==> CHROM: {} POS: {}, ALT: {}'.format(
                vcf_reference, vcf_position, vcf_alt_base, vcf_depth, vcf_alt_freq, reference_id, position, alt_snp))

            if reference_id != vcf_reference:
                logger.info('ERROR: References are different')
                sys.exit(1)
            elif (len(alt_snp) > 1 and str(position) == str(vcf_position)) or (len(vcf_reference) > 1 and str(position) == str(vcf_position)):
                if vcf_depth <= min_cov_low_freq and vcf_depth > 0:
                    logger.debug('Position: {} LOWDEPTH: {}'.format(
                        vcf_position, vcf_alt_freq))
                    return '?'
                else:
                    return vcf_alt_freq
            elif str(position) == str(vcf_position) and alt_snp == vcf_alt_base:
                if vcf_depth <= min_cov_low_freq and vcf_depth > 0:
                    logger.debug('Position: {} LOWDEPTH: {}'.format(
                        vcf_position, vcf_alt_freq))
                    return '?'
                if vcf_alt_freq > 0.1:
                    logger.debug('Position: {} HTZ with freq > 0.1: {}'.format(
                        vcf_position, vcf_alt_freq))
                    return vcf_alt_freq
            else:
                logger.debug('ELSE POS: {} SAMPLE: {}'.format(
                    vcf_position, sample))
                return 0


if __name__ == '__main__':
    # LOGGING
    # Create log file with date and time
    right_now = str(datetime.datetime.now())
    right_now_full = "_".join(right_now.split(" "))
    log_filename = "test" + "_" + right_now_full + ".log"
    log_folder = os.path.join(".", 'Logs')
    check_create_dir(log_folder)
    log_full_path = os.path.join(log_folder, log_filename)

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s:%(message)s')

    file_handler = logging.FileHandler(log_full_path)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    # stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)
    args = sys.argv
    recalibrated_snp_matrix_mpileup = recalibrate_ddbb_vcf_intermediate(
        args[1], args[2], min_cov_low_freq=10)
    recalibrated_snp_matrix_mpileup.to_csv(
        "compare_snp_matrix_recal_mpileup.tsv", sep="\t", index=False)
