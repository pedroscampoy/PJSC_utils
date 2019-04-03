#!/home/laura/env36/bin/python


import os
import argparse
import re
import shutil


def get_arguments():

    parser = argparse.ArgumentParser(prog = 'bowtie_mapper.py', description= 'Creates an index and map reads')
    

    parser.add_argument('-i', '--input', dest="source_dir", type=str, required=True, help='Input directory') #metavar="source_dir"
    parser.add_argument('-o', '--output', dest="output_dir", type=str, required=True, help='Output directory')
    parser.add_argument('-d', '--dict', dest="dictionary", type=str, required=True, help='dictionaru file with paired terms')

    arguments = parser.parse_args()

    return arguments

args = get_arguments()

def import_to_dict(paider_file):
    dictionary_pairs = {}
    with open(paider_file) as dictionary_file:
        for line in dictionary_file:
            line = line.strip()
            pairs = re.split(" +", line)
            name1 = pairs[0]
            name2 = pairs[1]
            dictionary_pairs[name1] = name2
    
    #for k,v in dictionary_pairs.items():
    #    print("%s:%s" % (k,v))

    return dictionary_pairs




def filter_list_re(regex, list_to_filter):
    """
    https://stackoverflow.com/questions/53129958/python-regex-extract-list-elements-each-of-which-matches-multiple-patterns

    To avoid error: 'filter' object is not subscriptable

    https://stackoverflow.com/questions/15876259/typeerror-filter-object-is-not-subscriptable
    """
    regex = re.compile(regex)
    filtered_list = list(filter(lambda x: re.search(regex, x), list_to_filter))
    
    return filtered_list

def rename_fasta_files(dictionary, source_dir, output_dir):
    final_list_to_replace = []
    dictionary = import_to_dict(dictionary)
    full_output_dir = os.path.abspath(output_dir)
    full_source_dir = os.path.abspath(source_dir)
    all_files_in_folder = os.listdir(full_source_dir)
    for name1,name2 in dictionary.items():
                
        R1_source = filter_list_re(name1 + '.*R1.*' + 'fastq.gz$', all_files_in_folder)[0]
        R2_source = filter_list_re(name1 + '.*R2.*' + 'fastq.gz$', all_files_in_folder)[0]

        R1_source_path = os.path.join(full_source_dir, R1_source)
        R2_source_path = os.path.join(full_source_dir, R2_source)

        R1_new = name2 + "_R1.fastq.gz"
        R2_new = name2 + "_R2.fastq.gz"

        R1_new_path = os.path.join(full_output_dir, R1_new)
        R2_new_path = os.path.join(full_output_dir, R2_new)

        final_list_to_replace.append([R1_source_path, R1_new_path])
        final_list_to_replace.append([R2_source_path, R2_new_path])

    return final_list_to_replace
        #print("%s\t%s\n%s\t%s" % (R1_source_path, R1_new_path, R2_source_path, R2_new_path))

def copy_files_paired(list_paired_files):
    """
    For future intention of just renaming"
    """
    for i in list_paired_files:
        shutil.copy2(i[0], i[1])





#Create the list
list_replace = rename_fasta_files(args.dictionary,args.source_dir, args.output_dir)
#Copy files
copy_files_paired(list_replace)

