from __future__ import print_function

import distutils
import shutil

import pandas as pd
import os
import re
from collections import defaultdict

import six
import sys
from Bio import SeqIO

from tracerlib.tracer_func import process_chunk, find_possible_alignments


def makeOutputDir(output_dir_path):
    if not os.path.exists(output_dir_path):
        os.mkdir(output_dir_path)


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def clean_file_list(file_list):
    return_list = []
    trinity_pattern = re.compile(r"(.+)\.Trinity\.fasta")
    for file_name in file_list:
        clean_name = os.path.split(file_name)[1]
        trinity_match = trinity_pattern.search(clean_name)
        if trinity_match:
            clean_name = trinity_match.group(1)
        return_list.append(clean_name)

    return (sorted(return_list))


def get_filename_and_locus(name):
    pattern = re.compile(r"(.+)_TCR([ABDG])")
    pattern_match = pattern.search(name)
    file = pattern_match.group(1)
    locus = pattern_match.group(2)
    return ([file, locus])


def sort_locus_names(dictionary_to_sort):
    for key, value in six.iteritems(dictionary_to_sort):
        sorted_value = sorted(value)
        dictionary_to_sort[key] = sorted_value
    return (dictionary_to_sort)


def load_IMGT_seqs(file):
    seqs = {}
    with open(file, 'rU') as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            seqs[record.id] = str(record.seq)
    return (seqs)


def parse_IgBLAST(locus_names, output_dir, cell_name, imgt_seq_location, species, seq_method, const_seq_file):
    segment_names = ['TRAJ', 'TRAV', 'TRBD', 'TRBJ', 'TRBV']
    IMGT_seqs = dict()
    for segment in segment_names:
        IMGT_seqs[segment] = load_IMGT_seqs("{}/{}.fa".format(imgt_seq_location, segment))

    all_locus_data = defaultdict(dict)
    for locus in locus_names:
        file = "{output_dir}/IgBLAST_output/{cell_name}_{locus}.IgBLASTOut".format(output_dir=output_dir,
                                                                                   cell_name=cell_name, locus=locus)
        if os.path.isfile(file):
            igblast_result_chunks = split_igblast_file(file)

            for chunk in igblast_result_chunks:
                (query_name, chunk_details) = process_chunk(chunk)

                all_locus_data[locus][query_name] = chunk_details
        else:
            all_locus_data[locus] = None

    constant_seqs = pd.read_csv(const_seq_file, index_col=0)['sequence'].to_dict()

    cell = find_possible_alignments(all_locus_data, locus_names, cell_name, IMGT_seqs, output_dir, species, seq_method,
                                    constant_seqs)
    return (cell)


def split_igblast_file(filename):
    # code adapted from http://stackoverflow.com/questions/19575702/pythonhow-to-split-file-into-chunks-by-the-occurrence-of-the-header-word
    token = '# IGBLASTN'
    chunks = []
    current_chunk = []

    with open(filename) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(token) and current_chunk:
                # if line starts with token and the current chunk is not empty
                chunks.append(current_chunk[:])  # add not empty chunk to chunks
                current_chunk = []  # make current chunk blank
            # just append a line to the current chunk on each iteration
            current_chunk.append(line)

        chunks.append(current_chunk)  # append the last chunk outside the loop
    return (chunks)


def check_binary(name, user_path=None):
    if user_path:
        if not is_exe(user_path):
            print("The user provided path for {name} is not executable {user_path}. "
                  "Checking PATH for alternative... ".format(name=name, user_path=user_path))
        else:
            return user_path
    if sys.version_info[0] < 3:
        binary_path = distutils.spawn.find_executable(name)
    else:
        binary_path = shutil.which(name)
    if not binary_path:
        raise OSError("Required binary not find: {name}. Please add to PATH or specify location in config file."
                      .format(name=name))
    else:
        print("Binary for {name} found at {binary_path}.".format(name=name, binary_path=binary_path))
        return binary_path
