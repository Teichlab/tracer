import os
import re
from collections import defaultdict

import six
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
    for record in SeqIO.parse(open(file, 'rU'), 'fasta'):
        seqs[record.id] = str(record.seq)
    return (seqs)


def parse_IgBLAST(locus_names, output_dir, cell_name, imgt_seq_location, species, seq_method):
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

    cell = find_possible_alignments(all_locus_data, locus_names, cell_name, IMGT_seqs, output_dir, species, seq_method)
    return (cell)


def split_igblast_file(filename):
    # code adapted from http://stackoverflow.com/questions/19575702/pythonhow-to-split-file-into-chunks-by-the-occurrence-of-the-header-word
    token = '# IGBLASTN'
    chunks = []
    current_chunk = []

    for line in open(filename):
        line = line.rstrip()
        if line.startswith(token) and current_chunk and not line.startswith("Total queries"):
            # if line starts with token and the current chunk is not empty
            chunks.append(current_chunk[:])  # add not empty chunk to chunks
            current_chunk = []  # make current chunk blank
        # just append a line to the current chunk on each iteration
        if not line.startswith("Total queries"):
            current_chunk.append(line)

    chunks.append(current_chunk)  # append the last chunk outside the loop
    return (chunks)