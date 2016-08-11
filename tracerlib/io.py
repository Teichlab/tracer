from __future__ import print_function

import csv
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
from tracerlib.core import Invar_cell
import glob
import pdb

import json

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
    name = name.split("_")
    cell = name[0]
    locus = "_".join(name[1:3])
    return ([cell, locus])


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


def parse_IgBLAST(receptor, loci, output_dir, cell_name, raw_seq_dir, species,
                  seq_method, max_junc_len=50, invariant_seqs=None):
    
    IMGT_seqs = dict()
    #expecting_D = dict()
    
    loci_for_segments = defaultdict(list)
    
    #for locus in loci:
    #    expecting_D[locus] = False
    for locus in loci:
        seq_files = glob.glob(os.path.join(raw_seq_dir, "{receptor}_{locus}_*.fa".format(receptor=receptor, 
                                                                                    locus=locus)))
        for f in seq_files:
            #if not f.endswith("_C.fa"):
                segment_name = os.path.splitext(os.path.split(f)[1])[0]
                IMGT_seqs[segment_name] = load_IMGT_seqs(f)
                #if segment_name.split("_")[2] == 'D':
                #    expecting_D[locus] = True
                loci_for_segments[segment_name.split("_")[2]].append(locus)
                    
    #segment_names = ['TRAJ', 'TRAV', 'TRBD', 'TRBJ', 'TRBV']
    #IMGT_seqs = dict()
    #for segment in segment_names:
    #    IMGT_seqs[segment] = load_IMGT_seqs("{}/{}.fa".format(imgt_seq_location, segment))
    
    locus_names = ["_".join([receptor,x]) for x in loci]
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
    cell = find_possible_alignments(all_locus_data, locus_names, cell_name, IMGT_seqs, output_dir, species, seq_method,
                                     invariant_seqs, loci_for_segments, receptor, loci, max_junc_len)
    return (cell)


def split_igblast_file(filename):
    # code adapted from http://stackoverflow.com/questions/19575702/pythonhow-to-split-file-into-chunks-by-the-occurrence-of-the-header-word
    token = '# IGBLASTN'
    chunks = []
    current_chunk = []

    with open(filename) as fh:
        for line in fh:
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


# def check_binary(name, user_path=None):
#     if user_path:
#         if not is_exe(user_path):
#             print("The user provided path for {name} is not executable {user_path}. "
#                   "Checking PATH for alternative... ".format(name=name, user_path=user_path))
#         else:
#             return user_path
#     if sys.version_info[0] < 3:
#         binary_path = distutils.spawn.find_executable(name)
#     else:
#         binary_path = shutil.which(name)
#     if not binary_path:
#     else:
#         print("Binary for {name} found at {binary_path}.".format(name=name, binary_path=binary_path))
#         return binary_path

def check_binary(name, user_path=None):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    if user_path:
        if is_exe(user_path):
            return user_path
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, name)
            if is_exe(exe_file):
                return exe_file

    raise OSError("Required binary not found: {name}. Please add to PATH or specify location in config file."
                  .format(name=name))


def parse_invariant_cells(filename):

    invariant_cells = []
    with open(filename) as fh:
        json_dict = json.load(fh)
        for c in json_dict:
            invariant_cells.append(Invar_cell(c))
    return invariant_cells

def read_colour_file(filename, return_used_list=False, receptor_name=None):
    colour_map = dict()
    used_colours = set()
    with open(filename) as f:
        for line in f:
            line = line.rstrip()
            receptor, locus, prod_colour, nonprod_colour = line.split(",")
            d = {locus : (prod_colour, nonprod_colour)}
            
            if receptor in colour_map:
                colour_map[receptor].update(d)
            else:
                colour_map[receptor] = d
            if receptor_name is not None and receptor == receptor_name:
                used_colours.add(prod_colour)
            elif receptor_name is None:
                used_colours.add(prod_colour)
    if return_used_list:
        t = (colour_map, used_colours)
        return t
    else:
        return colour_map
            
def write_colour_file(filename, colour_map):
    sorted_receptors = sorted(colour_map.keys())
    with open(filename, 'w') as f:
        for receptor in sorted_receptors:
            sorted_loci = sorted(colour_map[receptor].keys())
            for l in sorted_loci:
                colours = colour_map[receptor][l]
                f.write("{},{},{},{}\n".format(receptor, l, colours[0], colours[1]))
