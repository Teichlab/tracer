from __future__ import print_function

import csv

import six
import matplotlib as mpl
import pandas as pd
import numpy as np
from pprint import pprint

from tracerlib import io, core
from tracerlib.io import check_binary

mpl.use('pdf')
import re
import seaborn as sns
from matplotlib import pyplot as plt
from tracerlib import base_dir
from tracerlib import tracer_func
try:
    from configparser import ConfigParser, NoOptionError
except ImportError:
    from ConfigParser import ConfigParser, NoOptionError
import argparse
import sys
import os
import subprocess
import glob
import shutil
from collections import defaultdict, Counter
from time import sleep
import warnings
import pickle
from prettytable import PrettyTable
from Bio.Seq import Seq
from Bio import SeqIO
import itertools
import pdb
from numpy import percentile, array
from matplotlib.colors import hex2color, rgb2hex
import random
import copy
import colorsys


class TracerTask(object):
    base_parser = argparse.ArgumentParser(add_help=False,
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    base_parser.add_argument('--ncores', '-p', metavar="<CORES>",
                             help='number of processor cores to use', type=int,
                             default=1)
    base_parser.add_argument('--config_file', '-c', metavar="<CONFIG_FILE>",
                             help='config file to use',
                             default=None)
    base_parser.add_argument('--resource_dir', metavar="<RESOURCE_DIR>",
                             help='root directory for resources', default=None)

    config = None

    def run(self):
        pass

    def get_binary(self, name):
        tool_key = name.lower() + '_path'
        user_path = None
        if self.config.has_option('tool_locations', tool_key):
            user_path = self.resolve_relative_path(
                self.config.get('tool_locations', tool_key))
        return check_binary(name, user_path)

    def get_tracer_path(self):
        tracer_path = None
        if self.config.has_option('tracer_location', 'tracer_path'):
            path = self.config.get('tracer_location', 'tracer_path')
            if os.path.exists(path):
                tracer_path = path
            else:
                print("Please specify the path to where you originally"
                    " downloaded TraCeR in the config file.")
        return tracer_path


    def read_config(self, config_file):
        ###first check if tracer is run within the docker container
        if os.path.isfile("/tracer/docker_helper_files/docker_tracer.conf") and not config_file:
            config_file="/tracer/docker_helper_files/docker_tracer.conf"
            print("Using default configuration file in Docker container...")
        elif not config_file:
            config_file = os.path.expanduser('~/.tracerrc')
            if not os.path.isfile(config_file):
                print("Config file not found at ~/.tracerrc."
                      " Using default tracer.conf in repo...")
                tracer_path = self.get_tracer_path()
                config_file = os.path.join(tracer_path, 'tracer.conf')
                if not os.path.isfile(config_file):
                    config_file = os.path.join(base_dir, 'tracer.conf')
        # then look for environmental variable
        elif not config_file and :
            config_file = os.environ.get('TRACER_CONF', None)
            if config_file is not None:
                config_file = os.path.expanduser(config_file)
                if not os.path.isfile(config_file):
                    config_file = None
        # Then check the default location
        
        tracer_func.check_config_file(config_file)
        config = ConfigParser()
        config.read(config_file)

        return config

    def resolve_relative_path(self, path):
        if not path.startswith("/"):
            base_directory = os.path.abspath(os.path.dirname(__file__))
            full_path = os.path.normpath(
                "/{}/../{}".format(base_directory, path))
        else:
            full_path = path
        return full_path

    def print_cell_summary(self, cell, output_file, receptor_name, loci):
        out_file = open(output_file, 'w')
        out_file.write(
            '------------------\n{name}\n------------------\n'.format(
                name=cell.name))

        # summarise the productive/total recombinants
        for l in loci:
            out_file.write(
                '{receptor}_{locus} recombinants: {summary}\n'.format(
                    receptor=receptor_name, locus=l,
                    summary=cell.summarise_productivity(receptor_name, l)))

        out_file.write('\n\n')

        for l in loci:
            out_file.write(
                "#{receptor}_{locus}#\n".format(receptor=receptor_name,
                                                locus=l))
            rs = cell.recombinants[receptor_name][l]
            if rs is None:
                out_file.write(
                    "No {receptor}_{locus} recombinants found\n\n".format(
                        receptor=receptor_name, locus=l))
            else:
                for r in rs:
                    out_file.write(r.get_summary())
                    out_file.write("\n\n")

        # out_file.write('#TCRA#\n')
        # if cell.A_recombinants is None:
        #    out_file.write("No TCRA recombinants found\n\n")
        # else:
        #    for rec in cell.A_recombinants:
        #        out_file.write(rec.get_summary())
        #        out_file.write("\n\n")
        # out_file.write('#TCRB#\n')
        #
        # if cell.B_recombinants is None:
        #    out_file.write("No TCRB recombinants found\n\n")
        # else:
        #    for rec in cell.B_recombinants:
        #        out_file.write(rec.get_summary())
        #        out_file.write("\n\n")
        out_file.close()

    def die_with_empty_cell(self, cell_name, output_dir, species):
        print("##No recombinants found##")
        cell = core.Cell(cell_name, None, is_empty=True, species=species,
                         receptor=self.receptor_name, loci=self.loci)

        self.print_cell_summary(
            cell,
            "{output_dir}/unfiltered_{receptor}_seqs/unfiltered_{receptor}s.txt".format(
                output_dir=self.output_dir,
                receptor=self.receptor_name),
            self.receptor_name, self.loci)

        # Save cell in a pickle
        with open(
                "{output_dir}/unfiltered_{receptor}_seqs/{cell_name}.pkl".format(
                    output_dir=self.output_dir,
                    cell_name=cell.name,
                    receptor=self.receptor_name), 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)
        cell.filter_recombinants()
        self.print_cell_summary(
            cell,
            "{output_dir}/filtered_{receptor}_seqs/filtered_{receptor}s.txt".format(
                output_dir=self.output_dir,
                receptor=self.receptor_name),
            self.receptor_name, self.loci)

        with open(
                "{output_dir}/filtered_{receptor}_seqs/{cell_name}.pkl".format(
                    output_dir=self.output_dir,
                    cell_name=cell.name,
                    receptor=self.receptor_name), 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)

        exit(0)

    def get_species_root(self, species, root=None, build_mode = False):
        if root is None:
            if os.path.isfile("/tracer/docker_helper_files/docker_tracer.conf") and tracer_path is None:
                tracer_path="/tracer"
            else:
            tracer_path = self.get_tracer_path()
            if tracer_path is not None:
                resources_root = os.path.join(tracer_path, 'resources', species)
            else: 
                # Look for resources in base directory if tracer_path is not specified
                resources_root = os.path.join(base_dir, 'resources', species)
        else:
            resources_root = os.path.join(root, species)
        if not build_mode:
            assert os.path.isdir(resources_root), "Species not found in resources"
        return (resources_root)


    # def get_available_species(self, root=None):
    #    if root is None:
    #        resources_dir = os.path.join(base_dir, 'resources')
    #    else:
    #        resources_dir = root
    #    species_dirs = next(os.walk(resources_dir))[1]
    #    return (species_dirs)


class Assembler(TracerTask):
    def __init__(self, **kwargs):
        if not kwargs:

            # get list of all available species in resources

            parser = argparse.ArgumentParser(
                description="Reconstruct TCR sequences from RNAseq reads for a single cell",
                parents=[self.base_parser],
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
            parser.add_argument('--resume_with_existing_files', '-r',
                                help='look for existing intermediate files and use those instead of starting from scratch',
                                action="store_true")
            parser.add_argument('--species', '-s',
                                help='Species to use for reconstruction',
                                default='Mmus')
            parser.add_argument('--receptor_name',
                                help="Name of receptor to reconstruct",
                                default='TCR')
            parser.add_argument('--loci',
                                help="Space-separated list of loci to reconstruct for receptor",
                                default=['A', 'B'], nargs='+')
            parser.add_argument('--seq_method', '-m',
                                help='Method for constructing sequence to assess productivity, \
                                quantify expression and for output reporting. See README for details.',
                                choices=['imgt', 'assembly'], default='imgt')
            parser.add_argument('--quant_method', '-q',
                                help='Method for expression quantification. See README for details.',
                                choices=['kallisto', 'salmon'],
                                default='kallisto')
            parser.add_argument('--small_index',
                                help='Set this to speed up expression quantification by using a smaller index. See README for details.',
                                action="store_true")
            parser.add_argument('--single_end',
                                help='set this if your sequencing data are single-end reads',
                                action="store_true")
            parser.add_argument('--fragment_length',
                                help='Estimated average fragment length in the sequencing library.'
                                     ' Used for Kallisto quantification. REQUIRED for single-end data.',
                                default=False)
            parser.add_argument('--fragment_sd',
                                help='Estimated standard deviation of average fragment length in the sequencing library.'
                                     ' Used for Kallisto quantification. REQUIRED for single-end data.',
                                default=False)
            parser.add_argument('--max_junc_len',
                                help="Maximum permitted length of junction string in recombinant identifier. "
                                     "Used to filter out artefacts. May need to be longer for TCRdelta.",
                                default=50, type=int)

            parser.add_argument('fastq1', metavar="<FASTQ1>",
                                help='first fastq file')
            parser.add_argument('fastq2', metavar="<FASTQ2>",
                                help='second fastq file', nargs='?')
            parser.add_argument('cell_name', metavar="<CELL_NAME>",
                                help='name of cell for file labels')
            parser.add_argument('output_dir', metavar="<OUTPUT_DIR>",
                                help='directory for output as <output_dir>/<cell_name>')

            args = parser.parse_args(sys.argv[2:])

            resource_dir = args.resource_dir
            self.cell_name = args.cell_name
            self.fastq1 = args.fastq1
            self.single_end = args.single_end
            self.fastq2 = args.fastq2
            self.ncores = str(args.ncores)
            self.species = args.species
            self.seq_method = args.seq_method
            self.quant_method = args.quant_method
            self.small_index = args.small_index
            self.resume_with_existing_files = args.resume_with_existing_files
            self.fragment_length = args.fragment_length
            self.fragment_sd = args.fragment_sd
            self.output_dir = args.output_dir
            self.receptor_name = args.receptor_name
            self.loci = args.loci
            self.max_junc_len = args.max_junc_len
            config_file = args.config_file

        else:
            resource_dir = kwargs.get('resource_dir')
            self.cell_name = kwargs.get('cell_name')
            self.fastq1 = kwargs.get('fastq1')
            self.fastq2 = kwargs.get('fastq2')
            self.ncores = kwargs.get('ncores')
            self.species = kwargs.get('species')
            self.seq_method = kwargs.get('seq_method')
            self.quant_method = kwargs.get('quant_method')
            self.small_index = kwargs.get('small_index')
            self.resume_with_existing_files = kwargs.get(
                'resume_with_existing_files')
            self.output_dir = kwargs.get('output_dir')
            self.single_end = kwargs.get('single_end')
            self.fragment_length = kwargs.get('fragment_length')
            self.fragment_sd = kwargs.get('fragment_sd')
            self.receptor_name = kwargs.get('receptor_name')
            self.loci = kwargs.get('loci')
            self.max_junc_len = kwargs.get('max_junc_len')
            config_file = kwargs.get('config_file')

        self.config = self.read_config(config_file)

        self.species_root = self.get_species_root(self.species,
                                                  root=resource_dir)
        # self.locus_names = ["TCRA", "TCRB"]

        # Check the fastq config is correct
        if not self.single_end:
            assert self.fastq2, "Only one fastq file specified. Either set --single_end or provide second fastq."
        else:
            self.fastq2 = None
            if self.fastq2:
                print(
                    "Two fastq files given with --single-end option. Ignoring second file.")
            assert self.fragment_length and self.fragment_sd, \
                'Must specify estimated average fragment length (--fragment_length)' \
                ' and standard deviation (--fragment_sd) for use with single-end data'
            assert self.fragment_length, \
                'Must specify estimated average fragment length (--fragment_length) for use with single-end data'
            assert self.fragment_sd, \
                'Must specify estimated fragment length standard deviation (--fragment_sd) for use with single-end data'

        # Check FASTQ files exist
        if not os.path.isfile(self.fastq1):
            raise OSError('2', 'FASTQ file not found', self.fastq1)
        if not self.single_end and self.fastq2:
            if not os.path.isfile(self.fastq2):
                raise OSError('2', 'FASTQ file not found', self.fastq2)

    def run(self, **kwargs):

        # Set-up output directories
        root_output_dir = os.path.abspath(self.output_dir)
        io.makeOutputDir(root_output_dir)
        self.output_dir = root_output_dir + "/" + self.cell_name

        io.makeOutputDir(self.output_dir)

        data_dirs = ['aligned_reads', 'Trinity_output', 'IgBLAST_output',
                     'unfiltered_{receptor}_seqs'.format(
                         receptor=self.receptor_name),
                     'expression_quantification',
                     'filtered_{receptor}_seqs'.format(
                         receptor=self.receptor_name)]
        for d in data_dirs:
            io.makeOutputDir("{}/{}".format(self.output_dir, d))

        # Perform TraCeR's core functions
        self.align()
        self.de_novo_assemble()
        cell = self.ig_blast()
        self.quantify(cell)

        fasta_filename = "{output_dir}/unfiltered_{receptor}_seqs/{cell_name}_{receptor}seqs.fa".format(
            output_dir=self.output_dir,
            cell_name=self.cell_name,
            receptor=self.receptor_name)
        fasta_file = open(fasta_filename, 'w')
        fasta_file.write(cell.get_fasta_string())
        fasta_file.close()

        self.print_cell_summary(
            cell,
            "{output_dir}/unfiltered_{receptor}_seqs/unfiltered_{receptor}s.txt".format(
                output_dir=self.output_dir,
                receptor=self.receptor_name),
            self.receptor_name, self.loci)

        # Save cell in a pickle
        with open(
                "{output_dir}/unfiltered_{receptor}_seqs/{cell_name}.pkl".format(
                    output_dir=self.output_dir,
                    cell_name=cell.name,
                    receptor=self.receptor_name), 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)
        print("##Filtering by read count##")
        cell.filter_recombinants()
        fasta_filename = "{output_dir}/filtered_{receptor}_seqs/{cell_name}_{receptor}seqs.fa".format(
            output_dir=self.output_dir,
            cell_name=self.cell_name,
            receptor=self.receptor_name)
        fasta_file = open(fasta_filename, 'w')
        fasta_file.write(cell.get_fasta_string())
        fasta_file.close()
        self.print_cell_summary(
            cell,
            "{output_dir}/filtered_{receptor}_seqs/filtered_{receptor}s.txt".format(
                output_dir=self.output_dir,
                receptor=self.receptor_name),
            self.receptor_name, self.loci)

        with open(
                "{output_dir}/filtered_{receptor}_seqs/{cell_name}.pkl".format(
                    output_dir=self.output_dir,
                    cell_name=cell.name,
                    receptor=self.receptor_name), 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)

    def align(self):
        bowtie2 = self.get_binary('bowtie2')

        synthetic_genome_path = os.path.join(self.species_root,
                                             'combinatorial_recombinomes')
        # Align with bowtie
        tracer_func.bowtie2_alignment(
            bowtie2, self.ncores, self.receptor_name, self.loci,
            self.output_dir, self.cell_name,
            synthetic_genome_path, self.fastq1, self.fastq2,
            self.resume_with_existing_files, self.single_end)
        print()

    def de_novo_assemble(self):

        try:
            trinity = self.get_binary('trinity')
        except OSError:
            trinity = self.get_binary('Trinity')

        # Trinity version
        if not self.config.has_option('trinity_options', 'trinity_version'):
            try:
                subprocess.check_output([trinity, '--version'])
            except subprocess.CalledProcessError as err:
                if re.search('v2', err.output.decode('utf-8')):
                    self.config.set('trinity_options', 'trinity_version', '2')
                else:
                    self.config.set('trinity_options', 'trinity_version', '1')

        if self.config.has_option('trinity_options', 'trinity_grid_conf'):
            trinity_grid_conf = self.resolve_relative_path(
                self.config.get('trinity_options', 'trinity_grid_conf'))
        else:
            trinity_grid_conf = False
            
        # Is Trinity version compatible with --no_normalize_reads argument
        no_normalise = False
        trinity_version = self.config.get('trinity_options', 'trinity_version')
        if trinity_version == '2':
            try:
                subprocess.check_output([trinity, '--version'])
            except subprocess.CalledProcessError as err:
                first_line = err.output.decode('utf-8').split("\n")[0]
                m = re.search(r'v(\d+\.\d+\.?\d*)', first_line)
                if m is not None:
                    minor_version = int(m.group()[1:].split(".")[1])
                    if minor_version >= 3:
                        no_normalise = True
                    else:
                        no_normalise = False
        

        # De novo assembly with trinity
        trinity_JM = self.config.get('trinity_options', 'max_jellyfish_memory')
        trinity_version = self.config.get('trinity_options', 'trinity_version')
        successful_files = tracer_func.assemble_with_trinity(
            trinity, self.receptor_name, self.loci, self.output_dir,
            self.cell_name, self.ncores, trinity_grid_conf,
            trinity_JM, trinity_version, self.resume_with_existing_files,
            self.single_end, self.species, no_normalise, self.config)
        if len(successful_files) == 0:
            print("No successful Trinity assemblies")
            self.die_with_empty_cell(self.cell_name, self.output_dir,
                                     self.species)

        print()

    def ig_blast(self):
        igblastn = self.get_binary('igblastn')

        # Reference data locations
        igblast_index_location = os.path.join(self.species_root, 'igblast_dbs')
        imgt_seq_location = os.path.join(self.species_root, 'raw_seqs')

        igblast_seqtype = self.config.get('IgBlast_options', 'igblast_seqtype')

        # IgBlast of assembled contigs - run twice. Once with output format 3 and once with output format 7
        for fmt in (str(3),str(7)):
            tracer_func.run_IgBlast(igblastn, self.receptor_name, self.loci,
                                self.output_dir, self.cell_name,
                                igblast_index_location,
                                igblast_seqtype, self.species,
                                self.resume_with_existing_files, fmt)
        
        print()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # cell = io.parse_IgBLAST(self.receptor_name, self.loci, self.output_dir, self.cell_name, imgt_seq_location,
            # self.species, self.seq_method, self.invariant_sequences)
            cell = io.parse_IgBLAST(self.receptor_name, self.loci,
                                    self.output_dir, self.cell_name,
                                    imgt_seq_location,
                                    self.species, self.seq_method,
                                    self.max_junc_len)
            if cell.is_empty:
                self.die_with_empty_cell(self.cell_name, self.output_dir,
                                         self.species)

        return cell

    def quantify(self, cell):

        if not self.config.has_option('base_transcriptomes', self.species):
            raise OSError("No transcriptome reference specified for {species}."
                          " Please specify location in config file."
                          .format(species=self.species))
        else:
            base_transcriptome = self.resolve_relative_path(
                self.config.get('base_transcriptomes', self.species))
            if not os.path.isfile(base_transcriptome):
                raise OSError(
                    '2', 'Transcriptome reference not found', base_transcriptome)

        # set up Salmon/Kallisto parameters

        if self.quant_method == 'salmon':
            salmon = self.get_binary('salmon')
            if self.config.has_option('salmon_options', 'libType'):
                salmon_libType = self.config.get('salmon_options', 'libType')
            else:
                print("No library type specified for salmon "
                      "in the configuration file. Using automatic detection "
                      "(--libType A).")
                salmon_libType = 'A'
            if self.config.has_option('salmon_options', 'kmerLen'):
                salmon_kmerLen = self.config.get('salmon_options', 'kmerLen')
            else:
                print("No kmer length specified for salmon "
                      "in the configuration file. Using default value of 31.")
                salmon_kmerLen = 31

        else:
            kallisto = self.get_binary('kallisto')

        # small-index method: filter base_transcriptome -> ref_transcriptome

        if self.small_index:

            if not self.config.has_option(
                    "{}{}".format(self.quant_method, '_base_indices'),
                    self.species):
                raise OSError(
                    "No {method} {species} reference index specified to be used"
                    " with option --small_index. "
                    "Please specify location in config file.".format(
                        method=self.quant_method, species=self.species))
            else:
                base_index = self.resolve_relative_path(
                    self.config.get(
                        "{}{}".format(self.quant_method, '_base_indices'),
                        self.species))
                if not os.path.exists(base_index):
                    raise OSError('2', 'Reference index not found', base_index)

            smind_dir = "{}/{}".format(self.output_dir,
                                       'expression_quantification_from_base_index')
            new_ref_transcriptome = "{}/{}".format(smind_dir, 'newref.fa')
            io.makeOutputDir(smind_dir)

            if self.quant_method == 'salmon':
                tracer_func.quantify_with_salmon_from_index(
                    salmon, cell, smind_dir, self.cell_name, base_index,
                    self.fastq1, self.fastq2, self.ncores,
                    self.resume_with_existing_files, self.single_end,
                    self.fragment_length, self.fragment_sd,
                    salmon_libType, salmon_kmerLen)

                quantfile = "{}/{}".format(smind_dir, 'quant.sf')
                tpmcol = 3  # TPM col = 3 in quant.sf

            else:  # use kallisto
                tracer_func.quantify_with_kallisto_from_index(
                    kallisto, cell, smind_dir, self.cell_name, base_index,
                    self.fastq1, self.fastq2, self.ncores,
                    self.resume_with_existing_files, self.single_end,
                    self.fragment_length, self.fragment_sd)

                quantfile = "{}/{}".format(smind_dir, 'abundance.tsv')
                tpmcol = 4  # TPM col = 4 in abundance.tsv

            tracer_func.extract_newref_from_quant(
                base_transcriptome, quantfile, tpmcol, new_ref_transcriptome)
            ref_transcriptome = new_ref_transcriptome

        else:
            # use entire base_transcriptome for index construction as usual
            ref_transcriptome = base_transcriptome

        # actual quantification of TCR Seq

        if self.quant_method == 'salmon':
            tracer_func.quantify_with_salmon(
                salmon, cell, self.output_dir, self.cell_name, ref_transcriptome,
                self.fastq1, self.fastq2, self.ncores,
                self.resume_with_existing_files,
                self.single_end, self.fragment_length,
                self.fragment_sd, salmon_libType, salmon_kmerLen)
            print()
            counts = tracer_func.load_salmon_counts(
                "{}/expression_quantification/quant.sf".format(self.output_dir))
        else:
            tracer_func.quantify_with_kallisto(
                kallisto, cell, self.output_dir, self.cell_name,
                ref_transcriptome, self.fastq1, self.fastq2, self.ncores,
                self.resume_with_existing_files, self.single_end,
                self.fragment_length, self.fragment_sd)
            print()
            counts = tracer_func.load_kallisto_counts(
                "{}/expression_quantification/abundance.tsv".format(self.output_dir))

        if self.small_index:
            os.remove(new_ref_transcriptome)  # remove filtered reference file

        for receptor, locus_dict in six.iteritems(cell.recombinants):
            for locus, recombinants in six.iteritems(locus_dict):
                if recombinants is not None:
                    for rec in recombinants:
                        tpm = counts[receptor][locus][rec.contig_name]
                        rec.TPM = tpm


class Summariser(TracerTask):
    def __init__(self, **kwargs):

        if not kwargs:
            parser = argparse.ArgumentParser(
                description="Summarise set of cells with reconstructed TCR sequences",
                parents=[self.base_parser])
            parser.add_argument('--species', '-s',
                                help='Species to use for reconstruction',
                                default='Mmus')
            parser.add_argument('--receptor_name',
                                help="Name of receptor to summarise",
                                default='TCR')
            parser.add_argument('--loci',
                                help="Space-separated list of loci to summarise for receptor",
                                default=['A', 'B'], nargs='+')
            parser.add_argument('--use_unfiltered', '-u',
                                help='use unfiltered recombinants',
                                action="store_true")
            parser.add_argument('--keep_invariant', '-i',
                                help='ignore invariant cells when constructing networks',
                                action="store_true")
            parser.add_argument('--graph_format', '-f',
                                metavar="<GRAPH_FORMAT>",
                                help='graphviz output format [pdf]',
                                default='pdf')
            parser.add_argument('--no_networks',
                                help='skip attempts to draw network graphs',
                                action="store_true")
            parser.add_argument('dir', metavar="<DIR>",
                                help='directory containing subdirectories for each cell to be summarised')
            args = parser.parse_args(sys.argv[2:])

            resource_dir = args.resource_dir
            self.root_dir = os.path.abspath(args.dir)
            self.graph_format = args.graph_format
            self.keep_invariant = args.keep_invariant
            self.use_unfiltered = args.use_unfiltered
            self.draw_graphs = not args.no_networks
            self.receptor_name = args.receptor_name
            self.loci = args.loci
            self.species = args.species
            config_file = args.config_file
        else:
            resource_dir = kwargs.get('resource_dir')
            self.use_unfiltered = kwargs.get('use_unfiltered')
            self.root_dir = os.path.abspath(kwargs.get('root_dir'))
            self.draw_graphs = not (kwargs.get('no_networks'))
            self.graph_format = kwargs.get('graph_format')
            self.keep_invariant = kwargs.get('keep_invariant')
            self.receptor_name = kwargs.get('receptor_name')
            self.loci = kwargs.get('loci')
            self.species = kwargs.get('species')
            config_file = kwargs.get('config_file')

        # Read config file
        self.config = self.read_config(config_file)
        self.species_dir = self.get_species_root(self.species,
                                                 root=resource_dir)

        invariant_cells = os.path.join(self.species_dir, 'invariant_cells.json')
        if os.path.isfile(invariant_cells):
            self.invariant_cells = io.parse_invariant_cells(invariant_cells)
        else:
            self.invariant_cells = None

    def run(self):

        if self.draw_graphs:
            dot = self.get_binary('dot')
            neato = self.get_binary('neato')
            #dot = self.resolve_relative_path(
            #    self.config.get('tool_locations', 'dot_path'))
            #neato = self.resolve_relative_path(
            #    self.config.get('tool_locations', 'neato_path'))
            #
            ## check that executables from config file can be used
            #not_executable = []
            #for name, x in six.iteritems({"dot": dot, "neato": neato}):
            #    if not io.is_exe(x):
            #        not_executable.append((name, x))
            #if len(not_executable) > 0:
            #    print()
            #    print("Could not execute the following required tools."
            #          " Check your configuration file.")
            #    for t in not_executable:
            #        print(t[0], t[1])
            #    print()
            #    exit(1)
        else:
            dot = ""
            neato = ""

        cells = {}
        empty_cells = []
        subdirectories = next(os.walk(self.root_dir))[1]

        if self.use_unfiltered:
            pkl_dir = "unfiltered_{}_seqs".format(self.receptor_name)
            outdir = "{}/unfiltered_{}_summary".format(
                self.root_dir, self.receptor_name + "".join(self.loci))
            # outfile = open("{root_dir}/unfiltered_TCR_summary.txt".format(root_dir=root_dir), 'w')
            # length_filename_root = "{}/unfiltered_reconstructed_lengths_TCR".format(root_dir)

        else:
            pkl_dir = "filtered_{}_seqs".format(self.receptor_name)
            outdir = "{}/filtered_{}_summary".format(
                self.root_dir, self.receptor_name + "".join(self.loci))
            # outfile = open("{root_dir}/filtered_TCR_summary.txt".format(root_dir=root_dir), 'w')
            # length_filename_root = "{}/filtered_reconstructed_lengths_TCR".format(root_dir)

        io.makeOutputDir(outdir)

        outfile = open("{}/{}_summary.txt".format(outdir, self.receptor_name),
                       'w')
        length_filename_root = "{}/reconstructed_lengths_{}".format(outdir,
                                                                    self.receptor_name)

        for d in subdirectories:
            cell_pkl = "{root_dir}/{d}/{pkl_dir}/{d}.pkl".format(
                pkl_dir=pkl_dir, d=d, root_dir=self.root_dir)
            if os.path.isfile(cell_pkl):
                with open(cell_pkl, 'rb') as pkl:
                    cl = pickle.load(pkl)
                cells[d] = cl
                if cl.is_empty or cl.missing_loci_of_interest(
                        self.receptor_name, self.loci):
                    empty_cells.append(d)

        cell_recovery_count = dict()
        # count cells with productive chains for each locus and for each
        # possible pair
        for l in self.loci:
            cell_recovery_count[l] = 0

        possible_pairs = ["".join(x) for x in
                          itertools.combinations(self.loci, 2)]

        for p in possible_pairs:
            cell_recovery_count[p] = 0

        for cell_name, cell in six.iteritems(cells):
            prod_counts = dict()
            for l in self.loci:
                prod_counts[l] = cell.count_productive_recombinants(
                    self.receptor_name, l)
                if prod_counts[l] > 0:
                    cell_recovery_count[l] += 1

            for p in possible_pairs:
                if prod_counts[p[0]] > 0 and prod_counts[p[1]] > 0:
                    cell_recovery_count[p] += 1

        total_cells = len(cells)

        for l in self.loci:
            count = cell_recovery_count[l]
            pc = round((count / float(total_cells)) * 100, 1)
            outfile.write(
                "{receptor}_{locus} reconstruction:\t{count} / {total} ({pc}%)\n".format(
                    receptor=self.receptor_name,
                    locus=l, count=count,
                    total=total_cells, pc=pc))
        outfile.write("\n")

        for p in possible_pairs:
            count = cell_recovery_count[p]
            pc = round((count / float(total_cells)) * 100, 1)
            outfile.write(
                "{p} productive reconstruction:\t{count} / {total} ({pc}%)\n".format(
                    p=p,
                    count=count,
                    total=total_cells, pc=pc))
        outfile.write("\n")

        all_counters = defaultdict(Counter)
        prod_counters = defaultdict(Counter)

        for cell in cells.values():
            for l in self.loci:
                all_counters[l].update(
                    {cell.count_total_recombinants(self.receptor_name, l): 1})
                prod_counters[l].update({cell.count_productive_recombinants(
                    self.receptor_name, l): 1})

        all_recombinant_counts = []
        for locus in all_counters:
            all_recombinant_counts = all_recombinant_counts + \
                list(all_counters[locus].keys())
        max_recombinant_count = max(all_recombinant_counts)

        # max_recombinant_count = max(list(counters['all_alpha'].keys()) + list(counters['all_beta'].keys()))
        table_header = ['', '0 recombinants',
                        '1 recombinant', '2 recombinants']
        recomb_range = range(0, 3)
        if max_recombinant_count > 2:
            extra_header = [str(x) + " recombinants" for x in
                            range(3, max_recombinant_count + 1)]
            table_header = table_header + extra_header
            recomb_range = range(0, max_recombinant_count + 1)

        t = PrettyTable(table_header)
        t.padding_width = 1
        t.align = "l"

        # make all recombinant table
        for counter_name in ['all_counters', 'prod_counters']:
            counter_type = counter_name.split("_")[0]
            counter_set = eval(counter_name)
            for l in self.loci:
                counter = counter_set[l]
                count_array = [counter[x] for x in recomb_range]
                total_with_at_least_one = sum(count_array[1:])
                if total_with_at_least_one > 0:
                    percentages = [''] + [" (" + str(round(
                        (float(x) / total_with_at_least_one) * 100)) + "%)" for
                        x in
                        count_array[1:]]
                else:
                    percentages = [''] + [" (N/A%)" for x in count_array[1:]]
                row = []
                for i in recomb_range:
                    row.append(str(count_array[i]) + percentages[i])
                label = '{} {}'.format(counter_type, l)
                t.add_row([label] + row)

        outfile.write(t.get_string())
        outfile.write("\n")

        # If using unfiltered, name cells with more than two recombinants#
        if self.use_unfiltered:
            outfile.write(
                "\n\n#Cells with more than two recombinants for a locus#\n")
            found_multi = False
            for cell in cells.values():
                # if cell.count_total_recombinants('A') > 2 or
                # cell.count_total_recombinants('B') > 2:
                if cell.has_excess_recombinants(2):
                    outfile.write("###{}###\n".format(cell.name))
                    for l in self.loci:
                        count = cell.count_total_recombinants(
                            self.receptor_name, l)
                        outfile.write("{receptor}_{l}:\t{count}\n".format(
                            receptor=self.receptor_name, l=l, count=count))
                    outfile.write("\n")
                    found_multi = True
            if not found_multi:
                outfile.write("None\n\n")

        # Reporting iNKT cells
        # iNKT_count = len(NKT_cells)
        # if iNKT_count == 1:
        #    cell_word = 'cell'
        # else:
        #    cell_word = 'cells'
        # outfile.write("\n\n#iNKT cells#\nFound {iNKT_count} iNKT {cell_word}\n".format(iNKT_count=iNKT_count,
        #                                                                               cell_word=cell_word))
        # if iNKT_count > 0:
        #    for cell_name, ids in six.iteritems(NKT_cells):
        #        outfile.write("###{cell_name}###\n".format(cell_name=cell_name))
        #        outfile.write("TCRA:\t{}\nTCRB\t{}\n\n".format(ids[0], ids[1]))
        #

        # reporting invariant cells
        invariant_cells = []

        if self.invariant_cells is not None:
            for ivc in self.invariant_cells:
                ivc_loci = []
                found_ivcs = {}
                defining_locus = ivc.defining_locus
                if defining_locus in self.loci:
                    ivc_loci.append(defining_locus)
                    for cell in cells.values():
                        found_idents = {}
                        found_defining_locus, defining_id = ivc.check_for_match(
                            cell, defining_locus)
                        if found_defining_locus:
                            found_idents[ivc.defining_locus] = defining_id

                            for l in ivc.invariant_recombinants.keys():
                                if not l == defining_locus:
                                    ivc_loci.append(l)
                                    if l in cell.recombinants[
                                            ivc.receptor_type] and \
                                            cell.recombinants[
                                            ivc.receptor_type][
                                            l] is not None:
                                        found_other_locus, invar_id = ivc.check_for_match(
                                            cell, l)
                                        if found_other_locus:
                                            pass
                                        else:
                                            invar_id = "Invariant recombinant not found for {}_{}. {} found in total ({})".format(
                                                ivc.receptor_type, l, len(
                                                    cell.recombinants[
                                                        ivc.receptor_type][l]),
                                                cell.getMainRecombinantIdentifiersForLocus(
                                                    ivc.receptor_type, l))

                                    else:
                                        invar_id = "No sequences reconstructed for {}_{}".format(
                                            ivc.receptor_type, l)
                                    found_idents[l] = invar_id

                            found_ivcs[cell.name] = found_idents
                            invariant_cells.append(cell.name)

                if len(found_ivcs) > 0:
                    outfile.write("\n#{} cells#\n".format(ivc.name))

                    outfile.write("Expected: {}\n".format(ivc.expected_string))
                    outfile.write(
                        "Found {} possible cells.\n\n".format(len(found_ivcs)))

                    sorted_names = sorted(list(found_ivcs.keys()))
                    for n in sorted_names:
                        outfile.write("### {} ###\n".format(n))
                        ivc_details = found_ivcs[n]
                        for l in ivc_loci:
                            outfile.write(
                                "{}_{}: {}\n".format(ivc.receptor_type, l,
                                                     ivc_details[l]))
                    outfile.write("\n")

        # plot lengths of reconstructed sequences
        lengths = defaultdict(list)
        for cell in cells.values():
            for l in self.loci:
                lengths[l].extend(
                    cell.get_trinity_lengths(self.receptor_name, l))

        # plot length distributions
        quartiles = dict()
        for l in self.loci:
            q = self.get_quartiles(self.receptor_name, l)
            quartiles[l] = q

        for l in self.loci:
            q = quartiles[l]
            lns = lengths[l]
            try:
                if len(lns) > 1:
                    plt.figure()
                    plt.axvline(q[0], linestyle="--", color='k')
                    plt.axvline(q[1], linestyle="--", color='k')
                    sns.distplot(lns)
                    sns.despine()
                    plt.xlabel(
                        "{receptor}_{locus} reconstructed length (bp)".format(
                            receptor=self.receptor_name,
                            locus=l))
                    plt.ylabel("Density")
                    plt.savefig("{}_{}.pdf".format(length_filename_root, l))
                if len(lns) > 0:
                    with open("{}_{}.txt".format(length_filename_root, l),
                              'w') as f:
                        for l in sorted(lns):
                            f.write("{}\n".format(l))
            except:
                print(self.receptor_name)
                print(type(lns))
                print(lns)

        for cell_name in empty_cells:
            del cells[cell_name]

        if not self.keep_invariant:
            for cell_name in invariant_cells:
                del cells[cell_name]

        recombinant_data = []
        # Write out recombinant details for each cell
        with open("{}/recombinants.txt".format(outdir), 'w') as f:
            f.write(
                "cell_name\tlocus\trecombinant_id\tproductive\treconstructed_length\tCDR3aa\tCDR3nt\n")
            sorted_cell_names = sorted(list(cells.keys()))
            for cell_name in sorted_cell_names:
                cell = cells[cell_name]
                cell_data = {"cell_name": cell_name}
                for locus in self.loci:
                    cell_data.update({locus + "_unproductive": None,
                                      locus + "_productive": None})
                    recombinants = cell.recombinants[self.receptor_name][locus]
                    if recombinants is not None:
                        for r in recombinants:

                            # check for cdr3nt attribute (to make backwards compatible)
                            if hasattr(r, 'cdr3nt'):
                                cdr3nt = r.cdr3nt
                                cdr3 = r.cdr3

                                # lowercase CDR3 sequences if non-productive
                                if not r.productive:
                                    cdr3 = cdr3.lower()
                                    cdr3nt = cdr3nt.lower()

                            else:
                                cdr3nt = 'N/A'
                                cdr3 = 'N/A'

                            f.write(
                                "{name}\t{locus}\t{ident}\t{productive}\t{length}\t{cdr3}\t{cdr3nt}\n".format(
                                    name=cell_name, locus=locus,
                                    ident=r.identifier,
                                    productive=r.productive,
                                    length=len(r.trinity_seq),
                                    cdr3=cdr3,
                                    cdr3nt=cdr3nt))
                            if r.productive:
                                cell_data[locus + "_productive"] = r.identifier
                            else:
                                cell_data[
                                    locus + "_unproductive"] = r.identifier
                f.write("\n")
                recombinant_data.append(cell_data)
            f.write("\n\n")
            for cell_name in empty_cells:
                f.write(
                    "{cell_name}\tNo seqs found for {receptor}_{loci}\n".format(
                        cell_name=cell_name,
                        receptor=self.receptor_name,
                        loci=self.loci))
        #pdb.set_trace()
        recombinant_data = pd.DataFrame(recombinant_data)

        # make clonotype networks
        network_colours = io.read_colour_file(
            os.path.join(self.species_dir, "colours.csv"))
        component_groups = tracer_func.draw_network_from_cells(cells, outdir,
                                                               self.graph_format,
                                                               dot, neato,
                                                               self.draw_graphs,
                                                               self.receptor_name,
                                                               self.loci,
                                                               network_colours)

        # Print component groups to the summary#
        outfile.write(
            "\n#Clonotype groups#\n"
            "This is a text representation of the groups shown in clonotype_network_with_identifiers.pdf.\n"
            "It does not exclude cells that only share beta and not alpha.\n\n")
        for g in component_groups:
            outfile.write(", ".join(g))
            outfile.write("\n\n")

        # Build group membership dictionary
        group_membership = []
        for index, group in enumerate(component_groups):
            group_len = len(group)
            for cell in group:
                group_membership.append({"cell_name": cell,
                                         "clonal_group": index,
                                         "group_size": group_len})

        group_membership = pd.DataFrame(group_membership)
        
        if not group_membership.empty:
            cell_data = pd.merge(recombinant_data, group_membership, how='outer',
                                 on='cell_name')
        else:
            cell_data = recombinant_data.copy()
            cell_data['clonal_group'] = ''
            cell_data['group_size'] = ''
            
        cell_data.set_index("cell_name", inplace=True)
        cell_data.to_csv(os.path.join(outdir, "cell_data.csv"))
        #pdb.set_trace()
        # plot clonotype sizes
        plt.figure()
        clonotype_sizes = tracer_func.get_component_groups_sizes(cells,
                                                                 self.receptor_name,
                                                                 self.loci)
        w = 0.85
        x_range = range(1, len(clonotype_sizes) + 1)
        plt.bar(x_range, height=clonotype_sizes, width=w, color='black',
                align='center')
        plt.gca().set_xticks(x_range)
        plt.xlabel("Clonotype size")
        plt.ylabel("Clonotype count")
        plt.savefig("{}/clonotype_sizes.pdf".format(outdir))

        # write clonotype sizes to text file
        with open("{}/clonotype_sizes.txt".format(outdir), 'w') as f:
            data = zip(x_range, clonotype_sizes)
            f.write("clonotype_size\tclonotype_count\n")
            for t in data:
                f.write("{}\t{}\n".format(t[0], t[1]))

        outfile.close()

    def get_quartiles(self, receptor, locus):

        fasta = os.path.join(self.species_dir, 'combinatorial_recombinomes',
                             '{receptor}_{locus}.fa'.format(
                                 receptor=receptor, locus=locus))

        # need to remove the start N padding and end C sequence from the
        # lengths
        constant_fasta = os.path.join(self.species_dir, 'raw_seqs',
                                      '{receptor}_{locus}_C.fa'.format(
                                          receptor=receptor, locus=locus))
        C_len = len(next(SeqIO.parse(constant_fasta, "fasta")))

        seq = str(next(SeqIO.parse(fasta, "fasta")).seq)
        if seq.startswith("N"):
            r = re.compile(r"^N+")
            N_len = r.search(seq).end()
        else:
            N_len = 0

        adj = C_len + N_len

        lengths = array(
            [len(rec) - adj for rec in SeqIO.parse(fasta, "fasta")])
        quartiles = (percentile(lengths, 25), percentile(lengths, 75))
        return (quartiles)


class Tester(TracerTask):
    def __init__(self, **kwargs):
        if not kwargs:
            parser = argparse.ArgumentParser(
                description="Test TraCeR installation with small dataset",
                parents=[self.base_parser])
            parser.add_argument('--graph_format', '-f',
                                metavar="<GRAPH_FORMAT>",
                                help='graphviz output format [pdf]',
                                default='pdf')
            parser.add_argument('--no_networks',
                                help='skip attempts to draw network graphs',
                                action="store_true")
            parser.add_argument('--resume_with_existing_files', '-r',
                                help='look for existing intermediate files and use those instead of starting from scratch',
                                action="store_true")
            parser.add_argument('--output', '-o',
                                help='directory for output data of test')
            args = parser.parse_args(sys.argv[2:])

            self.resource_dir = args.resource_dir
            self.output_dir = args.output
            self.ncores = args.ncores
            self.config_file = args.config_file
            self.graph_format = args.graph_format
            self.no_networks = args.no_networks
            self.resume = args.resume_with_existing_files
        else:
            self.resource_dir = kwargs.get('resource_dir')
            self.output_dir = kwargs.get('output')
            self.ncores = kwargs.get('ncores')
            self.config_file = kwargs.get('config_file')
            self.graph_format = kwargs.get('graph_format', 'pdf')
            self.no_networks = kwargs.get('no_networks')
            self.resume = kwargs.get('resume_with_existing_files')

        self.config = self.read_config(self.config_file)

    def run(self):

        # test_dir = self.resolve_relative_path("test_data")
        test_dir = os.path.join(base_dir, 'test_data')
        if not os.path.exists(test_dir):
            tracer_dir = self.get_tracer_path()
            if not tracer_dir is None:
                test_dir = os.path.join(tracer_dir, 'test_data')
        test_names = ['cell1']
        if self.output_dir:
            out_dir = os.path.join(self.output_dir, 'results')
        else:
            out_dir = os.path.join(test_dir, 'results')
        for name in test_names:
            f1 = "{}/{}_1.fastq".format(test_dir, name)
            f2 = "{}/{}_2.fastq".format(test_dir, name)

            Assembler(resource_dir=self.resource_dir, ncores=str(self.ncores),
                      config_file=self.config_file,
                      resume_with_existing_files=self.resume,
                      species='Mmus', seq_method='imgt', fastq1=f1, fastq2=f2,
                      cell_name=name, output_dir=out_dir,
                      single_end=False, fragment_length=False,
                      fragment_sd=False, receptor_name='TCR',
                      loci=['A', 'B'], max_junc_len=50).run()

        Summariser(resource_dir=self.resource_dir, config_file=self.config_file, use_unfiltered=False,
                   keep_invariant=False,
                   graph_format=self.graph_format, no_networks=self.no_networks,
                   root_dir=out_dir, receptor_name='TCR',
                   loci=['A', 'B'], species='Mmus').run()


class Builder(TracerTask):
    """ Build Combinatorial Recombinomes for a given species """

    def __init__(self, **kwargs):
        self.leader_padding = 20
        if not kwargs:
            parser = argparse.ArgumentParser(
                description="Build resources from sequences",
                parents=[self.base_parser])
            parser.add_argument('--force_overwrite', '-f',
                                help='force overwrite of existing resources',
                                action='store_true')
            parser.add_argument('species', metavar="<SPECIES>",
                                help='species (eg Mmus)')
            parser.add_argument('receptor_name', metavar="<RECEPTOR_NAME>",
                                help='name of receptor (eg TCR)')
            parser.add_argument('locus_name', metavar="<LOCUS_NAME>",
                                help='name of locus (eg A)')
            parser.add_argument('N_padding', metavar="<N_PADDING>",
                                help='number of ambiguous N nucleotides between V and J',
                                type=int)
            parser.add_argument('colour', metavar="<COLOUR>", default='random',
                                help='colour for productive recombinants. Specify as HTML (eg E41A1C)\
                                or use "random"', type=self.check_colour)
            parser.add_argument('V_seqs', metavar="<V_SEQS>",
                                help='fasta file containing V gene sequences')
            parser.add_argument('J_seqs', metavar="<J_SEQS>",
                                help='fasta file containing J gene sequences')
            parser.add_argument('C_seq', metavar="<C_SEQ>",
                                help='fasta file containing single constant region sequence')
            parser.add_argument('D_seqs', metavar="<D_SEQS>", nargs='?',
                                default=False,
                                help='fasta file containing D gene sequences (optional)')
            parser.add_argument('--output_dir', '-o', metavar="<OPTIONAL_OUTPUT_DIR>", default=None,
                                help='optional output directory for built resources')

            args = parser.parse_args(sys.argv[2:])

            self.ncores = args.ncores
            self.force_overwrite = args.force_overwrite
            self.species = args.species
            self.receptor_name = args.receptor_name
            self.locus_name = args.locus_name
            self.N_padding = args.N_padding
            self.raw_seq_files = {}
            self.raw_seq_files['V'] = args.V_seqs
            self.raw_seq_files['J'] = args.J_seqs
            self.raw_seq_files['C'] = args.C_seq
            self.prod_colour = args.colour
            if args.D_seqs:
                self.raw_seq_files['D'] = args.D_seqs
            config_file = args.config_file
            self.output = args.output_dir

        else:
            self.ncores = kwargs.get('ncores')
            self.force_overwrite = kwargs.get('force_overwrite')
            self.species = kwargs.get('species')
            self.receptor_name = kwargs.get('receptor_name')
            self.locus_name = kwargs.get('locus_name')
            self.N_padding = kwargs.get('N_padding')
            self.raw_seq_files = {}
            self.raw_seq_files['V'] = kwargs.get('V_seqs')
            self.raw_seq_files['J'] = kwargs.get('J_seqs')
            self.raw_seq_files['C'] = kwargs.get('C_seq')
            self.prod_colour = kwargs.get('colour')
            if kwargs.get('D_seqs'):
                self.raw_seq_files['D'] = kwargs.get('D_seqs')

            config_file = kwargs.get('config_file')
            self.output = kwargs.get('output_dir')

        self.config = self.read_config(config_file)

        self.species_dir = self.get_species_root(self.species, root=self.output, build_mode=True)

    def run(self):

        # Check that there will not be git conflicts with inbuilt species
        # assert self.species not in ('Mmus', 'Hsap'), \
        #    "Cannot overwrite inbuilt species. Please choose a unique name " \
        #    "e.g. 'Mmus_1'"

        self.init_dirs()

        self.calculate_colours(self.prod_colour)
        VDJC_files = self.copy_raw_files()
        recombinome_fasta = self.make_recombinomes(VDJC_files)
        self.make_bowtie2_index(recombinome_fasta)
        missing_dbs = self.make_igblast_db(VDJC_files)
        for s in missing_dbs:
            print(
                "\nIMPORTANT: there is no IgBLAST database for {receptor}_{segment}\n"
                "Run build with {segment} segments for {receptor} before using tracer assemble\n"
                .format(receptor=self.receptor_name, segment=s))

    def check_colour(self, c):
        if c == 'random':
            return (c)
        else:
            try:
                if not c.startswith("#"):
                    c = "#" + c
                hex2color(c)
                return (c)
            except ValueError:
                msg = "{c} is not a valid html colour. Specify as xxXXxx".format(
                    c=c)
                raise argparse.ArgumentTypeError(msg)

    def calculate_colours(self, c):
        c = c.lower()
        pal = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33',
               '#a65628', '#f781bf']
        allowed_pal = copy.copy(pal)
        colour_file = os.path.join(self.species_dir, "colours.csv")
        if os.path.exists(colour_file):
            colour_map, used_colours = io.read_colour_file(colour_file,
                                                           return_used_list=True,
                                                           receptor_name=self.receptor_name)
            if len(used_colours) < len(pal):
                for uc in used_colours:
                    if uc in pal:
                        allowed_pal.remove(uc)
        else:
            used_colours = None
            colour_map = {}
        if c == 'random':
            prod_colour = random.choice(allowed_pal)

        else:
            if used_colours is not None:
                if (self.receptor_name in colour_map and
                    self.locus_name in colour_map[self.receptor_name]
                    and not colour_map[self.receptor_name][self.locus_name][
                        0] == c):
                    if c in used_colours and c not in allowed_pal:
                        msg = "{c} already in use. Please specify a different colour.".format(
                            c=c)
                        raise argparse.ArgumentTypeError(msg)

            prod_colour = c

        prod_rgb = hex2color(prod_colour)
        h, s, v = colorsys.rgb_to_hsv(*prod_rgb)

        nonprod_rgb = colorsys.hsv_to_rgb(h, s * 0.5, self.adj_v(v))
        nonprod_colour = str(rgb2hex(nonprod_rgb))

        d1 = {self.locus_name: (prod_colour, nonprod_colour)}

        if self.receptor_name in colour_map:
            colour_map[self.receptor_name].update(d1)
        else:
            colour_map[self.receptor_name] = d1

        io.write_colour_file(colour_file, colour_map)

    def adj_v(self, v):
        new_v = v * 1.3
        if new_v > 1:
            new_v = 1
        return (new_v)

    def check_duplicate(self, new_path, segment=None, descriptor="Resource"):
        error_string = "{descriptor} already exists for {receptor}_{locus}" \
            .format(descriptor=descriptor, receptor=self.receptor_name,
                    locus=self.locus_name)
        if segment:
            error_string += "_" + segment
        error_string += ". Use --force_overwrite to replace existing file."
        if os.path.isfile(new_path):
            assert self.force_overwrite, error_string

    def init_dirs(self):

        # Set up output directories
        subdirs = ['igblast_dbs', 'combinatorial_recombinomes', 'raw_seqs']

        io.makeOutputDir(self.species_dir)
        for d in subdirs:
            io.makeOutputDir(os.path.join(self.species_dir, d))

    def copy_raw_files(self):
        """ Move user-specified files to internal resources file system """

        gene_segs = 'VJC'

        VDJC_files = {}

        if 'D' in self.raw_seq_files:
            gene_segs += 'D'

        for s in gene_segs:
            fn = "{receptor}_{locus}_{s}.fa".format(receptor=self.receptor_name,
                                                    locus=self.locus_name, s=s)
            out_file = os.path.join(self.species_dir, 'raw_seqs', fn)
            VDJC_files[s] = out_file
            self.check_duplicate(out_file, segment=s,
                                 descriptor="Sequence File")
            shutil.copy(self.raw_seq_files[s], out_file)

        return VDJC_files

    def load_segment_seqs(self, filename):
        seqs = {}
        with open(filename, 'rU') as fn:
            for record in SeqIO.parse(fn, 'fasta'):
                seqs[record.id] = str(record.seq)

        return seqs

    def make_recombinomes(self, VDJC_files):

        out_fasta = os.path.join(
            self.species_dir, 'combinatorial_recombinomes',
            '{receptor}_{locus}.fa'.format(receptor=self.receptor_name,
                                           locus=self.locus_name))

        self.check_duplicate(out_fasta, descriptor="Combinatorial recombinome")

        seqs = {}

        for s in 'VJC':
            in_file = VDJC_files[s]
            seqs[s] = self.load_segment_seqs(in_file)

        # Logical check for C region
        if len(seqs['C']) > 1:
            print(
                "\nMore than one constant region sequence included in {C_file}."
                .format(C_file=self.raw_seq_files['C']))
            print("Please only provide one constant sequence.\n")
            sys.exit(1)

        const_seq = list(seqs['C'].values())[0].upper()
        N_junction_string = "N" * self.N_padding
        N_leader_string = "N" * self.leader_padding

        seqs_to_write = []

        # Compile sequences to write
        for V_name, V_seq in six.iteritems(seqs['V']):
            for J_name, J_seq in six.iteritems(seqs['J']):
                chr_name = ">chr={V_name}_{J_name}".format(J_name=J_name,
                                                           V_name=V_name)
                seq = N_leader_string + V_seq.lower() + N_junction_string + \
                    J_seq.lower() + const_seq
                seqs_to_write.append("{chr_name}\n{seq}\n"
                                     .format(seq=seq, chr_name=chr_name))

        with open(out_fasta, 'w') as f:
            for seq in seqs_to_write:
                f.write(seq)

        return out_fasta

    def make_bowtie2_index(self, recombinome_fasta):

        bowtie2_build = self.get_binary('bowtie2-build')
        index_base = os.path.join(
            self.species_dir, 'combinatorial_recombinomes',
            '{receptor}_{locus}'.format(receptor=self.receptor_name,
                                        locus=self.locus_name))

        self.check_duplicate(index_base + ".1.bt2", descriptor="Bowtie2 index")

        command = [bowtie2_build, '-q', recombinome_fasta, index_base]
        try:
            subprocess.check_call(command)
        except subprocess.CalledProcessError:
            print("bowtie2-build failed")

    def make_igblast_db(self, VDJC_files):

        igblast_dir = os.path.join(self.species_dir, 'igblast_dbs')

        makeblastdb = self.get_binary('makeblastdb')
        missing_dbs = []
        for s in 'VDJ':
            fn = "{receptor}_{segment}.fa".format(receptor=self.receptor_name,
                                                  segment=s)
            fasta_file = os.path.join(igblast_dir, fn)

            # Create file if it doesn't already exist
            open(fasta_file, 'a').close()

            # pdb.set_trace()
            if s in VDJC_files:
                with open(fasta_file) as e:
                    existing_seqs = SeqIO.to_dict(SeqIO.parse(e, "fasta"))
                with open(VDJC_files[s]) as n:
                    new_seqs = SeqIO.to_dict(SeqIO.parse(n, "fasta"))

                non_overwritten_seqs = []

                for seq_name, seq in six.iteritems(new_seqs):
                    if seq_name in existing_seqs:
                        if not self.force_overwrite:
                            non_overwritten_seqs.append(seq_name)
                        else:
                            existing_seqs.update({seq_name: seq})
                    else:
                        existing_seqs.update({seq_name: seq})
                with open(fasta_file, 'w') as f:
                    SeqIO.write(existing_seqs.values(), f, "fasta")

                if len(existing_seqs) == 0:
                    missing_dbs.append(s)

                if len(non_overwritten_seqs) > 0:
                    print('The follwing IgBLAST DB sequences for '
                          '{receptor}_{segment} already found in {file}.'
                          .format(receptor=self.receptor_name, segment=s,
                                  file=fasta_file))
                    print('These sequences were not overwritten. '
                          'Use --force_overwrite to replace with new ones')
                    for seq in non_overwritten_seqs:
                        print(seq)

                command = [makeblastdb, '-parse_seqids', '-dbtype', 'nucl',
                           '-in', fasta_file]
                try:
                    subprocess.check_call(command)
                except subprocess.CalledProcessError:
                    print("makeblastdb failed for {receptor}_{segment}"
                          .format(receptor=self.receptor_name, segment=s))

            else:
                with open(fasta_file) as e:
                    existing_seqs = SeqIO.to_dict(SeqIO.parse(e, "fasta"))
                if len(existing_seqs) == 0:
                    missing_dbs.append(s)

        return missing_dbs
