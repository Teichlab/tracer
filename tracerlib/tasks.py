from __future__ import print_function

import csv

import six
import matplotlib as mpl

from tracerlib import io, core
from tracerlib.io import check_binary

mpl.use('pdf')
import re
import seaborn as sns
from matplotlib import pyplot as plt
from tracerlib import base_dir
from tracerlib import tracer_func
from configparser import ConfigParser, NoOptionError
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

import pdb


class TracerTask(object):

    base_parser = argparse.ArgumentParser(add_help=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    base_parser.add_argument('--ncores', '-p', metavar="<CORES>", help='number of processor cores to use', type=int,
                             default=1)
    base_parser.add_argument('--config_file', '-c', metavar="<CONFIG_FILE>", help='config file to use',
                             default='~/.tracerrc')

    config = None

    def run(self):
        pass

    def get_binary(self, name):
        tool_key = name.lower() + '_path'
        user_path = None
        if self.config.has_option('tool_locations', tool_key):
            user_path = self.resolve_relative_path(self.config.get('tool_locations', tool_key))
        return check_binary(name, user_path)

    def read_config(self, config_file):
        # Read config file
        if not config_file:
            config_file = '~/.tracerrc'
        config_file = os.path.expanduser(config_file)
        if not os.path.isfile(config_file):
            print("Config file not found at ~/.tracerrc. Using default tracer.conf in repo...")
            config_file = os.path.join(base_dir, 'tracer.conf')
        tracer_func.check_config_file(config_file)
        config = ConfigParser()
        config.read(config_file)
        return config

    def resolve_relative_path(self, path):
        if not path.startswith("/"):
            base_directory = os.path.abspath(os.path.dirname(__file__))
            full_path = os.path.normpath("/{}/../{}".format(base_directory, path))
        else:
            full_path = path
        return full_path

    def print_cell_summary(self, cell, output_file):
        out_file = open(output_file, 'w')
        out_file.write('------------------\n{name}\n------------------\n'.format(name=cell.name))
        out_file.write('TCRA recombinants: {}\n'.format(cell.summarise_productivity('A')))
        out_file.write('TCRB recombinants: {}\n'.format(cell.summarise_productivity('B')))
        out_file.write('\n\n')
        out_file.write('#TCRA#\n')
        if cell.A_recombinants is None:
            out_file.write("No TCRA recombinants found\n\n")
        else:
            for rec in cell.A_recombinants:
                out_file.write(rec.get_summary())
                out_file.write("\n\n")
        out_file.write('#TCRB#\n')

        if cell.B_recombinants is None:
            out_file.write("No TCRB recombinants found\n\n")
        else:
            for rec in cell.B_recombinants:
                out_file.write(rec.get_summary())
                out_file.write("\n\n")
        out_file.close()

    def die_with_empty_cell(self, cell_name, output_dir, species):
        print("##No TCR recombinants found##")
        cell = core.Cell(cell_name, None, None, None, None, is_empty=True, species=species)
        self.print_cell_summary(cell,
                                "{output_dir}/unfiltered_TCR_seqs/unfiltered_TCRs.txt".format(output_dir=output_dir))
        with open("{output_dir}/unfiltered_TCR_seqs/{cell_name}.pkl".format(output_dir=output_dir,
                                                                            cell_name=cell.name), 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)
        cell.filter_recombinants()
        self.print_cell_summary(cell, "{output_dir}/filtered_TCR_seqs/filtered_TCRs.txt".format(output_dir=output_dir))
        with open("{output_dir}/filtered_TCR_seqs/{cell_name}.pkl".format(output_dir=output_dir,
                                                                          cell_name=cell.name), 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)
        exit(0)
    
    def get_resources_root(self, species):
        resources_dir = os.path.join(base_dir, 'resources')
        resources_root = os.path.join(resources_dir, species)
        return(resources_root)


class Assembler(TracerTask):

    def __init__(self, **kwargs):
        if not kwargs:
            
            # get list of all available species in resources
            
            parser = argparse.ArgumentParser(
                description="Reconstruct TCR sequences from RNAseq reads for a single cell", 
                parents=[self.base_parser], formatter_class=argparse.ArgumentDefaultsHelpFormatter)
            parser.add_argument('--resume_with_existing_files', '-r',
                                help='look for existing intermediate files and use those instead of starting from scratch',
                                action="store_true")
            parser.add_argument('--species', '-s',
                                help='species from which T cells were isolated - important to determination of iNKT cells',
                                choices=self.get_available_species(), default='Mmus')
            parser.add_argument('--receptor_name',
                                help="Name of receptor to reconstruct", default='TCR')
            parser.add_argument('--loci',
                                help="Space-separated list of loci to reconstruct for receptor", 
                                default=['A','B'], nargs = '+')
            parser.add_argument('--seq_method', '-m',
                                help='Method for constructing sequence to assess productivity, \
                                quantify expression and for output reporting. See README for details.',
                                choices=['imgt', 'assembly'], default='imgt')
            parser.add_argument('--single_end', help='set this if your sequencing data are single-end reads',
                                action="store_true")
            parser.add_argument('--fragment_length',
                                help='Estimated average fragment length in the sequencing library.'
                                     ' Used for Kallisto quantification. REQUIRED for single-end data.',
                                default=False)
            parser.add_argument('--fragment_sd',
                                help='Estimated standard deviation of average fragment length in the sequencing library.'
                                     ' Used for Kallisto quantification. REQUIRED for single-end data.',
                                default=False)
            parser.add_argument('--invariant_sequences',
                                help="Custom invariant sequence file. "
                                     "Use the example in 'resources/Mmus/invariant_seqs.csv'")

            parser.add_argument('fastq1', metavar="<FASTQ1>", help='first fastq file')
            parser.add_argument('fastq2', metavar="<FASTQ2>", help='second fastq file', nargs='?')
            parser.add_argument('cell_name', metavar="<CELL_NAME>", help='name of cell for file labels')
            parser.add_argument('output_dir', metavar="<OUTPUT_DIR>",
                                help='directory for output as <output_dir>/<cell_name>')

            args = parser.parse_args(sys.argv[2:])

            self.cell_name = args.cell_name
            self.fastq1 = args.fastq1
            self.single_end = args.single_end
            self.fastq2 = args.fastq2
            self.ncores = str(args.ncores)
            self.species = args.species
            self.seq_method = args.seq_method
            self.resume_with_existing_files = args.resume_with_existing_files
            self.fragment_length = args.fragment_length
            self.fragment_sd = args.fragment_sd
            self.output_dir = args.output_dir
            self.receptor_name = args.receptor_name
            self.loci = args.loci
            invariant_sequences = args.invariant_sequences
            config_file = args.config_file

        else:
            self.cell_name = kwargs.get('cell_name')
            self.fastq1 = kwargs.get('fastq1')
            self.fastq2 = kwargs.get('fastq2')
            self.ncores = kwargs.get('ncores')
            self.species = kwargs.get('species')
            self.seq_method = kwargs.get('seq_method')
            self.resume_with_existing_files = kwargs.get('resume_with_existing_files')
            self.output_dir = kwargs.get('output_dir')
            self.single_end = kwargs.get('single_end')
            self.fragment_length = kwargs.get('fragment_length')
            self.fragment_sd = kwargs.get('fragment_sd')
            self.receptor_name = kwargs.get('receptor_name')
            self.loci = kwargs.get('loci')
            invariant_sequences = kwargs.get('invariant_sequences')
            config_file = kwargs.get('config_file')

        self.config = self.read_config(config_file)
        self.locus_names = ["TCRA", "TCRB"]

        # Check the fastq config is correct
        if not self.single_end:
            assert self.fastq2, "Only one fastq file specified. Either set --single_end or provide second fastq."
        else:
            self.fastq2 = None
            if self.fastq2:
                print("Two fastq files given with --single-end option. Ignoring second file.")
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

        # Get Invariant Sequences
        if not invariant_sequences:
            invariant_sequences = self.resolve_relative_path(os.path.join('resources', self.species,
                                                                          'invariant_seqs.csv'))
        if not os.path.isfile(invariant_sequences):
            raise OSError('2', 'Invariant Sequence file not found', invariant_sequences)

        self.invariant_sequences = io.parse_invariant_seqs(invariant_sequences)
        
    def get_available_species(self):
        resources_dir = os.path.join(base_dir, 'resources')
        species_dirs = next(os.walk(resources_dir))[1]
        return(species_dirs)
        
    def run(self, **kwargs):

        # Set-up output directories
        root_output_dir = os.path.abspath(self.output_dir)
        io.makeOutputDir(root_output_dir)
        self.output_dir = root_output_dir + "/" + self.cell_name

        io.makeOutputDir(self.output_dir)

        data_dirs = ['aligned_reads', 'Trinity_output', 'IgBLAST_output', 
                     'unfiltered_{receptor}_seqs'.format(receptor = self.receptor_name),'expression_quantification', 
                     'filtered_{receptor}_seqs'.format(receptor = self.receptor_name)]
        for d in data_dirs:
            io.makeOutputDir("{}/{}".format(self.output_dir, d))

        # Perform TraCeR's core functions
        self.align()
        self.de_novo_assemble()
        cell = self.ig_blast()
        self.quantify(cell)

        self.print_cell_summary(
            cell, "{output_dir}/unfiltered_{receptor}_seqs/unfiltered_{receptor}s.txt".format(output_dir=self.output_dir))

        # Save cell in a pickle
        with open("{output_dir}/unfiltered_{receptor}_seqs/{cell_name}.pkl".format(output_dir=self.output_dir,
                                                                            cell_name=cell.name, receptor=self.receptor_name), 
                                                                            'wb') as pf:
            pickle.dump(cell, pf, protocol=0)
        print("##Filtering by read count##")
        cell.filter_recombinants()
        fasta_filename = "{output_dir}/filtered_{receptor}seqs/{cell_name}_{receptor}seqs.fa".format(output_dir=self.output_dir,
                                                                                        cell_name=self.cell_name,
                                                                                        receptor=self.receptor_name)
        fasta_file = open(fasta_filename, 'w')
        fasta_file.write(cell.get_fasta_string())
        fasta_file.close()
        self.print_cell_summary(cell, "{output_dir}/filtered_{receptor}_seqs/filtered_{receptor}s.txt".format(
                                                                                            output_dir=self.output_dir),
                                                                                            receptor=self.receptor_name)
        with open("{output_dir}/filtered_{receptor}_seqs/{cell_name}.pkl".format(output_dir=self.output_dir,
                                                                          cell_name=cell.name,
                                                                          receptor=self.receptor_name), 'wb') as pf:
            pickle.dump(cell, pf, protocol=0)

    def get_index_location(self, name):
        location = os.path.join(base_dir, 'resources', self.species, name)

        return location

    def align(self):
        bowtie2 = self.get_binary('bowtie2')

        synthetic_genome_path = self.get_index_location('combinatorial_recombinomes')
        # Align with bowtie
        tracer_func.bowtie2_alignment(
            bowtie2, self.ncores, self.receptor_name, self.loci, self.output_dir, self.cell_name, 
            synthetic_genome_path, self.fastq1, self.fastq2, self.resume_with_existing_files, self.single_end)
        print()

    def de_novo_assemble(self):

        trinity = self.get_binary('trinity')

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
            trinity_grid_conf = self.resolve_relative_path(self.config.get('trinity_options', 'trinity_grid_conf'))
        else:
            trinity_grid_conf = False

        # De novo assembly with trinity
        trinity_JM = self.config.get('trinity_options', 'max_jellyfish_memory')
        trinity_version = self.config.get('trinity_options', 'trinity_version')
        successful_files = tracer_func.assemble_with_trinity(
            trinity, self.receptor_name, self.loci, self.output_dir, self.cell_name, self.ncores, trinity_grid_conf, 
            trinity_JM, trinity_version, self.resume_with_existing_files, self.single_end, self.species)
        if len(successful_files) == 0:
            print("No successful Trinity assemblies")
            self.die_with_empty_cell(self.cell_name, self.output_dir, self.species)

        print()

    def ig_blast(self):
        igblastn = self.get_binary('igblastn')

        # Reference data locations
        igblast_index_location = self.get_index_location('igblast_dbs')
        imgt_seq_location = self.get_index_location('raw_seqs')

        igblast_seqtype = self.config.get('IgBlast_options', 'igblast_seqtype')


        # IgBlast of assembled contigs
        tracer_func.run_IgBlast(igblastn, self.receptor_name, self.loci, self.output_dir, self.cell_name, igblast_index_location,
                                igblast_seqtype, self.species, self.resume_with_existing_files)
        print()
        
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cell = io.parse_IgBLAST(self.receptor_name, self.loci, self.output_dir, self.cell_name, imgt_seq_location, 
                                    self.species, self.seq_method, self.invariant_sequences)
            if cell.is_empty:
                self.die_with_empty_cell(self.cell_name, self.output_dir, self.species)

        return cell

    def quantify(self, cell):
        kallisto = self.get_binary('kallisto')
        
        if not self.config.has_option('kallisto_transcriptomes', self.species):
            raise OSError("No transcriptome reference specified for {species}. Please specify location in config file."
                          .format(species = self.species))
        else:
            kallisto_base_transcriptome = self.resolve_relative_path(self.config.get('kallisto_transcriptomes',
                                                                                 self.species))

        # Quantification with kallisto
        tracer_func.quantify_with_kallisto(
            kallisto, cell, self.output_dir, self.cell_name, kallisto_base_transcriptome, self.fastq1, self.fastq2,
            self.ncores, self.resume_with_existing_files, self.single_end, self.fragment_length, self.fragment_sd)
        print()

        counts = tracer_func.load_kallisto_counts("{}/expression_quantification/abundance.tsv".format(self.output_dir))

        for locus, recombinants in six.iteritems(cell.all_recombinants):
            if recombinants is not None:
                for rec in recombinants:
                    tpm = counts[locus][rec.contig_name]
                    rec.TPM = tpm


class Summariser(TracerTask):

    def __init__(self, **kwargs):

        if not kwargs:
            parser = argparse.ArgumentParser(description="Summarise set of cells with reconstructed TCR sequences",
                                             parents=[self.base_parser])
            parser.add_argument('--use_unfiltered', '-u', help='use unfiltered recombinants', action="store_true")
            parser.add_argument('--keep_inkt', '-i', help='ignore iNKT cells when constructing networks',
                                action="store_true")
            parser.add_argument('--graph_format', '-f', metavar="<GRAPH_FORMAT>", help='graphviz output format [pdf]',
                                default='pdf')
            parser.add_argument('--no_networks', help='skip attempts to draw network graphs', action = "store_true")
            parser.add_argument('dir', metavar="<DIR>",
                                help='directory containing subdirectories for each cell to be summarised')
            args = parser.parse_args(sys.argv[2:])

            self.root_dir = os.path.abspath(args.dir)
            self.graph_format = args.graph_format
            self.keep_inkt = args.keep_inkt
            self.use_unfiltered = args.use_unfiltered
            self.draw_graphs = not args.no_networks
            config_file = args.config_file
        else:
            self.use_unfiltered = kwargs.get('use_unfiltered')
            self.root_dir = os.path.abspath(kwargs.get('root_dir'))
            self.draw_graphs = not (kwargs.get('no_networks'))
            self.graph_format = kwargs.get('graph_format')
            self.keep_inkt = kwargs.get('keep_inkt')
            config_file = kwargs.get('config_file')

        # Read config file
        self.config = self.read_config(config_file)

    def run(self):

        if self.draw_graphs:
            dot = self.resolve_relative_path(self.config.get('tool_locations', 'dot_path'))
            neato = self.resolve_relative_path(self.config.get('tool_locations', 'neato_path'))

            # check that executables from config file can be used
            not_executable = []
            for name, x in six.iteritems({"dot": dot, "neato": neato}):
                if not io.is_exe(x):
                    not_executable.append((name, x))
            if len(not_executable) > 0:
                print()
                print("Could not execute the following required tools. Check your configuration file.")
                for t in not_executable:
                    print( t[0], t[1])
                print()
                exit(1)
        else:
            dot = ""
            neato = ""

        cells = {}
        empty_cells = []
        NKT_cells = {}
        subdirectories = next(os.walk(self.root_dir))[1]

        if self.use_unfiltered:
            pkl_dir = "unfiltered_TCR_seqs"
            outdir = "{}/unfiltered_TCR_summary".format(self.root_dir)
            # outfile = open("{root_dir}/unfiltered_TCR_summary.txt".format(root_dir=root_dir), 'w')
            # length_filename_root = "{}/unfiltered_reconstructed_lengths_TCR".format(root_dir)

        else:
            pkl_dir = "filtered_TCR_seqs"
            outdir = "{}/filtered_TCR_summary".format(self.root_dir)
            # outfile = open("{root_dir}/filtered_TCR_summary.txt".format(root_dir=root_dir), 'w')
            # length_filename_root = "{}/filtered_reconstructed_lengths_TCR".format(root_dir)

        io.makeOutputDir(outdir)

        outfile = open("{}/TCR_summary.txt".format(outdir), 'w')
        length_filename_root = "{}/reconstructed_lengths_TCR".format(outdir)

        for d in subdirectories:
            cell_pkl = "{root_dir}/{d}/{pkl_dir}/{d}.pkl".format(pkl_dir=pkl_dir, d=d, root_dir=self.root_dir)
            if os.path.isfile(cell_pkl):
                with open(cell_pkl, 'rb') as pkl:
                    cl = pickle.load(pkl)
                cells[d] = cl
                if cl.is_empty:
                    empty_cells.append(d)
                if cl.is_inkt:
                    NKT_cells[d] = (cl.is_inkt, cl.getMainRecombinantIdentifiersForLocus('B'))
        count_of_cells_with_alpha_recovered = 0
        count_of_cells_with_beta_recovered = 0
        count_of_cells_with_paired_recovered = 0
        for cell_name, cell in six.iteritems(cells):
            prod_a_count = cell.count_productive_recombinants('A')
            prod_b_count = cell.count_productive_recombinants('B')
            if prod_a_count > 0:
                count_of_cells_with_alpha_recovered += 1
            if prod_b_count > 0:
                count_of_cells_with_beta_recovered += 1
            if prod_a_count > 0 and prod_b_count > 0:
                count_of_cells_with_paired_recovered += 1

        total_cells = len(cells)

        outfile.write(
            "TCRA reconstruction:\t{count_of_cells_with_alpha_recovered} / {total_cells} ({alpha_percent}%)\n"
            "TCRB reconstruction:\t{count_of_cells_with_beta_recovered} / {total_cells} ({beta_percent}%)\n"
            "Paired productive chains\t{count_of_cells_with_paired_recovered} / {total_cells} ({paired_percent}%)\n\n"
            .format(
                paired_percent=round((count_of_cells_with_paired_recovered / float(total_cells)) * 100, 1),
                total_cells=total_cells,
                alpha_percent=round((count_of_cells_with_alpha_recovered / float(total_cells)) * 100, 1),
                beta_percent=round((count_of_cells_with_beta_recovered / float(total_cells)) * 100, 1),
                count_of_cells_with_beta_recovered=count_of_cells_with_beta_recovered,
                count_of_cells_with_paired_recovered=count_of_cells_with_paired_recovered,
                count_of_cells_with_alpha_recovered=count_of_cells_with_alpha_recovered))

        all_alpha_counter = Counter()
        all_beta_counter = Counter()
        prod_alpha_counter = Counter()
        prod_beta_count = Counter()

        counters = {'all_alpha': Counter(), 'all_beta': Counter(), 'prod_alpha': Counter(), 'prod_beta': Counter()}

        for cell in cells.values():
            counters['all_alpha'].update({cell.count_total_recombinants('A'): 1})
            counters['all_beta'].update({cell.count_total_recombinants('B'): 1})
            counters['prod_alpha'].update({cell.count_productive_recombinants('A'): 1})
            counters['prod_beta'].update({cell.count_productive_recombinants('B'): 1})

        max_recombinant_count = max(list(counters['all_alpha'].keys()) + list(counters['all_beta'].keys()))
        table_header = ['', '0 recombinants', '1 recombinant', '2 recombinants']
        recomb_range = range(0, 3)
        if max_recombinant_count > 2:
            extra_header = [str(x) + " recombinants" for x in range(3, max_recombinant_count + 1)]
            table_header = table_header + extra_header
            recomb_range = range(0, max_recombinant_count + 1)

        t = PrettyTable(table_header)
        t.padding_width = 1
        t.align = "l"
        for label in ['all_alpha', 'all_beta', 'prod_alpha', 'prod_beta']:
            counter = counters[label]
            count_array = [counter[x] for x in recomb_range]
            total_with_at_least_one = sum(count_array[1:])
            if total_with_at_least_one > 0:
                percentages = [''] + [" (" + str(round((float(x) / total_with_at_least_one) * 100)) + "%)" for x in
                                      count_array[1:]]
            else:
                percentages = [''] + [" (N/A%)" for x in count_array[1:]]
            row = []
            for i in recomb_range:
                row.append(str(count_array[i]) + percentages[i])

            t.add_row([label] + row)
        outfile.write(t.get_string())

        # If using unfiltered, name cells with more than two recombinants#
        if self.use_unfiltered:
            outfile.write("\n\n#Cells with more than two recombinants for a locus#\n")
            found_multi = False
            for cell in cells.values():
                if cell.count_total_recombinants('A') > 2 or cell.count_total_recombinants('B') > 2:
                    outfile.write("###{}###\n".format(cell.name))
                    outfile.write("TCRA:\t{}\nTCRB:\t{}\n\n".format(cell.count_total_recombinants('A'),
                                                                    cell.count_total_recombinants('B')))
                    found_multi = True
            if not found_multi:
                outfile.write("None\n\n")

        # Reporting iNKT cells
        iNKT_count = len(NKT_cells)
        if iNKT_count == 1:
            cell_word = 'cell'
        else:
            cell_word = 'cells'
        outfile.write("\n\n#iNKT cells#\nFound {iNKT_count} iNKT {cell_word}\n".format(iNKT_count=iNKT_count,
                                                                                       cell_word=cell_word))
        if iNKT_count > 0:
            for cell_name, ids in six.iteritems(NKT_cells):
                outfile.write("###{cell_name}###\n".format(cell_name=cell_name))
                outfile.write("TCRA:\t{}\nTCRB\t{}\n\n".format(ids[0], ids[1]))

        # plot lengths of reconstructed sequences
        lengths = {'A': [], 'B': []}
        for cell in cells.values():
            for locus in lengths.keys():
                lengths[locus] = lengths[locus] + cell.get_trinity_lengths(locus)

        # plot TCRA length distributions
        if len(lengths['A']) > 1:
            plt.figure()
            plt.axvline(334, linestyle="--", color='k')
            plt.axvline(344, linestyle="--", color='k')
            sns.distplot(lengths['A'])
            sns.despine()
            plt.xlabel("TCRa reconstructed length (bp)")
            plt.ylabel("Density")
            plt.savefig("{}A.pdf".format(length_filename_root))
        if len(lengths['A']) > 0:
            with open("{}A.txt".format(length_filename_root), 'w') as f:
                for l in sorted(lengths['A']):
                    f.write("{}\n".format(l))

        # plot TCRB length distributions
        if len(lengths['B']) > 1:
            plt.figure()
            plt.axvline(339, linestyle="--", color='k')
            plt.axvline(345, linestyle="--", color='k')
            sns.distplot(lengths['B'])
            sns.despine()
            plt.xlabel("TCRb reconstructed length (bp)")
            plt.ylabel("Density")
            plt.savefig("{}B.pdf".format(length_filename_root))
        if len(lengths['B']) > 0:
            with open("{}B.txt".format(length_filename_root), 'w') as f:
                for l in sorted(lengths['B']):
                    f.write("{}\n".format(l))

        for cell_name in empty_cells:
            del cells[cell_name]

        if not self.keep_inkt:
            for cell_name in NKT_cells.keys():
                del cells[cell_name]

        # make clonotype networks
        component_groups = tracer_func.draw_network_from_cells(cells, outdir, self.graph_format, dot, neato,
                                                               self.draw_graphs)

        # Print component groups to the summary#
        outfile.write(
            "\n###Clonotype groups###\n"
            "This is a text representation of the groups shown in clonotype_network_with_identifiers.pdf."
            " It does not exclude cells that only share beta and not alpha.\n\n")
        for g in component_groups:
            outfile.write(", ".join(g))
            outfile.write("\n\n")

        # plot clonotype sizes
        plt.figure()
        clonotype_sizes = tracer_func.get_component_groups_sizes(cells)
        w = 0.85
        x_range = range(1, len(clonotype_sizes) + 1)
        plt.bar(x_range, height=clonotype_sizes, width=w, color='black', align='center')
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

        # Write out recombinant details for each cell
        with open("{}/recombinants.txt".format(outdir), 'w') as f:
            f.write("cell_name\tlocus\trecombinant_id\tproductive\treconstructed_length\n")
            sorted_cell_names = sorted(list(cells.keys()))
            for cell_name in sorted_cell_names:
                cell = cells[cell_name]
                for locus in "AB":
                    recombinants = cell.all_recombinants[locus]
                    if recombinants is not None:
                        for r in recombinants:
                            f.write(
                                "{name}\t{locus}\t{ident}\t{productive}\t{length}\n".format(
                                    name=cell_name, locus=locus, ident=r.identifier,
                                    productive=r.productive, length=len(r.trinity_seq)))
                f.write("\n")
            f.write("\n\n")
            for cell_name in empty_cells:
                f.write("{}\tNo TCRs found\n".format(cell_name))

        outfile.close()


class Tester(TracerTask):

    def __init__(self, **kwargs):
        if not kwargs:
            parser = argparse.ArgumentParser(description="Test TraCeR installation with small dataset",
                                             parents=[self.base_parser])
            parser.add_argument('--graph_format', '-f', metavar="<GRAPH_FORMAT>", help='graphviz output format [pdf]',
                                default='pdf')
            parser.add_argument('--no_networks', help='skip attempts to draw network graphs', action="store_true")
            args = parser.parse_args(sys.argv[2:])

            self.ncores = args.ncores
            self.config_file = args.config_file
            self.graph_format = args.graph_format
            self.no_networks = args.no_networks
        else:
            self.ncores = kwargs.get('ncores')
            self.config_file = kwargs.get('config_file')
            self.graph_format = kwargs.get('graph_format', 'pdf')
            self.no_networks = kwargs.get('no_networks')

    def run(self):

        # test_dir = self.resolve_relative_path("test_data")
        test_dir = os.path.join(base_dir, 'test_data')
        test_names = ['cell1']
        out_dir = "{}/results".format(test_dir)
        for name in test_names:
            f1 = "{}/{}_1.fastq".format(test_dir, name)
            f2 = "{}/{}_2.fastq".format(test_dir, name)

            Assembler(ncores=str(self.ncores), config_file=self.config_file, resume_with_existing_files=True,
                      species='Mmus', seq_method='imgt', fastq1=f1, fastq2=f2, cell_name=name, output_dir=out_dir,
                      single_end=False, fragment_length=False, fragment_sd=False, receptor_name='TCR',
                      loci=['A', 'B']).run()

        Summariser(config_file=self.config_file, use_unfiltered=False, keep_inkt=False,
                   graph_format=self.graph_format, no_networks=self.no_networks, root_dir=out_dir) \
            .run()


class Builder(TracerTask):

    """ Build Combinatorial Recombinomes for a given species """

    def __init__(self, **kwargs):
        self.leader_padding = 20
        if not kwargs:
            parser = argparse.ArgumentParser(description="Build resources from sequences", parents=[self.base_parser])
            parser.add_argument('--force_overwrite', '-f', help='force overwrite of existing resources',
                                action='store_true')
            parser.add_argument('species', metavar="<SPECIES>", help='species (eg Mmus)')
            parser.add_argument('receptor_name', metavar="<RECEPTOR_NAME>", help='name of receptor (eg TCR)')
            parser.add_argument('locus_name', metavar="<LOCUS_NAME>", help='name of locus (eg A)')
            parser.add_argument('N_padding', metavar="<N_PADDING>", 
                                 help='number of ambiguous N nucleotides between V and J', type=int)
            parser.add_argument('V_seqs', metavar="<V_SEQS>", help='fasta file containing V gene sequences')
            parser.add_argument('J_seqs', metavar="<J_SEQS>", help='fasta file containing J gene sequences')
            parser.add_argument('C_seq', metavar="<C_SEQ>", 
                                help='fasta file containing single constant region sequence')
            parser.add_argument('D_seqs', metavar="<D_SEQS>", nargs='?', default=False,
                                help='fasta file containing D gene sequences (optional)')
            
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
            if args.D_seqs:
                self.raw_seq_files['D'] = args.D_seqs
            config_file = args.config_file
            
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
            if kwargs.get('D_seqs'):
                self.raw_seq_files['D'] = kwargs.get('D_seqs')
            
            config_file = kwargs.get('config_file')

        self.config = self.read_config(config_file)
        self.species_dir = self.get_resources_root(self.species)

    def run(self):

        # Check that there will not be git conflicts with inbuilt species
        assert self.species not in ('Mmus', 'Hsap'), \
            "Cannot overwrite inbuilt species. Please choose a unique name " \
            "e.g. 'Mmus_1'"

        self.init_dirs()
        VDJC_files = self.copy_raw_files()
        recombinome_fasta = self.make_recombinomes(VDJC_files)
        self.make_bowtie2_index(recombinome_fasta)
        missing_dbs = self.make_igblast_db(VDJC_files)
        for s in missing_dbs:
            print("\nIMPORTANT: there is no IgBLAST database for {receptor}_{segment}\n"\
                  "Run build with {segment} segments for {receptor} before using tracer assemble\n"
                  .format(receptor = self.receptor_name, segment=s))

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
            print("\nMore than one constant region sequence included in {C_file}." \
                  .format(self.raw_seq_files['C']))
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

            #pdb.set_trace()
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
