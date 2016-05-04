from __future__ import print_function
import six
import matplotlib as mpl

from tracer.tracerlib import io, core

mpl.use('pdf')
import re
import seaborn as sns
from matplotlib import pyplot as plt
from tracer import base_dir
from tracer.tracerlib import tracer_func
from configparser import ConfigParser
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


class Launcher(object):

    def __init__(self):
        parser = argparse.ArgumentParser(
            description='TraCeR: reconstruction of TCR sequences from single-cell RNAseq data',
            usage=''' tracer <mode> [<args>]

                  Modes are :

                  - assemble: assemble TCR sequences from single-cell RNA-sequencing reads
                  - summarise: summarise TCR sequences from set of cells, build clonotype networks
                  - test : use a small dataset from three cells to test TraCeR installation

                  use tracer <mode> -h for specific help
                  ''')
        parser.add_argument('mode', metavar="<MODE>", help='tracer mode (assemble, summarise or get_test_data)',
                            choices=['assemble', 'summarise', 'summarize', 'test'])
        args = parser.parse_args(sys.argv[1:2])

        if not hasattr(self, args.mode):
            print('Unrecognised mode')
            parser.print_help()
            exit(1)

        getattr(self, args.mode)()

    ##ASSEMBLE
    def assemble(self, **kwargs):
        if not kwargs:
            parser = argparse.ArgumentParser(
                description="Reconstruct TCR sequences from RNAseq reads for a single cell")
            parser.add_argument('--ncores', '-p', metavar="<CORES>", help='number of processor cores to use', type=int,
                                default=1)
            parser.add_argument('--config_file', '-c', metavar="<CONFIG_FILE>", help='config file to use [tracer.conf]',
                                default='tracer.conf')
            parser.add_argument('--resume_with_existing_files', '-r',
                                help='look for existing intermediate files and use those instead of starting from scratch',
                                action="store_true")
            parser.add_argument('--species', '-s',
                                help='species from which T cells were isolated - important to determination of iNKT cells',
                                choices=['Mmus', 'Hsap'], default='Mmus')
            parser.add_argument('--seq_method', '-m',
                                help='Method for constructing sequence to assess productivity, quantify expression and for output reporting. See README for details.',
                                choices=['imgt', 'assembly'], default='imgt')
            parser.add_argument('--single_end', help='set this if your sequencing data are single-end reads',
                                action="store_true")
            parser.add_argument('--fragment_length',
                                help='Estimated average fragment length in the sequencing library. Used for Kallisto quantification. REQUIRED for single-end data.',
                                default=False)
            parser.add_argument('--fragment_sd',
                                help='Estimated standard deviation of average fragment length in the sequencing library. Used for Kallisto quantification. REQUIRED for single-end data.',
                                default=False)
            parser.add_argument('fastq1', metavar="<FASTQ1>", help='first fastq file')
            parser.add_argument('fastq2', metavar="<FASTQ2>", help='second fastq file', nargs='?')
            parser.add_argument('cell_name', metavar="<CELL_NAME>", help='name of cell for file labels')
            parser.add_argument('output_dir', metavar="<OUTPUT_DIR>",
                                help='directory for output as <output_dir>/<cell_name>')

            args = parser.parse_args(sys.argv[2:])

            cell_name = args.cell_name
            fastq1 = args.fastq1
            single_end = args.single_end
            if not single_end:
                if not args.fastq2:
                    print("Only one fastq file specified. Either set --single_end or provide second fastq.")
                    exit(1)
                else:
                    fastq2 = args.fastq2
            else:
                fastq2 = None
                if args.fastq2:
                    print("Two fastq files given with --single-end option. Ignoring second file.")
                if not args.fragment_length and not args.fragment_sd:
                    print('Must specify estimated average fragment length (--fragment_length) and standard deviation (--fragment_sd) for use with single-end data')
                    exit(1)
                elif not args.fragment_length:
                    print('Must specify estimated average fragment length (--fragment_length) for use with single-end data')
                    exit(1)
                elif not args.fragment_sd:
                    print('Must specify estimated fragment length standard deviation (--fragment_sd) for use with single-end data')
                    exit(1)

            ncores = str(args.ncores)
            config_file = args.config_file
            species = args.species
            seq_method = args.seq_method
            resume_with_existing_files = args.resume_with_existing_files
            fragment_length = args.fragment_length
            fragment_sd = args.fragment_sd
            output_dir = args.output_dir

        else:
            cell_name = kwargs['cell_name']
            fastq1 = kwargs['fastq1']
            fastq2 = kwargs['fastq2']
            ncores = kwargs['ncores']
            config_file = kwargs['config_file']
            species = kwargs['species']
            seq_method = kwargs['seq_method']
            resume_with_existing_files = kwargs['resume_with_existing_files']
            output_dir = kwargs['output_dir']
            single_end = kwargs['single_end']
            fragment_length = kwargs['fragment_length']
            fragment_sd = kwargs['fragment_sd']

        # Read config file
        tracer_func.check_config_file(config_file)
        config = ConfigParser()
        config.read(config_file)

        bowtie2 = self.resolve_relative_path(config.get('tool_locations', 'bowtie2_path'))
        igblast = self.resolve_relative_path(config.get('tool_locations', 'igblast_path'))
        kallisto = self.resolve_relative_path(config.get('tool_locations', 'kallisto_path'))
        trinity = self.resolve_relative_path(config.get('tool_locations', 'trinity_path'))

        if config.has_option('trinity_options', 'trinity_grid_conf'):
            trinity_grid_conf = self.resolve_relative_path(config.get('trinity_options', 'trinity_grid_conf'))
        else:
            trinity_grid_conf = False

        # Trinity version
        if not config.has_option('trinity_options', 'trinity_version'):
            try:
                subprocess.check_output([trinity, '--version'])
            except subprocess.CalledProcessError as err:
                if re.search('v2', err.output.decode('utf-8')):
                    config.set('trinity_options', 'trinity_version', '2')
                else:
                    config.set('trinity_options', 'trinity_version', '1')

        synthetic_genome_path = self.resolve_relative_path(config.get('bowtie2_options', 'synthetic_genome_index_path'))
        igblast_index_location = self.resolve_relative_path(config.get('IgBlast_options', 'igblast_index_location'))
        igblast_seqtype = config.get('IgBlast_options', 'igblast_seqtype')
        imgt_seq_location = self.resolve_relative_path(config.get('IgBlast_options', 'imgt_seq_location'))

        kallisto_base_transcriptome = self.resolve_relative_path(config.get('kallisto_options', 'base_transcriptome'))

        # check that executables from config file can be used
        not_executable = []
        for name, x in six.iteritems({"bowtie2": bowtie2, "igblast": igblast, "kallisto": kallisto, "trinity": trinity}):
            if not io.is_exe(x):
                not_executable.append((name, x))
        if len(not_executable) > 0:
            print()
            print("Could not execute the following required tools. Check your configuration file.")
            for t in not_executable:
                print(t[0], t[1])
            print()
            exit(1)

        # set-up output directories
        root_output_dir = os.path.abspath(output_dir)
        io.makeOutputDir(root_output_dir)
        output_dir = root_output_dir + "/" + cell_name

        io.makeOutputDir(output_dir)

        data_dirs = ['aligned_reads', 'Trinity_output', 'IgBLAST_output', 'unfiltered_TCR_seqs',
                     'expression_quantification', 'filtered_TCR_seqs']
        for d in data_dirs:
            io.makeOutputDir("{}/{}".format(output_dir, d))

        locus_names = ["TCRA", "TCRB"]

        should_resume = resume_with_existing_files

        self.bowtie2_alignment(bowtie2, ncores, locus_names, output_dir, cell_name, synthetic_genome_path, fastq1,
                               fastq2, should_resume, single_end)
        print()
        trinity_JM = config.get('trinity_options', 'max_jellyfish_memory')
        trinity_version = config.get('trinity_options', 'trinity_version')
        self.assemble_with_trinity(trinity, locus_names, output_dir, cell_name, ncores, trinity_grid_conf, trinity_JM,
                                   trinity_version, should_resume, single_end, species)
        print()
        self.run_IgBlast(igblast, locus_names, output_dir, cell_name, igblast_index_location, igblast_seqtype, species,
                         should_resume)
        print()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cell = io.parse_IgBLAST(locus_names, output_dir, cell_name, imgt_seq_location, species, seq_method)
            if cell.is_empty:
                self.die_with_empty_cell(cell_name, output_dir, species)

        self.quantify_with_kallisto(kallisto, cell, output_dir, cell_name, kallisto_base_transcriptome, fastq1, fastq2,
                                    ncores, should_resume, single_end, fragment_length, fragment_sd)

        print()

        counts = tracer_func.load_kallisto_counts("{}/expression_quantification/abundance.tsv".format(output_dir))

        # pdb.set_trace():
        for locus, recombinants in six.iteritems(cell.all_recombinants):
            if recombinants is not None:
                for rec in recombinants:
                    tpm = counts[locus][rec.contig_name]
                    rec.TPM = tpm

        self.print_cell_summary(cell,
                                "{output_dir}/unfiltered_TCR_seqs/unfiltered_TCRs.txt".format(output_dir=output_dir))
        with open("{output_dir}/unfiltered_TCR_seqs/{cell_name}.pkl".format(output_dir=output_dir,
                                                                            cell_name=cell.name), 'wb') as pf:
            pickle.dump(cell, pf)
        print("##Filtering by read count##")
        cell.filter_recombinants()
        fasta_filename = "{output_dir}/filtered_TCR_seqs/{cell_name}_TCRseqs.fa".format(output_dir=output_dir,
                                                                                        cell_name=cell_name)
        fasta_file = open(fasta_filename, 'w')
        fasta_file.write(cell.get_fasta_string())
        fasta_file.close()
        self.print_cell_summary(cell, "{output_dir}/filtered_TCR_seqs/filtered_TCRs.txt".format(output_dir=output_dir))
        with open("{output_dir}/filtered_TCR_seqs/{cell_name}.pkl".format(output_dir=output_dir,
                                                                          cell_name=cell.name), 'wb') as pf:
            pickle.dump(cell, pf)

    def resolve_relative_path(self, path):
        if not path.startswith("/"):
            base_directory = os.path.abspath(os.path.dirname(__file__))
            full_path = os.path.normpath("{}/{}".format(base_directory, path))
        else:
            full_path = path
        return full_path

    def bowtie2_alignment(self, bowtie2, ncores, locus_names, output_dir, cell_name, synthetic_genome_path, fastq1,
                          fastq2, should_resume, single_end):
        print("##Finding TCR-derived reads##")

        if should_resume:
            for locus in locus_names:
                aligned_read_path = "{}/aligned_reads/{}_{}_".format(output_dir, cell_name, locus)
                fastq1_out = "{}1.fastq".format(aligned_read_path)
                fastq2_out = "{}2.fastq".format(aligned_read_path)
                if os.path.isfile(fastq1_out) and os.path.isfile(fastq2_out):
                    print("Resuming with existing TCRA and B reads")
                    return

        for locus in locus_names:
            print("##{}##".format(locus))
            sam_file = "{}/aligned_reads/{}_{}.sam".format(output_dir, cell_name, locus)
            if not single_end:
                fastq_out_1 = open("{}/aligned_reads/{}_{}_1.fastq".format(output_dir, cell_name, locus), 'w')
                fastq_lines_1 = []
                fastq_out_2 = open("{}/aligned_reads/{}_{}_2.fastq".format(output_dir, cell_name, locus), 'w')
                fastq_lines_2 = []

                command = [bowtie2, '--no-unal', '-p', ncores, '-k', '1', '--np', '0', '--rdg', '1,1', '--rfg', '1,1',
                           '-x', "/".join([synthetic_genome_path, locus]), '-1', fastq1, '-2', fastq2, '-S', sam_file]

                subprocess.check_call(command)

                # now to split the sam file for Trinity.

                with open(sam_file) as sam_in:
                    for line in sam_in:
                        if not line.startswith("@"):

                            line = line.rstrip()
                            line = line.split("\t")
                            name = line[0]
                            seq = line[9]
                            qual = line[10]
                            flag = int(line[1])
                            mate_flag = "{0:b}".format(flag)[-7]
                            mate_mapped_flag = "{0:b}".format(flag)[-4]
                            revcomp_flag = "{0:b}".format(flag)[-5]

                            if revcomp_flag == "1":
                                seq = str(Seq(seq).reverse_complement())
                                qual = qual[::-1]
                            if mate_mapped_flag == "0":
                                if mate_flag == "1":
                                    name_ending = "/1"
                                    fastq_lines_1.append(
                                        "@{name}{name_ending}\n{seq}\n+\n{qual}\n".format(name=name, seq=seq,
                                                                                          name_ending=name_ending,
                                                                                          qual=qual))
                                else:
                                    name_ending = "/2"
                                    fastq_lines_2.append(
                                        "@{name}{name_ending}\n{seq}\n+\n{qual}\n".format(name=name, seq=seq,
                                                                                          name_ending=name_ending,
                                                                                          qual=qual))

                for line in fastq_lines_1:
                    fastq_out_1.write(line)
                for line in fastq_lines_2:
                    fastq_out_2.write(line)

                fastq_out_1.close()
                fastq_out_2.close()
            else:
                fastq_out = open("{}/aligned_reads/{}_{}.fastq".format(output_dir, cell_name, locus), 'w')
                command = [bowtie2, '--no-unal', '-p', ncores, '-k', '1', '--np', '0', '--rdg', '1,1', '--rfg', '1,1',
                           '-x', "/".join([synthetic_genome_path, locus]), '-U', fastq1, '-S', sam_file]

                subprocess.check_call(command)

                with open(sam_file) as sam_in:
                    for line in sam_in:
                        if not line.startswith("@"):

                            line = line.rstrip()
                            line = line.split("\t")
                            name = line[0]
                            seq = line[9]
                            qual = line[10]
                            flag = int(line[1])
                            if not flag == 0:
                                revcomp_flag = "{0:b}".format(flag)[-5]
                            else:
                                revcomp_flag = "0"

                            if revcomp_flag == "1":
                                seq = str(Seq(seq).reverse_complement())
                                qual = qual[::-1]
                            fastq_out.write("@{name}\n{seq}\n+\n{qual}\n".format(name=name, seq=seq, qual=qual))
                    fastq_out.close()

    def assemble_with_trinity(self, trinity, locus_names, output_dir, cell_name, ncores, trinity_grid_conf, JM,
                              version, should_resume, single_end, species):
        print("##Assembling Trinity Contigs##")

        if should_resume:
            trinity_report_successful = "{}/Trinity_output/successful_trinity_assemblies.txt".format(output_dir)
            trinity_report_unsuccessful = "{}/Trinity_output/unsuccessful_trinity_assemblies.txt".format(output_dir)
            if (os.path.isfile(trinity_report_successful) and os.path.isfile(trinity_report_unsuccessful)) and (
                            os.path.getsize(trinity_report_successful) > 0 or os.path.getsize(
                        trinity_report_unsuccessful) > 0):
                print("Resuming with existing Trinity output")
                return

        command = [trinity]
        if trinity_grid_conf:
            command = command + ['--grid_conf', trinity_grid_conf]

        memory_string = '--max_memory' if (version == '2') else '--JM'
        command = command + ['--seqType', 'fq', memory_string, JM, '--CPU', ncores, '--full_cleanup']

        for locus in locus_names:
            print("##{}##".format(locus))
            trinity_output = "{}/Trinity_output/{}_{}".format(output_dir, cell_name, locus)
            aligned_read_path = "{}/aligned_reads/{}_{}".format(output_dir, cell_name, locus)
            if not single_end:
                file1 = "{}_1.fastq".format(aligned_read_path)
                file2 = "{}_2.fastq".format(aligned_read_path)
                command = command + ["--left", file1, "--right", file2, "--output",
                                     '{}/Trinity_output/Trinity_{}_{}'.format(output_dir, cell_name, locus)]
            else:
                file = "{}.fastq".format(aligned_read_path)
                command = command + ["--single", file, "--output",
                                     '{}/Trinity_output/Trinity_{}_{}'.format(output_dir, cell_name, locus)]
            try:
                subprocess.check_call(command)
                shutil.move('{}/Trinity_output/Trinity_{}_{}.Trinity.fasta'.format(output_dir, cell_name, locus),
                            '{}/Trinity_output/{}_{}.Trinity.fasta'.format(output_dir, cell_name, locus))
            except (subprocess.CalledProcessError, IOError):
                print("Trinity failed for locus")

        # clean up unsuccessful assemblies
        sleep(10)  # this gives the cluster filesystem time to catch up and stops weird things happening
        successful_files = glob.glob("{}/Trinity_output/*.fasta".format(output_dir))
        unsuccessful_directories = next(os.walk("{}/Trinity_output".format(output_dir)))[1]
        for directory in unsuccessful_directories:
            shutil.rmtree("{}/Trinity_output/{}".format(output_dir, directory))
        successful_file_summary = "{}/Trinity_output/successful_trinity_assemblies.txt".format(output_dir)
        unsuccessful_file_summary = "{}/Trinity_output/unsuccessful_trinity_assemblies.txt".format(output_dir)

        successful_files = io.clean_file_list(successful_files)
        unsuccessful_directories = io.clean_file_list(unsuccessful_directories)

        success_out = open(successful_file_summary, "w")
        fail_out = open(unsuccessful_file_summary, "w")

        successful = defaultdict(list)
        unsuccessful = defaultdict(list)

        successful_ordered_files = set()
        unsuccessful_ordered_files = set()

        for filename in successful_files:
            # success_out.write("{}\n".format(filename))
            parsed_name = io.get_filename_and_locus(filename)
            successful[parsed_name[0]].append(parsed_name[1])
            successful_ordered_files.add(parsed_name[0])
        successful_ordered_files = sorted(list(successful_ordered_files))

        for filename in unsuccessful_directories:
            # fail_out.write("{}\n".format(filename))
            parsed_name = io.get_filename_and_locus(filename)
            unsuccessful[parsed_name[0]].append(parsed_name[1])
            unsuccessful_ordered_files.add(parsed_name[0])
        unsuccessful_ordered_files = sorted(list(unsuccessful_ordered_files))

        successful = io.sort_locus_names(successful)
        unsuccessful = io.sort_locus_names(unsuccessful)

        for file in successful_ordered_files:
            success_out.write("{}\t{}\n".format(file, successful[file]))

        for file in unsuccessful_ordered_files:
            fail_out.write("{}\t{}\n".format(file, unsuccessful[file]))

        success_out.close()
        fail_out.close()

        # remove pointless .readcount files
        readcount_files = glob.glob("{}/aligned_reads/*.readcount".format(output_dir))
        for f in readcount_files:
            os.remove(f)

        if len(unsuccessful_directories) == 2:
            print("No successful Trinity assemblies")
            self.die_with_empty_cell(cell_name, output_dir, species)

    def run_IgBlast(self, igblast, locus_names, output_dir, cell_name, index_location, ig_seqtype, species,
                    should_resume):
        print("##Running IgBLAST##")

        igblast_species_lookup = {'Hsap': 'human', 'Mmus': 'mouse'}
        igblast_species = igblast_species_lookup[species]

        if should_resume:
            igblast_out_A = "{output_dir}/IgBLAST_output/{cell_name}_TCRA.IgBLASTOut".format(output_dir=output_dir,
                                                                                             cell_name=cell_name)
            igblast_out_B = "{output_dir}/IgBLAST_output/{cell_name}_TCRB.IgBLASTOut".format(output_dir=output_dir,
                                                                                             cell_name=cell_name)
            if (os.path.isfile(igblast_out_A) and os.path.getsize(igblast_out_A) > 0) or (
                        os.path.isfile(igblast_out_B) and os.path.getsize(igblast_out_B) > 0):
                print("Resuming with existing IgBLAST output")
                return

        databases = {}
        for segment in ['v', 'd', 'j']:
            databases[segment] = "{}/imgt_tcr_db_{}.fa".format(index_location, segment)

        # Lines below suppress Igblast warning about not having an auxliary file.
        # Taken from http://stackoverflow.com/questions/11269575/how-to-hide-output-of-subprocess-in-python-2-7
        DEVNULL = open(os.devnull, 'wb')

        for locus in locus_names:
            print("##{}##".format(locus))
            trinity_fasta = "{}/Trinity_output/{}_{}.Trinity.fasta".format(output_dir, cell_name, locus)
            if os.path.isfile(trinity_fasta):
                command = [igblast, '-germline_db_V', databases['v'], '-germline_db_D', databases['d'],
                           '-germline_db_J', databases['j'], '-domain_system', 'imgt', '-organism', igblast_species,
                           '-ig_seqtype', ig_seqtype, '-show_translation', '-num_alignments_V', '5',
                           '-num_alignments_D', '5', '-num_alignments_J', '5', '-outfmt', '7', '-query', trinity_fasta]
                igblast_out = "{output_dir}/IgBLAST_output/{cell_name}_{locus}.IgBLASTOut".format(output_dir=output_dir,
                                                                                                  cell_name=cell_name,
                                                                                                  locus=locus)
                out = open(igblast_out, 'w')
                # print(" ").join(pipes.quote(s) for s in command)
                subprocess.check_call(command, stdout=out, stderr=DEVNULL)
                out.close()

        DEVNULL.close()

    def quantify_with_kallisto(self, kallisto, cell, output_dir, cell_name, kallisto_base_transcriptome, fastq1, fastq2,
                               ncores, should_resume, single_end, fragment_length, fragment_sd):
        print("##Running Kallisto##")
        if should_resume:
            if os.path.isfile("{}/expression_quantification/abundance.tsv".format(output_dir)):
                print("Resuming with existing Kallisto output")
                return

        print("##Making Kallisto indices##")
        kallisto_dirs = ['kallisto_index']
        for d in kallisto_dirs:
            io.makeOutputDir("{}/expression_quantification/{}".format(output_dir, d))
        fasta_filename = "{output_dir}/unfiltered_TCR_seqs/{cell_name}_TCRseqs.fa".format(output_dir=output_dir,
                                                                                          cell_name=cell_name)
        fasta_file = open(fasta_filename, 'w')
        fasta_file.write(cell.get_fasta_string())
        fasta_file.close()

        output_transcriptome = "{}/expression_quantification/kallisto_index/{}_transcriptome.fa".format(output_dir,
                                                                                                        cell_name)
        with open(output_transcriptome, 'w') as outfile:
            for fname in [kallisto_base_transcriptome, fasta_filename]:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)

        idx_file = "{}/expression_quantification/kallisto_index/{}_transcriptome.idx".format(output_dir, cell_name)

        index_command = [kallisto, 'index', '-i', idx_file, output_transcriptome]
        subprocess.check_call(index_command)
        print("##Quantifying with Kallisto##")

        if not single_end:
            if not fragment_length:
                kallisto_command = [kallisto, 'quant', '-i', idx_file, '-t', ncores, '-o',
                                    "{}/expression_quantification".format(output_dir), fastq1, fastq2]
            else:
                kallisto_command = [kallisto, 'quant', '-i', idx_file, '-t', ncores, '-l', fragment_length, '-o',
                                    "{}/expression_quantification".format(output_dir), fastq1, fastq2]
        else:
            kallisto_command = [kallisto, 'quant', '-i', idx_file, '-t', ncores, '--single', '-l', fragment_length,
                                '-s', fragment_sd, '-o', "{}/expression_quantification".format(output_dir), fastq1]
        subprocess.check_call(kallisto_command)

        # delete index file because it's huge and unecessary. Delete transcriptome file
        # os.remove(idx_file)
        # os.remove(output_transcriptome)
        shutil.rmtree("{}/expression_quantification/kallisto_index/".format(output_dir))

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
            pickle.dump(cell, pf)
        cell.filter_recombinants()
        self.print_cell_summary(cell, "{output_dir}/filtered_TCR_seqs/filtered_TCRs.txt".format(output_dir=output_dir))
        with open("{output_dir}/filtered_TCR_seqs/{cell_name}.pkl".format(output_dir=output_dir,
                                                                          cell_name=cell.name), 'wb') as pf:
            pickle.dump(cell, pf)
        exit(0)

    ##SUMMARISE

    def summarise(self, **kwargs):

        if not kwargs:
            parser = argparse.ArgumentParser(description="Summarise set of cells with reconstructed TCR sequences")
            parser.add_argument('--config_file', '-c', metavar="<CONFIG_FILE>", help='config file to use [tracer.conf]',
                                default='tracer.conf')
            parser.add_argument('--use_unfiltered', '-u', help='use unfiltered recombinants', action="store_true")
            parser.add_argument('--keep_inkt', '-i', help='ignore iNKT cells when constructing networks',
                                action="store_true")
            parser.add_argument('--graph_format', '-f', metavar="<GRAPH_FORMAT>", help='graphviz output format [pdf]',
                                default='pdf')
            parser.add_argument('dir', metavar="<DIR>",
                                help='directory containing subdirectories for each cell to be summarised')
            args = parser.parse_args(sys.argv[2:])

            root_dir = os.path.abspath(args.dir)
            graph_format = args.graph_format
            config_file = args.config_file
            keep_inkt = args.keep_inkt
            use_unfiltered = args.use_unfiltered
        else:
            config_file = kwargs['config_file']
            use_unfiltered = kwargs['use_unfiltered']
            keep_inkt = kwargs['keep_inkt']
            graph_format = kwargs['graph_format']
            root_dir = os.path.abspath(kwargs['root_dir'])

        # Read config file
        tracer_func.check_config_file(config_file)
        config = ConfigParser()
        config.read(config_file)

        dot = self.resolve_relative_path(config.get('tool_locations', 'dot_path'))
        neato = self.resolve_relative_path(config.get('tool_locations', 'neato_path'))

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

        cells = {}
        empty_cells = []
        NKT_cells = {}
        subdirectories = next(os.walk(root_dir))[1]

        if use_unfiltered:
            pkl_dir = "unfiltered_TCR_seqs"
            outdir = "{}/unfiltered_TCR_summary".format(root_dir)
            # outfile = open("{root_dir}/unfiltered_TCR_summary.txt".format(root_dir=root_dir), 'w')
            # length_filename_root = "{}/unfiltered_reconstructed_lengths_TCR".format(root_dir)

        else:
            pkl_dir = "filtered_TCR_seqs"
            outdir = "{}/filtered_TCR_summary".format(root_dir)
            # outfile = open("{root_dir}/filtered_TCR_summary.txt".format(root_dir=root_dir), 'w')
            # length_filename_root = "{}/filtered_reconstructed_lengths_TCR".format(root_dir)

        io.makeOutputDir(outdir)

        outfile = open("{}/TCR_summary.txt".format(outdir), 'w')
        length_filename_root = "{}/reconstructed_lengths_TCR".format(outdir)

        for d in subdirectories:
            cell_pkl = "{root_dir}/{d}/{pkl_dir}/{d}.pkl".format(pkl_dir=pkl_dir, d=d, root_dir=root_dir)
            if os.path.isfile(cell_pkl):
                cl = pickle.load(open(cell_pkl, 'rb'))
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
            "TCRA reconstruction:\t{count_of_cells_with_alpha_recovered} / {total_cells} ({alpha_percent}%)\nTCRB reconstruction:\t{count_of_cells_with_beta_recovered} / {total_cells} ({beta_percent}%)\nPaired productive chains\t{count_of_cells_with_paired_recovered} / {total_cells} ({paired_percent}%)\n\n".format(
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
        if use_unfiltered:
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

        if not keep_inkt:
            for cell_name in NKT_cells.keys():
                del cells[cell_name]

        # make clonotype networks
        component_groups = tracer_func.draw_network_from_cells(cells, outdir, graph_format, dot, neato)

        # Print component groups to the summary#
        outfile.write(
            "\n###Clonotype groups###\nThis is a text representation of the groups shown in clonotype_network_with_identifiers.pdf. It does not exclude cells that only share beta and not alpha.\n\n")
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

    summarize = summarise

    def test(self):
        parser = argparse.ArgumentParser(description="Test TraCeR installation with small dataset")
        parser.add_argument('--ncores', '-p', metavar="<CORES>", help='number of processor cores to use', type=int,
                            default=1)
        parser.add_argument('--config_file', '-c', metavar="<CONFIG_FILE>", help='config file to use [tracer.conf]',
                            default='tracer.conf')
        args = parser.parse_args(sys.argv[2:])

        # test_dir = self.resolve_relative_path("test_data")
        test_dir = os.path.join(base_dir, 'test_data')
        test_names = ['cell1']
        out_dir = "{}/results".format(test_dir)
        for name in test_names:
            f1 = "{}/{}_1.fastq".format(test_dir, name)
            f2 = "{}/{}_2.fastq".format(test_dir, name)

            self.assemble(ncores=str(args.ncores), config_file=args.config_file, resume_with_existing_files=False,
                          species='Mmus', seq_method='imgt', fastq1=f1, fastq2=f2, cell_name=name, output_dir=out_dir,
                          single_end=False, fragment_length=False, fragment_sd=False)

        self.summarise(config_file=args.config_file, use_unfiltered=False, keep_inkt=False, graph_format='pdf',
                       root_dir=out_dir)
