##############################################################################
#                         Functions for use with                             #
# TraCeR - a tool to reconstruct TCR sequences from single-cell RNA-seq data #    
#                                                                            #
# Please see README and LICENCE for details of use and licence conditions.   #
# This software was written by Mike Stubbington (ms31@sanger.ac.uk) from the #
# Teichmann Lab, EMBL-EBI and WTSI (www.teichlab.org). Latest versions are   #
# available for download at www.github.com/teichlab/tracer.                  #
#                                                                            #
#      Copyright (c) 2015, 2016 EMBL - European Bioinformatics Institute     #
#      Copyright (c) 2016 Genome Research Ltd.                               #
#      Author: M.J.T. Stubbington ms31@sanger.ac.uk                          #
##############################################################################

from __future__ import print_function

import glob
import os
import re
import shutil
import subprocess
from collections import defaultdict, Counter
import csv
from time import sleep

import Levenshtein
import networkx as nx
import six
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from tracerlib.core import Cell, Recombinant
from tracerlib import io


def process_chunk(chunk):
    store_VDJ_rearrangement_summary = False
    store_junction_details = False
    store_alignment_summary = False
    store_hit_table = False
    alignment_summary = []
    hit_table = []
    looking_for_end = False
    return_dict = defaultdict(list)
    for line_x in chunk:

        if store_VDJ_rearrangement_summary:
            VDJ_rearrangement_summary = line_x.split("\t")
            for i in VDJ_rearrangement_summary:
                return_dict['VDJ_rearrangement_summary'].append(i)
            store_VDJ_rearrangement_summary = False

        elif store_junction_details:
            junction_details = line_x.split("\t")
            for i in junction_details:
                return_dict["junction_details"].append(i)
            store_junction_details = False

        elif store_alignment_summary:
            if not line_x.startswith("#"):
                if line_x.startswith("Total"):
                    store_alignment_summary = False
                else:
                    return_dict['alignment_summary'].append(line_x)

        elif store_hit_table:
            if not looking_for_end:
                if not line_x.startswith("#"):
                    return_dict['hit_table'].append(line_x)
                    looking_for_end = True
            else:
                if line_x.startswith("#"):
                    store_hit_table = False
                else:
                    return_dict['hit_table'].append(line_x)

        elif line_x.startswith('# Query'):
            query_name = line_x.split(" ")[2]
            query_length = line_x.split(" ")[3]
            return_dict['query_length'] = int(query_length.split("=")[1])
            # return_dict['query_name'] = query_name

        elif line_x.startswith('# V-(D)-J rearrangement summary'):
            store_VDJ_rearrangement_summary = True

        elif line_x.startswith('# V-(D)-J junction details'):
            store_junction_details = True

        elif line_x.startswith('# Alignment summary'):
            store_alignment_summary = True

        elif line_x.startswith('# Hit table'):
            store_hit_table = True

    return (query_name, return_dict)


def find_possible_alignments(sample_dict, locus_names, cell_name, IMGT_seqs, output_dir, species, seq_method,
                             constant_seqs):
    alignment_dict = defaultdict(dict)
    recombinants = {'TCRA': [], 'TCRB': []}
    for locus in locus_names:
        data_for_locus = sample_dict[locus]
        if data_for_locus is not None:
            for query_name, query_data in six.iteritems(data_for_locus):
                processed_hit_table = process_hit_table(query_name, query_data, locus)

                if processed_hit_table is not None:
                    (returned_locus, good_hits, rearrangement_summary) = processed_hit_table
                    junction_list = query_data['junction_details']

                    best_V = remove_allele_stars(rearrangement_summary[0].split(",")[0])

                    junc_string = "".join(junction_list)
                    junc_string = remove_NA(junc_string)

                    if returned_locus in "BD":
                        best_J = remove_allele_stars(rearrangement_summary[2].split(",")[0])
                    elif returned_locus in "AG":
                        best_J = remove_allele_stars(rearrangement_summary[1].split(",")[0])

                    identifier = best_V + "_" + junc_string + "_" + best_J

                    ##line attempting to add alignment summary to data for use with PCR comparisons
                    alignment_summary = query_data['alignment_summary']

                    all_V_names = [remove_allele_stars(x) for x in rearrangement_summary[0].split(',')]

                    if locus == "TCRB":
                        all_J_names = [remove_allele_stars(x) for x in rearrangement_summary[2].split(',')]
                    elif locus == "TCRA":
                        all_J_names = [remove_allele_stars(x) for x in rearrangement_summary[1].split(',')]

                    all_poss_identifiers = set()
                    for V in all_V_names:
                        for J in all_J_names:
                            i = V + "_" + junc_string + "_" + J
                            all_poss_identifiers.add(i)

                    # get original sequence from Trinity file - needed for summary of reconstructed lengths.
                    # Only use the VDJ portion found by IgBLAST
                    trinity_file = "{output_dir}/Trinity_output/{cell_name}_{locus}.Trinity.fasta".format(
                        locus=locus, output_dir=output_dir, cell_name=cell_name)
                    with open(trinity_file, 'rU') as tf:
                        for record in SeqIO.parse(tf, 'fasta'):
                            if query_name in record.id:
                                trinity_seq = record

                    if 'reversed' in good_hits[0][1]:
                        trinity_seq = trinity_seq.reverse_complement().seq
                    else:
                        trinity_seq = trinity_seq.seq
                    start_coord, end_coord = get_coords(good_hits)
                    trinity_seq = str(trinity_seq[start_coord:end_coord])

                    (imgt_reconstructed_seq, is_productive, bestVJNames) = get_fasta_line_for_contig_imgt(
                        rearrangement_summary, junction_list, good_hits, returned_locus, IMGT_seqs, cell_name,
                        query_name, species, constant_seqs)
                    del (is_productive)
                    del (bestVJNames)

                    if seq_method == 'imgt':
                        (fasta_line_for_contig, is_productive, bestVJNames) = get_fasta_line_for_contig_imgt(
                            rearrangement_summary, junction_list, good_hits, returned_locus, IMGT_seqs, cell_name,
                            query_name, species, constant_seqs)

                    elif seq_method == 'assembly':
                        fasta_line_for_contig = trinity_seq
                        (is_productive, bestVJNames) = get_fasta_line_for_contig_assembly(trinity_seq, good_hits,
                                                                                          returned_locus, IMGT_seqs,
                                                                                          cell_name, query_name)

                    if len(junc_string) < 50:
                        rec = Recombinant(contig_name=query_name, locus=returned_locus, identifier=identifier,
                                          all_poss_identifiers=all_poss_identifiers, productive=is_productive[0],
                                          stop_codon=is_productive[1], in_frame=is_productive[2], TPM=0.0,
                                          dna_seq=fasta_line_for_contig, hit_table=good_hits,
                                          summary=rearrangement_summary, junction_details=junction_list,
                                          best_VJ_names=bestVJNames, alignment_summary=alignment_summary,
                                          trinity_seq=trinity_seq, imgt_reconstructed_seq=imgt_reconstructed_seq)
                        recombinants[locus].append(rec)

    if recombinants:
        for locus, rs in six.iteritems(recombinants):
            # Adding code to collapse sequences with very low Levenshtein distances caused by confusion between
            # TRAVxD and TRAVx segments with different alignment lengths from IgBlast.
            recombinants[locus] = collapse_close_sequences(rs, locus)

        # cell_name, A_recombinants, B_recombinants, G_recombinants, D_recombinants, is_empty=False, species="Mmus")
        cell = Cell(cell_name, recombinants['TCRA'], recombinants['TCRB'], None, None, species=species)
    else:
        cell = Cell(cell_name, None, None, None, None, species=species)

    # pdb.set_trace()
    return (cell)


def get_coords(hit_table):
    found_V = False
    found_J = False
    for entry in hit_table:
        if entry[0] == 'V':
            if not found_V:
                start = int(entry[8]) - 1
                found_V = True
        if entry[0] == 'J':
            if not found_J:
                end = int(entry[9])
                found_J = True
    return (start, end)


def remove_NA(junc_string):
    new_string = junc_string.replace("N/A", "")
    return (new_string)


def remove_allele_stars(segment):
    p = re.compile(r"(.+)\*\d+")
    m = p.search(segment)
    return (m.group(1))


def process_hit_table(query_name, query_data, locus):
    hit_table = query_data['hit_table']
    rearrangement_summary = query_data['VDJ_rearrangement_summary']

    e_value_cutoff = 5e-3

    found_V = set()
    found_D = set()
    found_J = set()

    good_hits = []

    segment_locus_pattern = re.compile(r"TRAV.+DV.+")

    for entry in hit_table:
        entry = entry.split("\t")
        segment = entry[2]
        if segment_locus_pattern.search(segment):
            segment_locus = "AD"
        else:
            segment_locus = segment[2]
        segment_type = segment[3]
        e_value = float(entry[12])

        if locus[3] in segment_locus:
            if e_value < e_value_cutoff:
                if segment_type == "V":
                    found_V.add(locus)
                    good_hits.append(entry)
                elif segment_type == "J":
                    found_J.add(locus)
                    good_hits.append(entry)
            else:
                if segment_type == "D":
                    percent_identity = float(entry[3])
                    if percent_identity == 100:
                        found_D.add(locus)
                        good_hits.append(entry)

    if locus == "TCRA":
        if "TCRA" in found_V and "TCRA" in found_J:
            return ("A", good_hits, rearrangement_summary)
        else:
            return (None)

    # elif locus == "D":
    #    if "D" in found_V and "D" in found_J:
    #        return("D", good_hits, rearrangement_summary)
    #    else:
    #        return(None)

    elif locus == "TCRB":
        if "TCRB" in found_V and "TCRB" in found_J:
            return ("B", good_hits, rearrangement_summary)
        else:
            return (None)

            # elif locus == "G":
            #    if "G" in found_V and "G" in found_J:
            #        return("G", good_hits, rearrangement_summary)
            #    else:
            #        return(None)


def get_fasta_line_for_contig_imgt(rearrangement_summary, junction_details, hit_table, locus, IMGT_seqs,
                                   sample_name, query_name, species, constant_seqs):

    # use first 258 bases of TRBC because they're the same between C1 and C2
    # for TRGC use first 150 bases. Found by aligning the 4 C region transcripts and taking consensus. Ignored start of TCRG-C4-201 because it's only in that one.
    # use first 360 nt of TRBC1 because they're very nearly the same between TRBC1 and TRBCC2

    found_best_V = False
    found_best_D = False
    found_best_J = False

    V_pattern = re.compile(r"TR[ABGD]V\d")
    D_pattern = re.compile(r"TR[BD]D\d")
    J_pattern = re.compile(r"TR[ABGD]J\d")

    for hit in hit_table:
        segment = hit[2]
        if V_pattern.search(segment) and not found_best_V:
            V_locus_key = "TR{}V".format(segment[2])
            best_V_name = segment
            # Remove forward slashes from shared A/D gene names to be the same as in the IMGT files.
            segment = segment.replace("/", "_")
            best_V_seq = IMGT_seqs[V_locus_key][segment]

            # hit[11] is the end of the V sequence
            best_V_seq = best_V_seq[0:int(hit[11])]
            found_best_V = True
        elif J_pattern.search(segment) and not found_best_J:
            J_locus_key = "TR{}J".format(segment[2])
            best_J_name = segment
            best_J_seq = IMGT_seqs[J_locus_key][segment]
            # hit 10 is the start of the J sequence
            best_J_seq = best_J_seq[int(hit[10]) - 1:]
            found_best_J = True

    junction = []

    parens_pattern = re.compile(r"\([CAGT]+\)")

    if locus == "B" or locus == "D":
        # junc_seqs = junction_details[1:3]
        VD_junc = junction_details[1]
        D_region = junction_details[2]
        DJ_junc = junction_details[3]
        if parens_pattern.search(VD_junc):
            VD_junc = re.sub(r'[\(\)]', '', VD_junc)
            length_in_parens = len(VD_junc)
            best_V_seq = best_V_seq[: -length_in_parens]
        if parens_pattern.search(DJ_junc):
            DJ_junc = re.sub(r'[\(\)]', '', DJ_junc)
            length_in_parens = len(DJ_junc)
            best_J_seq = best_J_seq[length_in_parens:]
        junc_seqs = [VD_junc, D_region, DJ_junc]


    elif locus == "A" or locus == "G":
        VJ_junc = junction_details[1]
        # junctions in parentheses are represented in the coordinates of the matched segments.
        # Need to trim them then include the NTs in the junction
        if parens_pattern.search(VJ_junc):
            VJ_junc = re.sub(r'[\(\)]', '', VJ_junc)
            length_in_parens = len(VJ_junc)
            best_V_seq = best_V_seq[: -length_in_parens]
            best_J_seq = best_J_seq[length_in_parens:]
        junc_seqs = [VJ_junc]

    for seq in junc_seqs:
        seq = re.sub(r'[\(\)]', '', seq)
        if seq != "N/A":
            junction.append(seq)

    junction = "".join(junction)

    constant_seq = constant_seqs[locus]

    # Editing IMGT V and J sequences to include any alterations from the junction details
    V_end_seq = junction_details[0]
    J_start_seq = junction_details[-1]
    best_V_seq = best_V_seq[:-(len(V_end_seq))]
    best_V_seq = best_V_seq + V_end_seq
    best_J_seq = best_J_seq[len(J_start_seq):]
    best_J_seq = J_start_seq + best_J_seq

    full_rearrangement = best_V_seq + junction + best_J_seq + constant_seq
    productive_rearrangement = is_rearrangement_productive(best_V_seq + junction + best_J_seq + constant_seq[0:2])
    # fasta_line = ">chr={}__TCR{}_{}\n{}\n".format(sample_name, locus, query_name, full_rearrangement)

    bestVJ = [best_V_name, best_J_name]

    return (full_rearrangement, productive_rearrangement, bestVJ)


def is_rearrangement_productive(seq):
    # returns a tuple of three true/false values (productive, contains stop, in-frame)
    seq_mod_3 = len(seq) % 3
    if seq_mod_3 == 0:
        in_frame = True
    else:
        in_frame = False

    seq = Seq(seq, IUPAC.unambiguous_dna)
    aa_seq = seq.translate()
    contains_stop = "*" in aa_seq

    if in_frame and not contains_stop:
        productive = True
    else:
        productive = False

    return (productive, contains_stop, in_frame)


def get_segment_name(name, pattern):
    match = pattern.search(name)
    number = match.group(1)
    if match.group(3):
        sub_number = match.group(3)
    else:
        sub_number = ""
    return (number)


def get_fasta_line_for_contig_assembly(trinity_seq, hit_table, locus, IMGT_seqs, sample_name, query_name):
    found_best_V = False
    found_best_D = False
    found_best_J = False

    V_pattern = re.compile(r"TR[ABGD]V\d")
    D_pattern = re.compile(r"TR[BD]D\d")
    J_pattern = re.compile(r"TR[ABGD]J\d")

    for hit in hit_table:
        segment = hit[2]
        if V_pattern.search(segment) and not found_best_V:
            V_locus_key = "TR{}V".format(segment[2])
            best_V_name = segment
            # Remove forward slashes from shared A/D gene names to be the same as in the IMGT files.
            segment = segment.replace("/", "_")
            ref_V_seq = IMGT_seqs[V_locus_key][segment]

            # hit[11] is the end of the V sequence
            # best_V_seq = best_V_seq[0:int(hit[11])]
            found_best_V = True
        elif J_pattern.search(segment) and not found_best_J:
            J_locus_key = "TR{}J".format(segment[2])
            best_J_name = segment
            ref_J_seq = IMGT_seqs[J_locus_key][segment]
            # hit 10 is the start of the J sequence
            # best_J_seq = best_J_seq[int(hit[10])-1 :]
            found_best_J = True

    # work out if sequence that exists is in frame

    found_V = False
    found_J = False
    for entry in hit_table:
        if entry[0] == 'V':
            if not found_V:
                ref_V_start = int(entry[10])
                found_V = True
        if entry[0] == 'J':
            if not found_J:
                ref_J_end = int(entry[11])
                found_J = True
    start_padding = ref_V_start - 1

    ref_J_length = len(ref_J_seq)
    end_padding = (ref_J_length - ref_J_end)

    full_effective_length = start_padding + len(
        trinity_seq) + end_padding + 2  # add two because need first two bases of constant region to put in frame

    if full_effective_length % 3 == 0:
        in_frame = True
    else:
        in_frame = False

    # remove the minimal nucleotides from the trinity sequence to check for stop codons
    start_base_removal_count = (3 - (ref_V_start - 1)) % 3
    end_base_removal_count = (1 - end_padding) % 3

    seq = trinity_seq[start_base_removal_count:-(end_base_removal_count)]
    seq = Seq(seq, IUPAC.unambiguous_dna)
    aa_seq = seq.translate()
    contains_stop = "*" in aa_seq

    if in_frame and not contains_stop:
        productive = True
    else:
        productive = False

    productive_rearrangement = (productive, contains_stop, in_frame)

    bestVJ = [best_V_name, best_J_name]

    return (productive_rearrangement, bestVJ)


def collapse_close_sequences(recombinants, locus):
    # pdb.set_trace()
    contig_names = [r.contig_name for r in recombinants]
    filtered_contig_names = [r.contig_name for r in recombinants]
    uncollapsible_contigs = []
    if len(recombinants) > 1:
        for i in range(len(recombinants) - 1):
            base_name = recombinants[i].contig_name
            base_seq = recombinants[i].imgt_reconstructed_seq
            base_V_segment = recombinants[i].best_VJ_names[0]
            base_J_segment = recombinants[i].best_VJ_names[1]

            base_id = recombinants[i].identifier
            base_junc = base_id.split("_")[1]
            base_e_value = float(recombinants[i].hit_table[0][-2])

            for j in range(i + 1, len(recombinants)):
                comp_name = recombinants[j].contig_name
                comp_seq = recombinants[j].imgt_reconstructed_seq
                comp_V_segment = recombinants[j].best_VJ_names[0]
                comp_J_segment = recombinants[j].best_VJ_names[1]
                comp_id = recombinants[j].identifier
                comp_junc = comp_id.split("_")[1]
                comp_e_value = float(recombinants[j].hit_table[0][-2])
                lev_dist = Levenshtein.distance(base_seq, comp_seq)
                # print("{}\t{}\t{}".format(base_id, comp_id, lev_dist))
                if lev_dist < 35 and not base_id == comp_id and base_name in filtered_contig_names \
                        and comp_name in filtered_contig_names:
                    # pdb.set_trace()
                    # define re pattern here to find TRAVx[DN] or TRDVx[DN] depending on locus
                    if locus == "TCRA":
                        duplicate_pattern = re.compile(r"TRAV\d+[DN]")
                        segment_pattern = re.compile(r"TRAV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    elif locus == "TCRD":
                        duplicate_pattern = re.compile(r"DV\d+[DN]")
                        segment_pattern = re.compile(r"DV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    else:
                        uncollapsible_contigs.append("{}_vs_{}".format(base_name, comp_name))
                        attempt_collapse = False
                    if attempt_collapse and (
                        duplicate_pattern.search(base_V_segment) or duplicate_pattern.search(comp_V_segment)):
                        base_segment = get_segment_name(base_V_segment, segment_pattern)
                        comp_segment = get_segment_name(comp_V_segment, segment_pattern)
                        if base_segment == comp_segment:
                            # find alignment with lowest E value for V match
                            if base_e_value <= comp_e_value:
                                filtered_contig_names.remove(comp_name)
                            else:
                                filtered_contig_names.remove(base_name)
                        else:
                            uncollapsible_contigs.append("{}_vs_{}".format(base_name, comp_name))

                    else:
                        uncollapsible_contigs.append("{}_vs_{}".format(base_name, comp_name))

                elif lev_dist < 75 and not base_id == comp_id and base_name in filtered_contig_names \
                        and comp_name in filtered_contig_names:
                    if locus == "TCRA":
                        duplicate_pattern = re.compile(r"TRAV\d+[DN]")
                        segment_pattern = re.compile(r"TRAV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    elif locus == "TCRD":
                        duplicate_pattern = re.compile(r"DV\d+[DN]")
                        segment_pattern = re.compile(r"DV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    else:
                        uncollapsible_contigs.append("{}_vs_{}".format(base_name, comp_name))
                        attempt_collapse = False
                    if attempt_collapse and (
                        duplicate_pattern.search(base_V_segment) or duplicate_pattern.search(comp_V_segment)):
                        base_segment = get_segment_name(base_V_segment, segment_pattern)
                        comp_segment = get_segment_name(comp_V_segment, segment_pattern)
                        if (base_segment == comp_segment) and (base_junc == comp_junc) and (
                            base_J_segment == comp_J_segment):
                            # find alignment with lowest E value for V match
                            if base_e_value <= comp_e_value:
                                filtered_contig_names.remove(comp_name)
                            else:
                                filtered_contig_names.remove(base_name)
                        else:
                            uncollapsible_contigs.append("{}_vs_{}".format(base_name, comp_name))

                elif base_id == comp_id and base_name in filtered_contig_names and comp_name in filtered_contig_names:
                    if base_e_value <= comp_e_value:
                        filtered_contig_names.remove(comp_name)
                    else:
                        filtered_contig_names.remove(base_name)

    recombinants_to_delete = []

    for r in recombinants:
        if not r.contig_name in filtered_contig_names:
            recombinants_to_delete.append(r)

    [recombinants.remove(r) for r in recombinants_to_delete]

    return (recombinants)


def load_kallisto_counts(tsv_file):
    counts = {'A': {}, 'B': {}, 'G': {}, 'D': {}}
    with open(tsv_file) as tsvh:
        for line in tsvh:
            if "TCR" in line:
                line = line.rstrip()
                line = line.split("\t")
                locus = line[0].split("|")[-1].split("_")[-1][2]
                name = line[0].split("|")[1]
                tpm = float(line[4])
                counts[locus][name] = tpm
    return counts


def make_cell_network_from_dna(cells, colorscheme, colours, keep_unlinked, shape, dot, neato):
    G = nx.MultiGraph()
    # initialise all cells as nodes

    if shape == 'circle':
        for cell in cells:
            G.add_node(cell, shape=shape, label=cell.html_style_label_for_circles(), sep=0.4, fontname="helvetica neue")
    else:
        for cell in cells:
            G.add_node(cell, shape=shape, label=cell.html_style_label_dna(), fontname="helvetica neue")
    # make edges:
    for i in range(len(cells)):
        current_cell = cells[i]
        comparison_cells = cells[i + 1:]

        for locus in ['A', 'B', 'D', 'G']:
            col = colours[locus]

            # current_identifiers = current_cell.getMainRecombinantIdentifiersForLocus(locus)
            for comparison_cell in comparison_cells:
                shared_identifiers = 0
                if current_cell.all_recombinants[locus] is not None:
                    for current_recombinant in current_cell.all_recombinants[locus]:
                        current_id_set = current_recombinant.all_poss_identifiers
                        if comparison_cell.all_recombinants[locus] is not None:
                            for comparison_recombinant in comparison_cell.all_recombinants[locus]:
                                comparison_id_set = comparison_recombinant.all_poss_identifiers
                                if len(current_id_set.intersection(comparison_id_set)) > 0:
                                    shared_identifiers += 1

                # comparison_identifiers = comparison_cell.getAllRecombinantIdentifiersForLocus(locus)
                # common_identifiers = current_identifiers.intersection(comparison_identifiers)
                if shared_identifiers > 0:
                    width = shared_identifiers * 2
                    G.add_edge(current_cell, comparison_cell, locus, penwidth=width, color=col,
                               weight=shared_identifiers, colorscheme=colorscheme)

    deg = G.degree()

    to_remove = [n for n in deg if deg[n] == 0]

    if len(to_remove) < len(G.nodes()):
        if not shape == 'circle':
            G.remove_nodes_from(to_remove)
            drawing_tool = [dot, '-Gsplines=true', '-Goverlap=false', '-Gsep=0.4']
        else:
            drawing_tool = [dot, '-Gsplines=true', '-Goverlap=false']
    else:
        drawing_tool = [neato, '-Gsplines=true', '-Goverlap=false']

    bgcolors = ['#8dd3c720', '#ffffb320', '#bebada20', '#fb807220', '#80b1d320', '#fdb46220', '#b3de6920', '#fccde520',
                '#d9d9d920', '#bc80bd20', '#ccebc520', '#ffed6f20']
    component_counter = 0
    component_groups = list()
    j = 0
    components = nx.connected_components(G)

    for component in components:
        members = list()
        if len(component) > 1:
            for cell in component:
                members.append(cell.name)
                G.node[cell]['style'] = 'filled'
                G.node[cell]['fillcolor'] = bgcolors[j]
                cell.bgcolor = bgcolors[j]

            if j < 11:
                j += 1
            else:
                component_counter += 1
                j = 0

        component_groups.append(members)

    return (G, drawing_tool, component_groups)


def draw_network_from_cells(cells, output_dir, output_format, dot, neato, draw_graphs):
    cells = list(cells.values())
    colorscheme = 'set15'
    colours = {'A': '1', 'B': '2', 'G': '3', 'D': '5', 'mean_both': '#a8a8a8bf'}
    network, draw_tool, component_groups = make_cell_network_from_dna(cells, colorscheme, colours, False, "box", dot,
                                                                      neato)
    network_file = "{}/clonotype_network_with_identifiers.dot".format(output_dir)
    try:
        nx.write_dot(network, network_file)
    except AttributeError:
        import pydotplus
        nx.drawing.nx_pydot.write_dot(network, network_file)
    if draw_graphs:
        command = draw_tool + ['-o', "{output_dir}/clonotype_network_with_identifiers.{output_format}".format(
            output_dir=output_dir, output_format=output_format), "-T", output_format, network_file]
        subprocess.check_call(command)

    network, draw_tool, cgx = make_cell_network_from_dna(cells, colorscheme, colours, False, "circle", dot, neato)
    network_file = "{}/clonotype_network_without_identifiers.dot".format(output_dir)
    try:
        nx.write_dot(network, network_file)
    except AttributeError:
        import pydotplus
        nx.drawing.nx_pydot.write_dot(network, network_file)
    if draw_graphs:
        command = draw_tool + ['-o', "{output_dir}/clonotype_network_without_identifiers.{output_format}".format(
            output_dir=output_dir, output_format=output_format), "-T", output_format, network_file]
        subprocess.check_call(command)
    return (component_groups)


def get_component_groups_sizes(cells):
    cells = list(cells.values())
    G = nx.MultiGraph()
    # initialise all cells as nodes
    for cell in cells:
        G.add_node(cell)
    # make edges:
    for i in range(len(cells)):
        current_cell = cells[i]
        comparison_cells = cells[i + 1:]

        for locus in ['A', 'B', 'D', 'G']:

            # current_identifiers = current_cell.getMainRecombinantIdentifiersForLocus(locus)
            for comparison_cell in comparison_cells:
                shared_identifiers = 0
                if current_cell.all_recombinants[locus] is not None:
                    for current_recombinant in current_cell.all_recombinants[locus]:
                        current_id_set = current_recombinant.all_poss_identifiers
                        if comparison_cell.all_recombinants[locus] is not None:
                            for comparison_recombinant in comparison_cell.all_recombinants[locus]:
                                comparison_id_set = comparison_recombinant.all_poss_identifiers
                                if len(current_id_set.intersection(comparison_id_set)) > 0:
                                    shared_identifiers += 1

                # comparison_identifiers = comparison_cell.getAllRecombinantIdentifiersForLocus(locus)
                # common_identifiers = current_identifiers.intersection(comparison_identifiers)
                if shared_identifiers > 0:
                    width = shared_identifiers * 2
                    G.add_edge(current_cell, comparison_cell, locus, penwidth=width, weight=shared_identifiers)

    deg = G.degree()

    to_remove = [n for n in deg if deg[n] == 0]

    # if len(to_remove) < len(G.nodes()):
    #    G.remove_nodes_from(to_remove)

    components = nx.connected_components(G)

    component_groups = list()

    singlets = []
    for component in components:
        members = list()
        if len(component) > 1:
            for cell in component:
                members.append(cell.name)
            component_groups.append(members)
        else:
            for cell in component:
                singlets.append(cell.name)

    clonotype_size_counter = Counter([len(x) for x in component_groups])
    clonotype_size_counter.update({1: len(singlets)})

    clonotype_sizes = []
    max_size = max(list(clonotype_size_counter.keys()))
    if max_size < 5:
        for x in range(1, max_size + 1):
            clonotype_sizes.append(clonotype_size_counter[x])
        zero_padding = 5 - len(clonotype_sizes)
        clonotype_sizes = clonotype_sizes + [0] * zero_padding
    else:
        for x in range(1, max_size + 1):
            clonotype_sizes.append(clonotype_size_counter[x])

    return (clonotype_sizes)


def check_config_file(filename):
    if not os.path.isfile(filename):
        print()
        print("Couldn't find config file: {}".format(filename))
        print()
        exit(1)


def bowtie2_alignment(bowtie2, ncores, locus_names, output_dir, cell_name, synthetic_genome_path, fastq1,
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


def assemble_with_trinity(trinity, locus_names, output_dir, cell_name, ncores, trinity_grid_conf, JM,
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

    # if len(unsuccessful_directories) == 2:

    return successful_files


def run_IgBlast(igblast, locus_names, output_dir, cell_name, index_location, ig_seqtype, species,
                should_resume):
    print("##Running IgBLAST##")

    species_mapper = {
        'Mmus': 'mouse',
        'Hsap': 'human'
    }

    igblast_species = species_mapper[species]

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
            with open(igblast_out, 'w') as out:
                # print(" ").join(pipes.quote(s) for s in command)
                subprocess.check_call(command, stdout=out, stderr=DEVNULL)

    DEVNULL.close()


def quantify_with_kallisto(kallisto, cell, output_dir, cell_name, kallisto_base_transcriptome, fastq1, fastq2,
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

