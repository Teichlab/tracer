
############################################################################## 
#                         Functions for use with                             #
# TraCeR - a tool to reconstruct TCR sequences from single-cell RNA-seq data #    
#                                                                            #
# Please see README and LICENCE for details of use and licence conditions.   #
# This software was written by Mike Stubbington (mstubb@ebi.ac.uk) from the  #
# Teichmann Lab, EMBL-EBI (www.teichlab.org). Latest versions are available  #
# for download at www.github.com/teichlab/tracer.                            #
#                                                                            #
#      Copyright (c) 2015 EMBL - European Bioinformatics Institute           # 
############################################################################## 

import sys
import os
import subprocess
import pipes
import glob
import shutil
import re
from collections import defaultdict, Counter
from time import sleep
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna
import pdb
import Levenshtein

##CLASSES##
class Cell:
    'Class to describe T cells containing A and B loci'
    
    def __init__(self, cell_name, A_recombinants, B_recombinants, G_recombinants, D_recombinants, expt_label):
        self.name = cell_name
        self.A_recombinants = A_recombinants
        self.B_recombinants = B_recombinants
        self.G_recombinants = G_recombinants
        self.D_recombinants = D_recombinants
        self.bgcolor = None
        self.all_recombinants = {'A' : A_recombinants, 'B' : B_recombinants, 'G' : G_recombinants, 'D' : D_recombinants}
        self.cdr3_comparisons = {'A' : None, 'B' : None, 'mean_both' : None}
        self.expt_label = expt_label
    
    def reset_cdr3_comparisons(self):
            self.cdr3_comparisons = {'A' : None, 'B' : None, 'mean_both' : None}
        
    def getAllRecombinantIdentifiersForLocus(self, locus):
        recombinants = self.all_recombinants[locus]
        identifier_list = set()
        if recombinants is not None:
            for recombinant in recombinants:
                all_possible_recombinant_identifiers = recombinant.all_poss_identifiers
                for identifier in all_possible_recombinant_identifiers:
                    identifier_list.add(identifier)
        return(identifier_list)
    
    def getMainRecombinantIdentifiersForLocus(self, locus):
        recombinants = self.all_recombinants[locus]
        identifier_list = set()
        if recombinants is not None:
            for recombinant in recombinants:
               identifier_list.add(recombinant.identifier)
        return(identifier_list)
    
    def getAllRecombinantCDR3ForLocus(self, locus):
        recombinants = self.all_recombinants[locus]
        identifier_list = set()
        if recombinants is not None:
            for recombinant in recombinants:
                cdr3 = str(recombinant.cdr3)
                if "Couldn't" not in cdr3:
                    identifier_list.add(cdr3)
        return(identifier_list)
   
    
    def html_style_label_dna(self, both_names):
        colours = {'A' : {'productive' : '#E41A1C', 'non-productive' : "#ff8c8e"}, 'B' : {'productive' : '#377eb8', 'non-productive' : "#95c1e5"}, 'G' : {'productive' : '#4daf4a', 'non-productive' : "#aee5ac"}, 'D' : {'productive' : '#984ea3', 'non-productive' : "#deace5"}}
        styles = {'Yes' : ('<I>', '</I>'), 'No' : ('', '')}
        locus_names = ['A','B','G','D']
        recombinants = dict()
        final_string = '<<FONT POINT-SIZE="16"><B>' + self.name + "</B></FONT>"
        for locus, recombinant_list in self.all_recombinants.iteritems():
            recombinant_set = set()
            if recombinant_list is not None:
                for recombinant in recombinant_list:
                    if recombinant.productive:
                        prod = "productive"
                    else:
                        prod = "non-productive"
                    if both_names:
                        recombinant_set.add("<BR/>" + '<FONT COLOR = "{}">'.format(colours[locus][prod]) + styles[recombinant.recovered_from_filter][0] + "{}<BR/>{}".format(recombinant.identifier, str(recombinant.cdr3)) + styles[recombinant.recovered_from_filter][1] + '</FONT>')
                    else:
                        recombinant_set.add("<BR/>" + '<FONT COLOR = "{}">'.format(colours[locus][prod]) + styles[recombinant.recovered_from_filter][0] + recombinant.identifier + styles[recombinant.recovered_from_filter][1] + '</FONT>')
                
                recombinants[locus] = recombinant_set
        for locus in locus_names:
            if locus in recombinants.keys():
                id_string = "".join(recombinants[locus])
                final_string = final_string + id_string
        final_string = final_string + ">"
        return(final_string)
        #return(self.name)
    
    def html_style_label_cdr3(self, both_names):
        colours = {'A' : {'productive' : '#E41A1C', 'non-productive' : "#ff8c8e"}, 'B' : {'productive' : '#377eb8', 'non-productive' : "#95c1e5"}, 'G' : {'productive' : '#4daf4a', 'non-productive' : "#aee5ac"}, 'D' : {'productive' : '#984ea3', 'non-productive' : "#deace5"}}
        styles = {'Yes' : ('<I>', '</I>'), 'No' : ('', '')}
        locus_names = ['A','B','G','D']
        recombinants = dict()
        final_string = '<<FONT POINT-SIZE="16"><B>' + self.name + "</B></FONT>"
        for locus, recombinant_list in self.all_recombinants.iteritems():
            recombinant_set = set()
            if recombinant_list is not None:
                for recombinant in recombinant_list:
                    if recombinant.productive:
                        prod = "productive"
                    else:
                        prod = "non-productive"
                    if both_names:
                        recombinant_set.add("<BR/>" + '<FONT COLOR = "{}">'.format(colours[locus][prod]) + styles[recombinant.recovered_from_filter][0] + "{}<BR/>{}".format(recombinant.identifier, str(recombinant.cdr3)) + styles[recombinant.recovered_from_filter][1] + '</FONT>')
                    else:
                        recombinant_set.add("<BR/>" + '<FONT COLOR = "{}">'.format(colours[locus][prod]) + styles[recombinant.recovered_from_filter][0] + str(recombinant.cdr3) + styles[recombinant.recovered_from_filter][1] + '</FONT>')
                
                recombinants[locus] = recombinant_set
        for locus in locus_names:
            if locus in recombinants.keys():
                id_string = "".join(recombinants[locus])
                final_string = final_string + id_string
        final_string = final_string + ">"
        return(final_string)
        #return(self.name)
    
    def __str__(self):
        return(self.name)
    
    def full_description(self):
        #pdb.set_trace()
        return_list = [self.name, '#TCRA#']
        
        if not self.A_recombinants is None:
            for recombinant in self.A_recombinants:
                return_list.append(str(recombinant))
        else:
            return_list.append("No TCRA recombinants")
        
        return_list.append('\n#TCRB#')
        if not self.B_recombinants is None:
            for recombinant in self.B_recombinants:
                return_list.append(str(recombinant))
        else:
            return_list.append("No TCRB recombinants")
        
        return_list.append('\n#TCRG#')
        if not self.G_recombinants is None:
            for recombinant in self.G_recombinants:
                return_list.append(str(recombinant))
        else:
            return_list.append("No TCRG recombinants")
            
        return_list.append('\n#TCRD#')
        if not self.D_recombinants is None:
            for recombinant in self.D_recombinants:
                return_list.append(str(recombinant))
        else:
            return_list.append("No TCRD recombinants")
            
            
        return("\n".join(return_list))
    
    def get_fasta_string(self):
        seq_string = []
        for locus, recombinants in self.all_recombinants.iteritems():
            if recombinants is not None:
                for rec in recombinants:
                    name = ">{contig_name}|{identifier}".format(contig_name=rec.contig_name, identifier=rec.identifier)
                    seq = rec.dna_seq
                    seq_string.append("\n".join([name, seq]))
        return("\n".join(seq_string + ["\n"]))    
        
class Recombinant:
    'Class to describe a recombined TCR locus as determined from the single-cell pipeline'
    def __init__(self, contig_name, locus, identifier, all_poss_identifiers, productive, stop_codon, in_frame, read_count, dna_seq, hit_table, summary, junction_details, best_VJ_names, alignment_summary):
        self.contig_name = contig_name
        self.locus = locus
        self.identifier = identifier
        self.all_poss_identifiers = all_poss_identifiers
        self.productive = productive
        self.read_count = read_count
        self.dna_seq = dna_seq
        self.cdr3 = self._get_cdr3(dna_seq)
        self.hit_table = hit_table
        self.summary = summary
        self.junction_details = junction_details
        self.best_VJ_names = best_VJ_names
        self.alignment_summary = alignment_summary
        self.in_frame = in_frame
        self.stop_codon = stop_codon

    def __str__(self):
        return("{} {} {} {}".format(self.identifier, self.productive, self.read_count, self.recovered_from_filter))
    
    def _get_cdr3(self, dna_seq):
        aaseq = Seq(str(dna_seq), generic_dna).translate()
        if re.findall('FG.G',str(aaseq)) and re.findall('C',str(aaseq)):
                        indices = [i for i, x in enumerate(aaseq) if x == 'C']
                        upper = str(aaseq).find(re.findall('FG.G',str(aaseq))[0])
                        for i in indices:
                            if i < upper:
                                lower = i
                        cdr3 = aaseq[lower:upper+4]
        elif re.findall('FG.G',str(aaseq)):
            cdr3 = "Couldn't find conserved cysteine"
        elif re.findall('C',str(aaseq)):
            cdr3 = "Couldn't find FGXG"
        else:
            cdr3 = "Couldn't find either conserved boundary"
        return(cdr3)



##FUNCTIONS##

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

    return(sorted(return_list))

def get_filename_and_locus(name):
    pattern = re.compile(r"(.+)_TCR([ABDG])")
    pattern_match = pattern.search(name)
    file = pattern_match.group(1)
    locus = pattern_match.group(2)
    return([file, locus])

def sort_locus_names(dictionary_to_sort):
    for key, value in dictionary_to_sort.iteritems():
        sorted_value = sorted(value)
        dictionary_to_sort[key] = sorted_value
    return(dictionary_to_sort)

def load_IMGT_seqs(file):
    seqs = {}
    for record in SeqIO.parse(open(file, 'rU'), 'fasta'):
        seqs[record.id] = str(record.seq)
    return(seqs)

def parse_IgBLAST(locus_names, output_dir, cell_name, expt_label, imgt_seq_location):
    segment_names = ['TRAJ', 'TRAV', 'TRBD', 'TRBJ', 'TRBV']
    IMGT_seqs = dict()
    for segment in segment_names:
        IMGT_seqs[segment] = load_IMGT_seqs("{}/{}.fa".format(imgt_seq_location, segment))
    
    all_locus_data = defaultdict(dict)
    for locus in locus_names:
        file = "{output_dir}/IgBLAST_output/{cell_name}_{locus}.IgBLASTOut".format(output_dir=output_dir, cell_name=cell_name, locus=locus)
        
        igblast_result_chunks = split_igblast_file(file)
        
        
        for chunk in igblast_result_chunks:
            (query_name, chunk_details) = process_chunk(chunk)
            
            all_locus_data[locus][query_name] = chunk_details
        
    cell = find_possible_alignments(all_locus_data, locus_names, cell_name, IMGT_seqs, expt_label)    
    return(cell)
    



def split_igblast_file(filename):
    # code adapted from http://stackoverflow.com/questions/19575702/pythonhow-to-split-file-into-chunks-by-the-occurrence-of-the-header-word
    token = '# IGBLASTN'
    chunks = []
    current_chunk = []
    
    for line in open(filename):
        line=line.rstrip()
        if line.startswith(token) and current_chunk: 
           # if line starts with token and the current chunk is not empty
           chunks.append(current_chunk[:]) #  add not empty chunk to chunks
           current_chunk = [] #  make current chunk blank
        # just append a line to the current chunk on each iteration
        current_chunk.append(line)

    chunks.append(current_chunk)  #  append the last chunk outside the loop
    return(chunks)

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
            #return_dict['query_name'] = query_name
            
        elif line_x.startswith('# V-(D)-J rearrangement summary'):
            store_VDJ_rearrangement_summary = True
        
        elif line_x.startswith('# V-(D)-J junction details'):
            store_junction_details = True
        
        elif line_x.startswith('# Alignment summary'):
            store_alignment_summary = True
        
        elif line_x.startswith('# Hit table'):
            store_hit_table = True
    
         


    return (query_name, return_dict)


def find_possible_alignments(sample_dict, locus_names, cell_name, IMGT_seqs, expt_label):
    #pdb.set_trace()
    alignment_dict = defaultdict(dict)
    recombinants = {'TCRA':[], 'TCRB':[]}
    for locus in locus_names:
        data_for_locus = sample_dict[locus]
        for query_name, query_data in data_for_locus.iteritems():
            #pdb.set_trace()
            processed_hit_table = process_hit_table(query_name, query_data, locus)
            
            if processed_hit_table is not None:
                (returned_locus, good_hits, rearrangement_summary) = processed_hit_table
                junction_list = query_data['junction_details']
                
                (fasta_line_for_contig, is_productive, bestVJNames) = get_fasta_line_for_contig(rearrangement_summary, junction_list, good_hits, returned_locus, IMGT_seqs, cell_name, query_name)

                
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
                
                #pdb.set_trace()
                rec = Recombinant(contig_name=query_name, locus=returned_locus, identifier=identifier, all_poss_identifiers=all_poss_identifiers, productive=is_productive[0], stop_codon=is_productive[1], in_frame=is_productive[2], read_count=0.0, dna_seq=fasta_line_for_contig, hit_table=good_hits, summary=rearrangement_summary, junction_details=junction_list, best_VJ_names=bestVJNames, alignment_summary=alignment_summary)
                recombinants[locus].append(rec)
              
    #pdb.set_trace()
    if recombinants:
        for locus, rs in recombinants.iteritems():
                #Adding code to collapse sequences with very low Levenshtein distances caused by confusion between TRAVxD and TRAVx segments with different alignment lengths from IgBlast.     
                recombinants[locus] = collapse_close_sequences(rs, locus)
        cell = Cell(cell_name, recombinants['TCRA'], recombinants['TCRB'], None, None, expt_label)
    else:
        cell = Cell(cell_name, None, None, None, None, expt_label)
    
    #pdb.set_trace()
    return(cell)

def remove_NA(junc_string):
    new_string = junc_string.replace("N/A", "")
    return(new_string)


def remove_allele_stars(segment):
    p = re.compile(r"(.+)\*\d+")
    m = p.search(segment)
    return(m.group(1))    

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
            return("A", good_hits, rearrangement_summary)
        else:
            return(None)
        
    #elif locus == "D":
    #    if "D" in found_V and "D" in found_J:
    #        return("D", good_hits, rearrangement_summary)
    #    else:
    #        return(None)
            
    elif locus == "TCRB":
        if "TCRB" in found_V and "TCRB" in found_J:
            return("B", good_hits, rearrangement_summary)
        else:
            return(None)
    
    #elif locus == "G":
    #    if "G" in found_V and "G" in found_J:
    #        return("G", good_hits, rearrangement_summary)
    #    else:
    #        return(None)

def get_fasta_line_for_contig(rearrangement_summary, junction_details, hit_table, locus, IMGT_seqs, sample_name, query_name):
    constant_seqs = dict()
    constant_seqs["A"] = "ACATCCAGAACCCAGAACCTGCTGTGTACCAGTTAAAAGATCCTCGGTCTCAGGACAGCACCCTCTGCCTGTTCACCGACTTTGACTCCCAAATCAATGTGCCGAAAACCATGGAATCTGGAACGTTCATCACTGACAAAACTGTGCTGGACATGAAAGCTATGGATTCCAAGAGCAATGGGGCCATTGCCTGGAGCAACCAGACAAGCTTCACCTGCCAAGATATCTTCAAAGAGACCAACGCCACCTACCCCAGTTCAGACGTTCCCTGTGATGCCACGTTGACTGAGAAAAGCTTTGAAACAGATATGAACCTAAACTTTCAAAACCTGTCAGTTATGGGACTCCGAATCCTCCTGCTGAAAGTAGCCGGATTTAACCTGCTCATGACGCTGAGGCTGTGGTCCAGTTGA"
    #use first 258 bases of TRBC because they're the same between C1 and C2
    constant_seqs["B"] = "AGGATCTGAGAAATGTGACTCCACCCAAGGTCTCCTTGTTTGAGCCATCAAAAGCAGAGATTGCAAACAAACAAAAGGCTACCCTCGTGTGCTTGGCCAGGGGCTTCTTCCCTGACCACGTGGAGCTGAGCTGGTGGGTGAATGGCAAGGAGGTCCACAGTGGGGTCAGCACGGACCCTCAGGCCTACAAGGAGAGCAATTATAGCTACTGCCTGAGCAGCCGCCTGAGGGTCTCTGCTACCTTCTGGCACAATCCTC"
    constant_seqs["D"] = "AAAGCCAGCCTCCGGCCAAACCATCTGTTTTCATCATGAAAAATGGAACAAATGTTGCTTGTCTGGTGAAAGATTTCTACCCTAAAGAGGTGACTATAAGTCTCAGATCATCCAAGAAGATTGTGGAATTCGACCCTGCTATAGTCATCTCCCCCAGCGGGAAGTACAGTGCTGTCAAGCTTGGTCAGTATGGAGATTCGAATTCAGTGACATGTTCAGTTCAGCACAACAGTGAAACTGTGCACTCGACTGACTTTGAACCATATGCAAATTCTTTCAATAATGAAAAACTACCAGAACCTGAAAATGACACACAAATTTCAGAGCCTTGCTATGGCCCAAGAGTCACAGTTCACACTGAGAAGGTAAACATGATGTCCCTCACGGTGCTGGGCCTACGACTGCTGTTTGCCAAGACCATTGCCATCAATTTTCTCTTGACTGTTAAGTTATTCTTTTAAGGGTGGGCTGACATGAGGAGACTACGGTTCCTGAAAGAAATCAAAAGCTTAGAAAGATGCTATATCCCAGGCTTCCAACTTCTCAGTGCTTCAGACTGACCCTTCACCACCACATTTAAACAGCTGCTAACAAAACCAGCTTTTCTGTGACAGCAACAAGCCTAGCTAATCCTCCAGTCTAGAAGAAAAGCAAAAGCCCTCGGGACCCCCGGCTTTACCTGCTGCTTTATAAAGGCATGGGAAGTTATGAAAACAGATCCATTTTATTTTGCCCCCATAATTGGTATACTTTGAAAATGGTGTTTCATCCTTCTTCATTTACCCAGAACTAGGAAGTGGGGACCAGCTTCATTATCCAGGAGGAAATAATCTTGAGAGAGAGAACCCGTATCTTTTTAGCTAAACATGGAAAGCTGTACTCAACTCATCCCTAGCCAGAGCCCCCTCCTCCTCTCCTGAGGCGAGCATGGCCCAGCCCCCCCCCCTTTGTATTTACTCCAATAGTCACACAGGAGAGTTTTCCTAGCAGCACTACGGTGTGAACAATTTTAGCACTTTCTGTTTCTCCTAATACTTTACAAACAAACTCACACTTGGCTTCCTTAATGCTCTCCAAGCAGACAATAAAGCTTCTAAGATCGCATC"
    #for TRGC use first 150 bases. Found by aligning the 4 C region transcripts and taking consensus. Ignored start of TCRG-C4-201 because it's only in that one. 
    constant_seqs["G"] = "GACAAAAGGCTTGATGCAGACATTTCCCCCAAGCCCACTATTTTCCTTCCTTCTGTTGCTGAAACAAATCTCCATAAGACTGGGACATACCTTTGTCTCCTTGAAAAGTTCTTTCCCGATGTCATAAGGGTGTATTGGAAAGAAAAGG"
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
            segment = segment.replace("/", "_") #remove forward slashes from shared A/D gene names to be the same as in the IMGT files.
            best_V_seq = IMGT_seqs[V_locus_key][segment]

            #hit[11] is the end of the V sequence
            best_V_seq = best_V_seq[0:int(hit[11])]
            found_best_V = True
        elif J_pattern.search(segment) and not found_best_J:
            J_locus_key = "TR{}J".format(segment[2])
            best_J_name = segment
            best_J_seq = IMGT_seqs[J_locus_key][segment]
            #hit 10 is the start of the J sequence
            best_J_seq = best_J_seq[int(hit[10])-1 :]
            found_best_J = True

    junction = []
    
    parens_pattern = re.compile(r"\([CAGT]+\)")
    
    if locus == "B" or locus == "D":
        #junc_seqs = junction_details[1:3]
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
            best_J_seq = best_J_seq[length_in_parens :]
        junc_seqs = [VD_junc, D_region, DJ_junc]

        
    elif locus == "A" or locus == "G":
        VJ_junc = junction_details[1]
        #junctions in parentheses are represented in the coordinates of the matched segments. Need to trim them then include the NTs in the junction
        if parens_pattern.search(VJ_junc):
            VJ_junc = re.sub(r'[\(\)]', '', VJ_junc)
            length_in_parens = len(VJ_junc)
            best_V_seq = best_V_seq[: -length_in_parens]
            best_J_seq = best_J_seq[length_in_parens :]
        junc_seqs = [VJ_junc]

    
    for seq in junc_seqs:
        seq = re.sub(r'[\(\)]', '', seq)
        if seq != "N/A":
            junction.append(seq)
    
    junction = "".join(junction)
    
    constant_seq = constant_seqs[locus]
    
    #editing IMGT V and J sequences to include any alterations from the junction details
    V_end_seq = junction_details[0]
    J_start_seq = junction_details[-1]
    best_V_seq = best_V_seq[:-(len(V_end_seq))]
    best_V_seq = best_V_seq + V_end_seq
    best_J_seq = best_J_seq[len(J_start_seq):]
    best_J_seq = J_start_seq + best_J_seq
    
    
    full_rearrangement = best_V_seq + junction + best_J_seq + constant_seq
    productive_rearrangement = is_rearrangement_productive(best_V_seq + junction + best_J_seq + constant_seq[0:2])
    #fasta_line = ">chr={}__TCR{}_{}\n{}\n".format(sample_name, locus, query_name, full_rearrangement)
    
    bestVJ = [best_V_name, best_J_name]
    
    return(full_rearrangement, productive_rearrangement, bestVJ)

def is_rearrangement_productive(seq):
    #returns a tuple of three true/false values (productive, contains stop, in-frame)
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
    
    return(productive, contains_stop, in_frame)
    
def collapse_close_sequences(recombinants, locus):
    
    contig_names = [r.contig_name for r in recombinants]
    filtered_contig_names = [r.contig_name for r in recombinants]
    uncollapsible_contigs = []
    if len(recombinants)>1:
        for i in range(len(recombinants)-1):
            base_name = recombinants[i].contig_name
            base_seq = recombinants[i].dna_seq
            base_V_segment = recombinants[i].best_VJ_names[0]
            base_id = recombinants[i].identifier
            base_e_value = float(recombinants[i].hit_table[0][-2])
            
            for j in range(i+1, len(recombinants)):
                comp_name = recombinants[j].contig_name
                comp_seq = recombinants[j].dna_seq
                comp_V_segment = recombinants[j].best_VJ_names[0]
                comp_id = recombinants[j].identifier
                comp_e_value = float(recombinants[j].hit_table[0][-2])
                lev_dist = Levenshtein.distance(base_seq, comp_seq)
                #print("{}\t{}\t{}".format(base_id, comp_id, lev_dist))
                if lev_dist < 35 and not base_id == comp_id and base_name in filtered_contig_names and comp_name in filtered_contig_names:   
                    #pdb.set_trace()             
                    #define re pattern here to find TRAVx[DN] or TRDVx[DN] depending on locus
                    if locus == "A":
                        duplicate_pattern = re.compile(r"TRAV\d+[DN]")
                        segment_pattern = re.compile(r"TRAV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    elif locus == "D":
                        duplicate_pattern = re.compile(r"DV\d+[DN]")
                        segment_pattern = re.compile(r"DV(\d+)([DN])?(-\d)?.+")
                        attempt_collapse = True
                    else:
                        uncollapsible_contigs.append("{}_vs_{}".format(base_name, comp_name))
                        attempt_collapse = False
                    if attempt_collapse and (duplicate_pattern.search(base_V_segment) or duplicate_pattern.search(comp_V_segment)):
                        base_segment = get_segment_name(base_V_segment, segment_pattern)
                        comp_segment = get_segment_name(comp_V_segment, segment_pattern)
                        if base_segment == comp_segment:
                            #find alignment with lowest E value for V match
                            if base_e_value <= comp_e_value:
                                filtered_contig_names.remove(comp_name)
                            else:
                                filtered_contig_names.remove(base_name) 
                        else:
                            uncollapsible_contigs.append("{}_vs_{}".format(base_name, comp_name))
                        
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
    
                        
    
    return(recombinants)        
                    