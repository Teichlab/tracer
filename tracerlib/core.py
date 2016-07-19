import re
from collections import Counter

import six
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq


class Cell(object):

    """Class to describe T cells containing A and B loci"""

    def __init__(self, cell_name, A_recombinants, B_recombinants, G_recombinants, D_recombinants, is_empty=False,
                 species="Mmus"):
        self.name = cell_name
        self.A_recombinants = A_recombinants
        self.B_recombinants = B_recombinants
        self.G_recombinants = G_recombinants
        self.D_recombinants = D_recombinants
        self.bgcolor = None
        self.all_recombinants = {'A': A_recombinants, 'B': B_recombinants, 'G': G_recombinants, 'D': D_recombinants}
        self.cdr3_comparisons = {'A': None, 'B': None, 'mean_both': None}
        self.is_empty = self._check_is_empty()
        self.is_inkt = self._check_if_inkt(species)

    def _check_is_empty(self):
        if (self.A_recombinants is None or len(self.A_recombinants) == 0) and (
                self.B_recombinants is None or len(self.B_recombinants) == 0):
            return (True)

    def _check_if_inkt(self, species):
        A_recombs = self.getMainRecombinantIdentifiersForLocus("A")
        inkt_ident = False
        for recomb in A_recombs:
            if species == "Mmus":
                if "TRAV11" in recomb and "TRAJ18" in recomb:
                    inkt_ident = recomb
            if species == "Hsap":
                if "TRAV10" in recomb and "TRAJ18" in recomb:
                    inkt_ident = recomb
        return (inkt_ident)

    def reset_cdr3_comparisons(self):
        self.cdr3_comparisons = {'A': None, 'B': None, 'mean_both': None}

    def getAllRecombinantIdentifiersForLocus(self, locus):
        recombinants = self.all_recombinants[locus]
        identifier_list = set()
        if recombinants is not None:
            for recombinant in recombinants:
                all_possible_recombinant_identifiers = recombinant.all_poss_identifiers
                for identifier in all_possible_recombinant_identifiers:
                    identifier_list.add(identifier)
        return (identifier_list)

    def getMainRecombinantIdentifiersForLocus(self, locus):
        recombinants = self.all_recombinants[locus]
        identifier_list = set()
        if recombinants is not None:
            for recombinant in recombinants:
                identifier_list.add(recombinant.identifier)
        return (identifier_list)

    def getAllRecombinantCDR3ForLocus(self, locus):
        recombinants = self.all_recombinants[locus]
        identifier_list = set()
        if recombinants is not None:
            for recombinant in recombinants:
                cdr3 = str(recombinant.cdr3)
                if "Couldn't" not in cdr3:
                    identifier_list.add(cdr3)
        return (identifier_list)

    def html_style_label_dna(self, transgenic_ids = False):
        colours = {'A': {'productive': '#E41A1C', 'non-productive': "#ff8c8e"},
                   'B': {'productive': '#377eb8', 'non-productive': "#95c1e5"},
                   'G': {'productive': '#4daf4a', 'non-productive': "#aee5ac"},
                   'D': {'productive': '#984ea3', 'non-productive': "#deace5"}}
        locus_names = ['A', 'B', 'G', 'D']
        recombinants = dict()
        final_string = '<<FONT POINT-SIZE="16"><B>' + self.name + "</B></FONT>"
        for locus, recombinant_list in six.iteritems(self.all_recombinants):
            recombinant_set = set()
            if recombinant_list is not None:
                for recombinant in recombinant_list:
                    if recombinant.productive:
                        prod = "productive"
                    else:
                        prod = "non-productive"
                    
                    if transgenic_ids and recombinant.identifier in transgenic_ids:
                        recombinant_set.add("<BR/>" + '<FONT COLOR = "#999999">'.format(colours[locus][prod])  + recombinant.identifier  + '</FONT>')
                        
                        
                        
                    else:    
                        recombinant_set.add("<BR/>" + '<FONT COLOR = "{}">'.format(colours[locus][prod])  + recombinant.identifier  + '</FONT>')
                    

                recombinants[locus] = recombinant_set
        for locus in locus_names:
            if locus in recombinants.keys():
                id_string = "".join(recombinants[locus])
                final_string = final_string + id_string
        final_string = final_string + ">"
        return (final_string)
        # return(self.name)

    def html_style_label_for_circles(self, transgenic_ids = False):
        colours = {'A': {'productive': '#E41A1C', 'non-productive': "#ff8c8e"},
                   'B': {'productive': '#377eb8', 'non-productive': "#95c1e5"},
                   'G': {'productive': '#4daf4a', 'non-productive': "#aee5ac"},
                   'D': {'productive': '#984ea3', 'non-productive': "#deace5"}}
        locus_names = ['A', 'B', 'G', 'D']
        recombinants = dict()
        final_string = '<<table cellspacing="6px" border="0" cellborder="0">'
        # final_string = "<"
        for locus, recombinant_list in six.iteritems(self.all_recombinants):
            recombinant_set = set()
            if recombinant_list is not None:
                for recombinant in recombinant_list:
                    if recombinant.productive:
                        prod = "productive"
                    else:
                        prod = "non-productive"
                    if transgenic_ids and recombinant.identifier in transgenic_ids:
                        recombinant_set.add('<tr><td height="10" width="40" border="2" color="{}"></td></tr>'.format(colours[locus][prod]))
                    else:
                        recombinant_set.add('<tr><td height="10" width="40" bgcolor="{}"></td></tr>'.format(colours[locus][prod]))
                    

                recombinants[locus] = recombinant_set
        strings = []
        for locus in locus_names:
            if locus in recombinants.keys():
                strings.append("".join(recombinants[locus]))

        id_string = "".join(strings)
        final_string = final_string + id_string
        final_string = final_string + "</table>>"
        return (final_string)

    def __str__(self):
        return (self.name)

    def full_description(self):
        # pdb.set_trace()
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

        return ("\n".join(return_list))

    def get_fasta_string(self):
        seq_string = []
        for locus, recombinants in six.iteritems(self.all_recombinants):
            if recombinants is not None:
                for rec in recombinants:
                    name = ">TCR|{contig_name}|{identifier}".format(contig_name=rec.contig_name,
                                                                    identifier=rec.identifier)
                    seq = rec.dna_seq
                    seq_string.append("\n".join([name, seq]))
        return ("\n".join(seq_string + ["\n"]))

    def summarise_productivity(self, locus):
        if self.all_recombinants[locus] is None:
            return ("0/0")
        else:
            recs = self.all_recombinants[locus]
            prod_count = 0
            total_count = len(recs)
            for rec in recs:
                if rec.productive:
                    prod_count += 1
            return ("{}/{}".format(prod_count, total_count))

    def filter_recombinants(self, transgenic_ids=False):
        for locus in ['A', 'B']:
            recs = self.all_recombinants[locus]
            if recs is not None:
                to_remove = []
                
                if not transgenic_ids: #this is the normal situation
                    if len(recs) > 2:
                        TPM_ranks = Counter()
                        for rec in recs:
                            TPM_ranks.update({rec.contig_name: rec.TPM})
                        two_most_common = [x[0] for x in TPM_ranks.most_common(2)]
                        for rec in recs:
                            if rec.contig_name not in two_most_common:
                                to_remove.append(rec)
                        
                else:
                    if len(recs) > 2:
                        TPM_ranks = Counter()
                        for rec in recs:
                            TPM_ranks.update({rec.contig_name: rec.TPM})
                        three_most_common = [x[0] for x in TPM_ranks.most_common(3)]
                        three_most_common_ids = []
                        for rec in recs:
                            if rec.contig_name in three_most_common:
                                three_most_common_ids.append(rec.identifier)
                        if len(set(transgenic_ids).intersection(set(three_most_common_ids))) > 0:
                            for rec in recs:
                                if rec.contig_name not in three_most_common:
                                    to_remove.append(rec)
                        else:
                            TPM_ranks = Counter()
                            for rec in recs:
                                TPM_ranks.update({rec.contig_name: rec.TPM})
                            two_most_common = [x[0] for x in TPM_ranks.most_common(2)]
                            for rec in recs:
                                if rec.contig_name not in two_most_common:
                                    to_remove.append(rec)
                
                for rec in to_remove:
                    self.all_recombinants[locus].remove(rec)
                
                
                
                #if len(recs) > 2:
                #    TPM_ranks = Counter()
                #    for rec in recs:
                #        TPM_ranks.update({rec.contig_name: rec.TPM})
                #    two_most_common = [x[0] for x in TPM_ranks.most_common(2)]
                #    to_remove = []
                #    for rec in recs:
                #        if rec.contig_name not in two_most_common:
                #            to_remove.append(rec)
                #    for rec in to_remove:
                #        self.all_recombinants[locus].remove(rec)

    def count_productive_recombinants(self, locus):
        recs = self.all_recombinants[locus]
        count = 0
        if recs is not None:
            for rec in recs:
                if rec.productive:
                    count += 1
        return (count)

    def count_total_recombinants(self, locus):
        recs = self.all_recombinants[locus]
        count = 0
        if recs is not None:
            count = len(recs)
        return (count)

    def get_trinity_lengths(self, locus):
        recs = self.all_recombinants[locus]
        lengths = []
        if recs is not None:
            for rec in recs:
                lengths.append(len(rec.trinity_seq))
        return (lengths)


class Recombinant(object):

    """Class to describe a recombined TCR locus as determined from the single-cell pipeline"""

    def __init__(self, contig_name, locus, identifier, all_poss_identifiers, productive, stop_codon, in_frame, TPM,
                 dna_seq, hit_table, summary, junction_details, best_VJ_names, alignment_summary, trinity_seq,
                 imgt_reconstructed_seq):
        self.contig_name = contig_name
        self.locus = locus
        self.identifier = identifier
        self.all_poss_identifiers = all_poss_identifiers
        self.productive = productive
        self.TPM = TPM
        self.dna_seq = dna_seq
        self.cdr3 = self._get_cdr3(dna_seq)
        self.hit_table = hit_table
        self.summary = summary
        self.junction_details = junction_details
        self.best_VJ_names = best_VJ_names
        self.alignment_summary = alignment_summary
        self.in_frame = in_frame
        self.stop_codon = stop_codon
        self.trinity_seq = trinity_seq
        self.imgt_reconstructed_seq = imgt_reconstructed_seq

    def __str__(self):
        return ("{} {} {} {}".format(self.identifier, self.productive, self.TPM))

    def _get_cdr3(self, dna_seq):
        aaseq = Seq(str(dna_seq), generic_dna).translate()
        if re.findall('FG.G', str(aaseq)) and re.findall('C', str(aaseq)):
            indices = [i for i, x in enumerate(aaseq) if x == 'C']
            upper = str(aaseq).find(re.findall('FG.G', str(aaseq))[0])
            for i in indices:
                if i < upper:
                    lower = i
            cdr3 = aaseq[lower:upper + 4]
        elif re.findall('FG.G', str(aaseq)):
            cdr3 = "Couldn't find conserved cysteine"
        elif re.findall('C', str(aaseq)):
            cdr3 = "Couldn't find FGXG"
        else:
            cdr3 = "Couldn't find either conserved boundary"
        return (cdr3)

    def get_summary(self):
        summary_string = "##{contig_name}##\n".format(contig_name=self.contig_name)
        if self.locus == 'A':
            V_segment = self.summary[0]
            J_segment = self.summary[1]
            segments_string = "V segment:\t{V_segment}\nJ segment:\t{J_segment}\n".format(V_segment=V_segment,
                                                                                          J_segment=J_segment)
        elif self.locus == 'B':
            V_segment = self.summary[0]
            D_segment = self.summary[1]
            J_segment = self.summary[2]
            segments_string = "V segment:\t{V_segment}\nD segment:\t{D_segment}\nJ segment:\t{J_segment}\n".format(
                V_segment=V_segment, D_segment=D_segment, J_segment=J_segment)
        summary_string += segments_string
        summary_string += "ID:\t{}\n".format(self.identifier)
        summary_string += "TPM:\t{TPM}\nProductive:\t{productive}\nStop codon:\t{stop_codon}\nIn frame:\t{in_frame}\n\n".format(
            TPM=self.TPM, productive=self.productive, stop_codon=self.stop_codon, in_frame=self.in_frame)

        summary_string += 'Segment\tquery_id\tsubject_id\t% identity\talignment length\tmismatches\tgap opens\tgaps\tq start\tq end\ts start\ts end\te value\tbit score\n'
        for line in self.hit_table:
            summary_string = summary_string + "\t".join(line) + "\n"
        return (summary_string)


