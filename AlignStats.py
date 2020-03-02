#!/usr/bin/env python

import itertools
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped
from Bio import AlignIO
import sys
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import SeqIO


class AlignStats(object):

	def __init__(self, Align_opt):
		self.Align_opt = Align_opt

	def variant_sites_ratio(self):
		multi_count = 0
		bi_count = 0
		for j in range(len(self.Align_opt[1])):
			col = set([i for i in list(self.Align_opt[:,j]) if i != '-'])
			if len(col) == 2:
				bi_count += 1
			elif len(col) > 2:
				multi_count += 1
		all_count = bi_count + multi_count
		bi_ratio = bi_count/len(self.Align_opt[1])
		multi_ratio = all_count/len(self.Align_opt[1])
		return (bi_ratio, multi_ratio, bi_count, all_count)

	def pairwise_distance(self):
		# seqs_lst = [self.Align_opt[i].id for i in range(len(self.Align_opt[:,1]))]
		# pair_wise_seqs = itertools.combinations(seqs_lst, 2)
		sub_aln = self.Align_opt[:,:]
		calculator = DistanceCalculator('identity')
		dm = calculator.get_distance(sub_aln)
		return dm

	def remove_gappy_taxa(self, gap = 0.1):

		filtered_taxa = []
		for record in self.Align_opt:
			gap_ratio = record.seq.count('-')/len(record)
			if gap_ratio <= gap:
				filtered_taxa.append(record)
		return filtered_taxa		
			

			
		
