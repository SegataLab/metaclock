#!/usr/bin/env python

import itertools
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped
from Bio import AlignIO
import sys
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import SeqIO
import pandas as pd

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

	def column_MissingValue_dist(self):
		misinfo = []
		for j in range(len(self.Align_opt[1])):
			col = list(self.Align_opt[:,j])
			misvalues = col.count('-')/len(col)
			misinfo.append(misvalues)
		return misinfo
	def column_MissingValue_bar(self):

		cutoffs_all = {}
		cutoffs_a = {}
		cutoffs_m = {}

		aln_len = self.Align_opt.get_alignment_length()
		aln_a = MultipleSeqAlignment([i for i in self.Align_opt if i.id.startswith('a__')])
		aln_m = MultipleSeqAlignment([i for i in self.Align_opt if not i.id.startswith('a__')])
		
		all_a_m_cols = []
		for j in range(len(self.Align_opt[1])):
			col = list(self.Align_opt[:,j])
			col_mis = col.count('-')/len(col)
			col_a = list(aln_a[:,j])
			col_a_mis = col_a.count('-')/len(col_a)
			col_m = list(aln_m[:,j])
			col_m_mis = col_m.count('-')/len(col_m)
			all_a_m_cols.append((col_mis, col_a_mis, col_m_mis))

		misv_all_col = [i[0] for i in all_a_m_cols]
		misv_a_col = [i[1] for i in all_a_m_cols]
		misv_m_col = [i[2] for i in all_a_m_cols]


		cutoffs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
		
		for c in cutoffs:
			cutoffs_all[c] = sum(i < c for i in misv_all_col)
			cutoffs_a[c] = sum(i < c for i in misv_a_col)
			cutoffs_m[c] = sum(i < c for i in misv_m_col)
		data = {'Samples': ['AllSamples', 'AncientSamples', 'ModernSamples']}
		for i in cutoffs_all:
			data[str(i)] = [cutoffs_all[i], cutoffs_a[i], cutoffs_m[i]]
		df = pd.DataFrame(data, columns = ['Samples', '0.1', '0.2', '0.3', '0.4', '0.5'\
			, '0.6', '0.7', '0.8', '0.9', '1'])

		return df.melt(id_vars = 'Samples', var_name = 'Cutoffs', value_name = 'N. columns')










		

			

			
		
