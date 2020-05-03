#!/usr/bin/env python

import sys
import math
import subprocess
import argparse
import os
import itertools
import shutil
import utils
from collections import defaultdict
from functools import partial
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, IUPAC, Gapped
from Bio.Seq import Seq
from AlignStats import AlignStats
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import pandas as pd
from operator import itemgetter


def read_args(args):

	parser = argparse.ArgumentParser()

	parser.add_argument('ga_file',
						nargs = "?",
						metavar = "genome_alignment_file",
						help = "Input a genome alignment needs to be tailored. [fasta] format",
						type = str)

	parser.add_argument('tailored_ga_file',
						nargs = "?",
						metavar = "tailored_genome_alignment_file",
						help = "Output a tailored genome alignment in fasta format",
						type = str)

	parser.add_argument('-a',
						"--automated_tailoring",
						help = "This option allows automated tailoring for minimizing missing information in the alignment",
						action = 'store_true')
	parser.add_argument('-sa',
						'--semi_automated_tailoring',
						help = 'This option allows semi-automated tailoring with customized quantile for gap score.\
						 e.g. -sa 0.75,0.6 [Modern samples with gap score outside uppper 0.75 quantile are out, \
						 ancient samples with gap score outside upper 0.60 quantile are out.]',
						type = str)
	parser.add_argument('-sl',
						"--short_list",
						help = "Input a list of taxa to keep for alignment tailoring. Optionally, the limit for missing info in each column can\
						 be given and it is delimited by comma, or column tailoring will be performed automatically. e.g. list.txt,0.1",
						type = str)
	parser.add_argument('-m',
						'--manual_tailoring',
						help = "Missing info for modern, ancient samples and each column is manually.\
						e.g. 0.1,0.2,0.3 [Modern taxa with missing info > 0.1 are removed, and ancient\
						 taxa with missing info >0.2 are removed. Afterwards, columns with missing info > 0.3 are removed.]",
						type = str)

	return vars(parser.parse_args())


def tailoring_cols(Seq_dict, c):
	
	core_col = []
	Seq_lst_lst = [list(Seq_dict[g].seq) for g in Seq_dict]
	Seq_lst = [key for key in Seq_dict]

	for c in range(len(list(Seq_dict[Seq_lst[0]].seq))):
		col_sites = [i[c] for i in Seq_lst_lst]
		tot_sites_num = len(col_sites)
		if col_sites.count('-')/tot_sites_num > c+0.00001:
			continue
		else:
			core_col.append(c)

	rec_list=[]
	for g in Seq_dict:
		seq_lst = list(Seq_dict[g].seq)
		core_seq = itemgetter(*core_col)(seq_lst)
		seq = "".join(core_seq)
		_id = g
		rec_list.append(SeqRecord(Seq(seq, generic_dna), id = _id, description = ''))

	return rec_list

class tailor(object):

	def __init__(self, aln):
		self.aln = aln
	def manual_tailor(self, m, a, c):
		#m: gap ratio for modern samples.
		#a: gap ratio for acnient samples.
		#c: gap ratio for columns after gappy taxa being removed.
		selected_aln = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
		for g in self.aln:
			g_stat = utils.stats(g.seq)
			ID = str(g.id)
			seq = str(g.seq)
			gap = g_stat.calc_GapRatio()
			if g.id.startswith('a__') and gap <= a:
				selected_aln.add_sequence(ID, seq)
			elif g.id.startswith('m__') and gap <= m:
				selected_aln.add_sequence(ID, seq)
			elif (not g.id.startswith('a__')) and (not g.id.startswith('m__')):
				selected_aln.add_sequence(ID, seq)

		selected_aln_dict = SeqIO.to_dict(selected_aln)
		return tailoring_cols(selected_aln_dict, c)

	def short_list_tailor(self, short_list, c):
		# c is the cutoff for column
		# If c is not specified Trimal automated trimming will be activated
		Seq_dict = SeqIO.to_dict(self.aln)
		selected_aln = []
		for i in Seq_dict:
			if i in short_list:
				selected_aln.append(SeqRecord(Seq(str(Seq_dict[i].seq), generic_dna), i, description = ''))
			else:
				continue
		if c:
			selected_aln_dict = SeqIO.to_dict(selected_aln)
			return tailoring_cols(selected_aln_dict, c)
		else:
			# Here trimal gappyout is implemented. 
			# Make sure trimal is set in your system
			SeqIO.write(selected_aln, 'tmp_feed_to_trimal.fna', 'fasta')
			cmd = 'trimal -gappyout -in tmp_feed_to_trimal.fna -out TrimalGappyout_opt.fna'
			subprocess.call(cmd, shell = True)

	def auto_tailoring(self, m_q, a_q):
		# For rows, gap score > 95% quantile will be removed
		# For columns, Trimal gappyout is used
		sample, gap_ratio, type_ = [], [], []
		for g in self.aln:
			g_stat = utils.stats(g.seq)
			gap_score = g_stat.calc_GapRatio()
			if g.id.startswith('a__'): 
				sample.append(g.id) 
				gap_ratio.append(gap_score)
				type_.append('Ancient')
			elif g.id.startswith('m__'):
				sample.append(g.id)
				gap_ratio.append(gap_score)
				type_.append('Modern')
			elif (not g.id.startswith('a__')) and (not g.id.startswith('m__')):
				sample.append(g.id)
				gap_ratio.append(gap_score)
				type_.append('RefSeq')
		df_ = pd.DataFrame.from_dict({'Sample': sample, 'Gap score': gap_ratio, 'Type': type_})
		
		m_df_q = df_[df_['Type'] == 'Modern']['Gap score'].quantile(m_q)
		selected_m = df_[(df_['Type'] == 'Modern') & (df_['Gap score'] < m_df_q)]
		s_m_samples = list(selected_m['Sample'])
		a_df_q = df_[df_['Type'] == 'Ancient']['Gap score'].quantile(a_q)
		selected_a = df_[(df_['Type'] == 'Ancient') & (df_['Gap score'] < a_df_q)]
		s_a_samples = list(selected_a['Sample'])
		RefSeq = [i for i in list(df_['Sample']) if (not i.startswith('a__')) and (not i.startswith('m__'))]
		selected_samples = s_m_samples + s_a_samples + RefSeq

		Seq_dict = SeqIO.to_dict(self.aln)
		select_taxa = []
		for s in selected_samples:
			record = Seq_dict[s]
			header = str(record.id)
			seq = str(record.seq)
			newR = SeqRecord(Seq(seq, generic_dna), id = str(header), description = '')
			select_taxa.append(newR)
		
		SeqIO.write(select_taxa, 'tmp_feed_to_trimal.fna', 'fasta')
		cmd = 'trimal -gappyout -in tmp_feed_to_trimal.fna -out TrimalGappyout_opt.fna'
		subprocess.call(cmd, shell = True)

		




def manual_tailoring_pars_check(par):
	par_split = par.split(',')
	if len(par_split) == 3:
		return (float(par_split[0]), float(par_split[1]), float(par_split[2]))
	else:
		sys.exit('Please input gap ratio limit for rows of modern and ancient taxa and columns. e.g. 0.1,0.2,0.3')

def short_list_tailaring_pars_check(par):
	if ',' in par:
		par_split = [p for p in par.split(',') if len(p) > 0]
		if len(par_split) == 2:
			return (par_split[0], float(par_split[1]))
		else:
			sys.exit('Please specify short list file and gap cutoff for each column. e.g. text.txt,0.1')
	else:
		return (par, None)

def semi_auto_tailoring_pars_check(par):

	if ',' in par:
		par_split = [p for p in par.split(',') if len(p) > 0]
		if len(par_split) == 2:
			return (float(par_split[0]), float(par_split[1]))
		else:
			sys.exit('Please specify upper quantile of gap score for modern and ancient samples. e.g. 0.75,0.75')
	else:
		sys.exit('Please specify upper quantile of gap score for modern and ancient samples. e.g. 0.75,0.75')

				


if __name__ == '__main__':
	pars = read_args(sys.argv)
	wga_aln = AlignIO.read(pars['ga_file'], 'fasta')
	aln_obj = tailor(wga_aln)
	if pars['manual_tailoring']:
		m_size = manual_tailoring_pars_check(pars['manual_tailoring'])[0]
		a_size = manual_tailoring_pars_check(pars['manual_tailoring'])[1]
		c_size = manual_tailoring_pars_check(pars['manual_tailoring'])[2]
		SeqIO.write(aln_obj.manual_tailor(m_size, a_size, c_size), pars['tailored_ga_file'] ,'fasta')

	elif pars['short_list']:

		sl_args = short_list_tailaring_pars_check(pars['short_list'])
		short_list = [i.rstrip() for i in open(sl_args[0]).readlines()]
		c = sl_args[1]

		if c:
			SeqIO.write(aln_obj.short_list_tailor(short_list, c), pars['tailored_ga_file'] ,'fasta')
		else:
			aln_obj.short_list_tailor(short_list, c)
			subprocess.call('mv TrimalGappyout_opt.fna {}'.format(pars['tailored_ga_file']), shell = True)
			subprocess.call('rm tmp_feed_to_trimal.fna', shell = True)
	elif pars['automated_tailoring']:
		aln_obj.auto_tailoring(0.75, 0.6)
		subprocess.call('mv TrimalGappyout_opt.fna {}'.format(pars['tailored_ga_file']), shell = True)
		subprocess.call('rm tmp_feed_to_trimal.fna', shell = True)

	elif pars['semi_automated_tailoring']:
		check_par = semi_auto_tailoring_pars_check(pars['semi_automated_tailoring'])
		m_q = check_par[0]
		a_q = check_par[1]
		aln_obj.auto_tailoring(m_q, a_q)
		subprocess.call('mv TrimalGappyout_opt.fna {}'.format(pars['tailored_ga_file']), shell = True)
		subprocess.call('rm tmp_feed_to_trimal.fna', shell = True)
	else:
		sys.exit('Please choose one of these tailor strategies: -a, -sa, -sl, -m')

	








