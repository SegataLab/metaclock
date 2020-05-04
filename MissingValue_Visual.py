#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import itertools
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import seaborn as sns


def read_args(args):

	parser = argparse.ArgumentParser()

	parser.add_argument('ipt_aln',
						nargs = "?",
						metavar = "input_alignment",
						help = "Input a MSA alignment in fasta",
						type = str)

	parser.add_argument('opt_name',
						nargs = "?",
						metavar = "output_name",
						help = "Specify the output name for PNG figure",
						type = str)
	parser.add_argument('-w',
						'--window_size',
						nargs = '?',
						help = "Sliding window size. default:[200bp]",
						type = int,
						default = 200)
	parser.add_argument('-s',
						'--step_size',
						nargs = '?',
						help = 'Moving step of sliding window. default:[50bp]',
						type = int,
						default = 50)

	return vars(parser.parse_args())




def GenomeBar(genome_dict, opt_png, f_name = None, min_v = 0, max_v = 1):

	fig = plt.figure(figsize=(15,20))	
	x_scale = np.linspace(0,1, 100)


	axprops = dict(xticks=[], yticks=[]) 
	barprops = dict(aspect='auto', cmap=plt.cm.Blues, interpolation='nearest',vmin = min_v, vmax = max_v)

	bot = 0.8
	plt.axes([0.85, 0.85, 0.1, 0.000001])
	plt.yticks([])
	plt.xticks([0, 0.25, 0.5, 0.75, 1])

	ax_scale = fig.add_axes([0.85, 0.85, 0.1, 0.01], **axprops)
	ax_scale.imshow(x_scale.reshape((1, -1)), **barprops)
	ax_scale.text(0.5, 1.2, 'Reconstruction score', ha = 'center', transform=ax_scale.transAxes)

	for g in genome_dict:
		header = g
		x = np.asarray(genome_dict[g])
		ax = fig.add_axes([0.25, bot, 0.7, 0.01], **axprops)
		ax.imshow(x.reshape(1,-1), **barprops)
		ax.text(-0.23, 0.35, header, ha = 'left', transform=ax.transAxes)
		bot -= 0.01

	return plt.savefig(opt_png)

def color_scale(sites_lst=range(20)):
	x = np.asarray(sites_lst)
	axprops = dict(xticks=[], yticks=[]) 
	barprops = dict(aspect='auto', cmap=plt.cm.seismic, interpolation='nearest',vmin = 0, vmax = 20)
	fig=plt.figure(figsize=(8,3))
	plt.axes([0.05,0.1,0.1,0.00001])
	plt.yticks([])
	plt.xticks(np.arange(0, 25, step=5))

	ax2 = fig.add_axes([0.05, 0.15, 0.1, 0.07], **axprops) # Draw a rectangle 
	ax2.imshow(x.reshape((1,-1)), **barprops)

	
	return plt.savefig('color_scale.png')
 

def win_counting_mv(g_seq, win_size, step_size):
	mv_lst = []
	for i in range(0, len(g_seq), step_size):
		stretch_mv = g_seq[i: i+win_size].count('-')
		mv_ratio = 1- (stretch_mv/win_size)
		mv_lst.append(mv_ratio)
	return mv_lst 

		
def window_sliding_mv(fna, win_size, step_size):
	# This function calculates missing value along genome in sliding window manner
	# and return a dictionary in which header is key and the list of windowed values is value
	Seq_dict = SeqIO.to_dict(SeqIO.parse(fna, 'fasta'))
	mv_dict = {}
	for i in Seq_dict:
		seq = Seq_dict[i].seq
		mv_dict[i] = win_counting_mv(seq, win_size, step_size)
	return mv_dict

class Aln_Gap_Score(object):

	def __init__(self, aln):
		self.aln = aln
	def samples_gap_score(self):
		Samples, gap_score, type_ = [], [], []
		for g in self.aln:
			ID = g.id
			gap_ratio = g.seq.count('-')/len(g.seq)
			if ID.startswith('a__'):
				Samples.append(ID)
				gap_score.append((1 - gap_ratio)*100)
				type_.append('Ancient')
			else:
				Samples.append(ID)
				gap_score.append((1 - gap_ratio)*100)
				type_.append('Modern')
		df_ = pd.DataFrame.from_dict({'Samples': Samples, 'Reconstruction score': gap_score, 'Type': type_})
		
		return df_.sort_values('Reconstruction score', ascending = False)




if __name__ == '__main__':
	pars = read_args(sys.argv)
	windowed_mv = window_sliding_mv(pars['ipt_aln'], pars['window_size'], pars['step_size'])
	GenomeBar(windowed_mv, pars['opt_name'])
	wga = AlignIO.read(pars['ipt_aln'], 'fasta')
	aln_obj = Aln_Gap_Score(wga)
	df = aln_obj.samples_gap_score()
	fig, ax = plt.subplots(figsize=(15,15))
	ax = sns.barplot(x="Reconstruction score", y="Samples", hue = 'Type', data=df,\
	 dodge=False, order = df['Samples'])
	ax.set_yticklabels(ax.get_yticklabels(), size = 10)
	ax.set_xticklabels(ax.get_xticks(), size = 15)	
	ax.set_xlabel('Genome reconstructed (%)', fontsize = 15)
	ax.set_ylabel('Samples', fontsize = 15)
	fig.savefig('Sample_gap_score.png', bbox_inches="tight")


