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




def GenomeBar(genome_dict, f_name = None, min_v = 0, max_v = 1):

	fig = plt.figure(figsize=(15,20))	
	x_scale = np.linspace(0,1, 100)


	axprops = dict(xticks=[], yticks=[]) 
	barprops = dict(aspect='auto', cmap=plt.cm.seismic,interpolation='nearest',vmin = min_v, vmax = max_v)

	bot = 0.8
	plt.axes([0.65, 0.85, 0.1, 0.000001])
	plt.yticks([])
	plt.xticks([0, 0.25, 0.5, 0.75, 1])

	ax_scale = fig.add_axes([0.65, 0.85, 0.1, 0.01], **axprops)
	ax_scale.imshow(x_scale.reshape((1, -1)), **barprops)

	for g in genome_dict:
		header = g
		x = np.asarray(genome_dict[g])
		ax = fig.add_axes([0.25, bot, 0.5, 0.02], **axprops)
		ax.imshow(x.reshape(1,-1), **barprops)
		ax.text(-0.23, 0.4, header, transform=ax.transAxes)
		bot -= 0.02

	
	return plt.savefig('test.png')

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



if __name__ == '__main__':
	pars = read_args(sys.argv)
	windowed_mv = window_sliding_mv(pars['ipt_aln'], pars['window_size'], pars['step_size'])
	# one_sample = windowed_mv['a__calc_2090']
	# genome_lst = [windowed_mv[i] for i in windowed_mv]
	GenomeBar(windowed_mv)
	# color_scale(range(20))


