#!/usr/bin/env python

"""
NAME: CoreMaker.py

Description: Make_core_on_aDNA.py is a python script to sample columns which are stricly core based on
certain genome you provided (They usually are ancient genomes).
  
"""

__author__="Kun D. Huang"
__version__="0.2"
__date__="12.04.2019"

import argparse
import sys
from Bio import SeqIO
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from operator import itemgetter



def read_args(args):

	parser = argparse.ArgumentParser()

	parser.add_argument('ipt_lst',
						nargs = "?",
						metavar = "input_list",
						help = "Input a file which has a list of genomes on which you want to\
						build core genome",
						type = str)

	parser.add_argument('ipt_aln',
						nargs = "?",
						metavar = "input_alignment",
						help = "Input the genome alignment",
						type = str)

	parser.add_argument('opt_aln',
						nargs = "?",
						metavar = "output_alignment",
						help = "Output the genome alignment optimized",
						type = str)
	parser.add_argument('gt',
						nargs = '?',
						metavar = 'gap_threshold',
						help = "Gap ratio allowed in columns of provided genomes (aDNA) for stricting coreness.",
						type = float)

	return vars(parser.parse_args())

def MakeCore(Seq_dict, Seq_lst, opt_file, gt):
	
	core_col = []
	Seq_lst_lst = [list(Seq_dict[g].seq) for g in Seq_lst]

	for c in range(len(list(Seq_dict[Seq_lst[0]].seq))):
		col_sites = [i[c] for i in Seq_lst_lst]
		tot_sites_num = len(col_sites)
		if col_sites.count('-')/tot_sites_num > gt+0.00001:
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

	SeqIO.write(rec_list, opt_file, 'fasta')


if __name__ == '__main__':
	pars = read_args(sys.argv)
	Seq_dict = SeqIO.to_dict(SeqIO.parse(pars['ipt_aln'], "fasta"))
	ipt_genomes = [i.rstrip() for i in open(pars['ipt_lst']).readlines()]
	MakeCore(Seq_dict, ipt_genomes, pars['opt_aln'], pars['gt'])




