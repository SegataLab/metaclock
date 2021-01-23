#!/usr/bin/env python

import subprocess
import sys
import shutil
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq


"""
This script is to automatically prepare gffs for Multiple RefSeqs merging mode.
"""

alns_folder = sys.argv[1]

all_alns_fna = ['alns/'+ i for i in subprocess.getoutput('ls {}'.format(alns_folder)).split('\n')]

def label_length_check(label):
	if len(label) >= 20:
		short_lab = 'Seq_00001'
	else:
		short_lab = label

	return short_lab

def extract_refseq(alns):
	for aln in alns:
		rec_lst = []
		label = aln.split('/')[-1]
		aln_dict = SeqIO.to_dict(SeqIO.parse(aln,'fasta'))
		aln_id = label_length_check(aln_dict[label].id)
		aln_seq = aln_dict[label].seq
		rec_lst.append(SeqRecord(Seq(str(aln_seq),\
		 generic_dna), id = aln_id, description = ''))
		SeqIO.write(rec_lst, label, 'fasta')

def run_prokka(alns):
	all_refseqs = []
	all_outdirs = []
	all_prefix = []
	for aln in alns:
		refseq = aln.split('/')[-1]
		all_refseqs.append(refseq)
		outdir = refseq.split('.')[0]
		all_outdirs.append(outdir)
		prefix = refseq

		cmd = "prokka --outdir {} --prefix {} {}".format(outdir, prefix, refseq)
		subprocess.call(cmd, shell = True)

	return all_outdirs, all_refseqs 


opt_dir = os.getcwd() + '/{}'.format('gffs_folder')
if os.path.exists(opt_dir):
	shutil.rmtree(opt_dir)
os.makedirs(opt_dir)


extract_refseq(all_alns_fna)
tmp_itm = run_prokka(all_alns_fna)

for i in range(len(tmp_itm[0])):
	outdir = tmp_itm[0][i]
	gff = outdir+'/'+ tmp_itm[1][i] + '.gff'
	cmd = 'cp {} {}/.'.format(gff, opt_dir)
	subprocess.call(cmd, shell = True)


for i in range(len(tmp_itm[0])):
	outdir = tmp_itm[0][i]
	prefix = tmp_itm[1][i]
	cmd_1 = 'rm -r {}'.format(outdir)
	cmd_2 = 'rm -r {}'.format(prefix)
	subprocess.call(cmd_1, shell = True)
	subprocess.call(cmd_2, shell = True)


