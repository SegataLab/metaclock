#!/usr/bin/env python
from __future__ import division
import re 
from multiprocessing import Pool
import numpy as np
import bz2
import sys
from Bio import SeqIO
from itertools import chain
import argparse

def pattern_filtering(p,l):
	"""
	pattern_filtering is a function taking a defined pattern and a list as arguments.
	It filter out elements in a list which does not match the given pattern 
	"""
	new_l = []
	for f in l:
		if re.match(p,f):
			new_l.append(f)
	return new_l		

def multi_map(w,f,l):
	"""
	multiprocessor is a function taking an integer(number of processor), a defined function,
	and a list of works as arguments. 
	"""
	pool = Pool(processes = w)
	return pool.map(f,l)



def unique(list_):
	"""
	It generates a list which only contains unique elements
	"""
	x = np.array(list_)
	return np.unique(x)

def find(s):
	"""
	It finds all indexes of gaps in a sequence
	"""

	return [i for i, ltr in enumerate(s) if ltr == '-']

def openr( fn, mode = 'r'):
	if fn is None:
		return sys.stdin
	return bz2.BZ2File(fn) if fn.endswith(".bz2") else open(fn, mode)


def openw( fn ):
	if fn is None:
		return sys.stdout
	return bz2.BZ2File(fn, 'w') if fn.endswith(".bz2") else open(fn, "w")

def is_numbers(s):
	try:
		int(s)
		return True
	except ValueError:
		return False

class stats(object):

	def __init__(self, sequence):
		self.sequence = sequence

	def calc_GapRatio(self):
		length=len(self.sequence)
		gaps = self.sequence.count('N') + self.sequence.count('-')
		r = round(gaps/length,2)
		return r

	def calc_genome(self):
		genome1 = self.sequence.replace("N","")
		genome2 = genome1.replace("-","")
		return len(self.sequence),len(genome2)
	
	def calc_CG(self):
		GC = self.sequence.count("C") +self.sequence.count("c")+self.sequence.count("g")+ self.sequence.count("G")
		if GC != 0:
		        return round(GC/len(self.sequence.replace("-","")),2)
		else:
		        return "0"

def out_stats(fna_file):
	seqio_dict = SeqIO.to_dict(SeqIO.parse(fna_file, "fasta"))
	opt = open("aln_stats.tsv","w")
	opt.write("SeqHeader\t"+ "GapRatio\t" + "AlignLength\t" + "SeqLength\t" +"GCcontent\n")
	for seq in seqio_dict:
		Seq = str(seqio_dict[seq].seq)
		S = stats(Seq)
		line = str(seq) + "\t" + str(S.calc_GapRatio())+ "\t" + str(S.calc_genome()[0])+ "\t" + str(S.calc_genome()[1]) + "\t" + str(S.calc_CG()) + '\n'
		opt.write(line)
	opt.close()	








