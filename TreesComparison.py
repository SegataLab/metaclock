#!/usr/bin/env python
import subprocess
import os
import sys
from ete3 import Tree
from itertools import combinations

all_alns_folder = sys.argv[1]

paths =[os.getcwd()+'/'+ i for i in subprocess.getoutput('ls {}/*'.format(all_alns_folder)).split('\n')]

group = combinations(paths,2)

def tree_reader(tree_file):
	t = Tree(open(tree_file).read().rstrip(), format = 5)
	return t

# t1 = tree_reader("aln_2_case_1.tre")
# t2 = tree_reader("aln_2_case_2.tre")	

def calc_rf(t1, t2):

	result = t1.robinson_foulds(t2, unrooted_trees=True)

	normalized_rf = int(result[0])/int(result[1])

	return normalized_rf


for i in group:
	t1 = tree_reader(i[0])
	t2 = tree_reader(i[1])
	print("1_Refs"+"\t"+str(calc_rf(t1, t2)))