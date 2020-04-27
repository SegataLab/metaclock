#!/usr/bin/env python 

import sys
import subprocess
from itertools import combinations
import os


all_alns_folder = sys.argv[1]
number_of_members = sys.argv[2]

paths =[os.getcwd()+'/'+ i for i in subprocess.getoutput('ls {}/*'.format(all_alns_folder)).split('\n')]

group = combinations(paths,int(number_of_members))

counter = 0
for i in group:
	counter += 1
	cwd = os.getcwd()
	# os.mkdir(cwd+ '/aln_'+number_of_members+"_case_"+str(counter))
	folder_name = cwd+ '/aln_'+number_of_members+"_case_"+str(counter)
	os.mkdir(folder_name)
	for a in i:
		subprocess.call("ln -s {} {}/.".format(a, folder_name), shell = True)
		print("moving {} to {}".format(a, folder_name))

