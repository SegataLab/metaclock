#!/usr/bin/env python
	
import subprocess
import argparse
import sys
from collections import defaultdict
from collections import Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio import AlignIO
from itertools import groupby
from operator import itemgetter
from BuildGenomeAln import make_blast_db
from BuildGenomeAln import blast_
from BuildGenomeAln import QC_on_blastn
from BuildGenomeAln import blastn_sort
from BuildGenomeAln import complement_seq
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment
from utils import multi_map
from utils import homo_site_mapper
from collections import ChainMap
import itertools

"""
NAME: HomoAlignsGenerator_dev_2.py

DESCRIPTION: HomoAlignsGenerator.py is a python program which samples alignments columns
			 based on coordinates of homolougous sites of references, and outputs homologous
			 columns into seperate files. Those sampled files are supposed to be homolougous 
			 alignments.

			 New feature: Solve the time-demanding issue by breaking down long sequences so as to
			 parallelize the whole work.  
"""

__author__ = "Kun D. Huang"
__version__ = "0.1"
__date__ = "19.02.2020"

def read_args(args):

	parser = argparse.ArgumentParser()
	parser.add_argument('alns_folder',
						nargs = "?",
						metavar = "alignments_folder",
						help = 'Specify the folder which contains alignments built using\
						 single reference. Warning: alignment file name should be consistent\
						 with reference genome name inside fasta',
						type = str)
	parser.add_argument('-l',
		                '--length',
		                help = 'minimum length of homologous sequence alignment to keep as a hit.',
		                type = str,
		                default = '500')
	parser.add_argument('-i',
						'--identity',
						help = 'minium identity of homologous sequences to keep as a hit.',
						 type = str,
						 default = '95.0')
	parser.add_argument('-p',
						'--nproc',
						help = 'Number of processors to use',
						default = 1,
						type = int)
	parser.add_argument('-hf',
						'--homo_site_in_refs',
						help = 'The number of references in which homologous sites can be found. defaul:[2]',
						default = 2.0,
						type = int)

	return vars(parser.parse_args())



def refs_collector(alns_folder):
	
	"""
	It takes a folder which contains single alignments and writes all reference
	genomes into one fasta file for blast procedure, meanwhile it returns a list
	of reference headers. 
	"""
	rec_lst = [] # a list of SeqRecords, each for one reference genome
	ref_header_lst = [] # a list of reference names
	all_alns = subprocess.getoutput('ls {}/*'.format(alns_folder)).split('\n')
	for aln in all_alns:
		RefSeq_name = aln.split('/')[1].split('-')[0]
		aln_dict = SeqIO.to_dict(SeqIO.parse(aln,'fasta'))
		rec_lst.append(SeqRecord(Seq(str(aln_dict[RefSeq_name].seq),\
		 generic_dna), id = RefSeq_name, description = ''))
		ref_header_lst.append(RefSeq_name)
	SeqIO.write(rec_lst, 'Refs.fna', 'fasta')
	return ref_header_lst



def partition_genomes(blastn_tab_cleaned):
	partitioned_dict = defaultdict(list)
	"""
	It partitions hits of genomes in pairwise style, excluding self comparison. It 
	returns a dictionary in which reference genome as key and matrix of its paired
	query as value. So whole matrix was partitioned into sub-matrices for each
	reference genome 
	"""
	ipt_mtx = (l.rstrip('\n').split('\t') for l in open(blastn_tab_cleaned).readlines())
	for l in ipt_mtx:
		if l[0] != l[1]:
			header = l[0] + '@' + l[1]
			partitioned_dict[header].append(l)
		else:
			continue
	return partitioned_dict

def SubWise_Refs(Refs_PairWise):
	"""
	It generates a dictionary in which subject reference name is the key
	and value is a list of query reference names. 
	"""
	RefSeqs_SubWise = defaultdict(list)
	for g in Refs_PairWise:
		RefSeqs_SubWise[g.split('@')[0]].append(g.split('@')[1])
	return RefSeqs_SubWise	


def merge(list_1, list_2, list_3, list_4, list_5, list_6):

	merged_list = tuple(zip(list_1, list_2, list_3, list_4, list_5, list_6))
	return merged_list


def find_coordinates_between_fragements(sub_matrix):

	"""
	It finds original postions of sseq and qseq for each pair of reference genomes
	"""

	s_pos_lst_genome = []
	q_pos_lst_genome = []
	s_site_lst_blast = []
	q_site_lst_blast = []
	s_flag_lst_blast = []
	q_flag_lst_blast = []
	for hit in sub_matrix:
		hit_object =blastn_sort(hit)
		hit_senario = [hit_object.sseq, hit_object.sstart, hit_object.send, hit_object.qseq\
		, hit_object.qstart, hit_object.qend]
		homo_site_mapper(hit_senario)
		q_pos = homo_site_mapper(hit_senario)[0]
		s_pos = homo_site_mapper(hit_senario)[1]
		q_site = homo_site_mapper(hit_senario)[2]
		s_site = homo_site_mapper(hit_senario)[3]
		q_flag = homo_site_mapper(hit_senario)[4]
		s_flag = homo_site_mapper(hit_senario)[5]

		s_pos_lst_genome.append(s_pos)
		q_pos_lst_genome.append(q_pos)
		s_site_lst_blast.append(s_site)
		q_site_lst_blast.append(q_site)
		s_flag_lst_blast.append(s_flag)
		q_flag_lst_blast.append(q_flag)

	return merge(q_pos_lst_genome, s_pos_lst_genome, q_site_lst_blast,\
	 s_site_lst_blast, q_flag_lst_blast, s_flag_lst_blast)


def reorder_aln(aln):
	
	"""
	It reorders order of taxa name in the alignment in alphabetical order. 
	"""

	align = MultipleSeqAlignment([], Gapped(IUPAC.unambiguous_dna, "-"))
	reorder = sorted(record.id for record in aln)
	
	for new_header in reorder:
		for old_aln_record in aln:
			r_id = str(old_aln_record.id)
			r_seq = str(old_aln_record.seq)
			if new_header == r_id:
				align.add_sequence(new_header, r_seq)
			else:
				continue

	return align


def col_finder_2(Refs_aln_dict, ref_name, pos, pos_lst):
	"""
	It finds column in an alignment giving the site coordinate

	Refs_aln_dict: a dictionary of alignments:, aln file name is the key and alignment is the value
	ref_name: subject sequence given for searching
	pos: the site position to smaple
	pos_lst: all positions in a list (for detecting the order of sequence)
	"""
	Ref_aln = Refs_aln_dict[ref_name]
	tax_number = len(Ref_aln[:, 2])

	pos_lst = list(filter(lambda a: a != '-', pos_lst))

	if pos_lst[-1] > pos_lst[1]:
		if pos != '-':
			return Ref_aln[:, pos]
		else:
			return '-'*tax_number
	else:
		if pos != '-':
			return complement_seq(Ref_aln[:, pos])
		else:
			return '-'*tax_number

def extract_column_flagged(aln, pos, flag):
	if flag == 'NO':
		col_seq = aln[:, pos]
	else:
		col_seq = complement_seq(aln[:, pos])
	return col_seq 	

def col_finder(Refs_aln_dict, que_homo_site, sub_homo_site):
	sub_aln = Refs_aln_dict[sub_homo_site.split('@')[0]]
	que_aln = Refs_aln_dict[que_homo_site.split('@')[0]]
	taxa_number = len(sub_aln[:, 2])

	sub_pos = sub_homo_site.split('@')[1]
	que_pos = que_homo_site.split('@')[1]
	s_flag = sub_homo_site.split('@')[3]
	q_flag = que_homo_site.split('@')[3]

	if (sub_pos != '-') & (que_pos != '-'):
		sub_homo_col = extract_column_flagged(sub_aln, int(sub_pos), s_flag)		
		que_homo_col = extract_column_flagged(que_aln, int(que_pos), q_flag)
	elif (sub_pos != '-') & (que_pos == '-'):
		que_homo_col = '-'*taxa_number
		sub_homo_col = extract_column_flagged(sub_aln, int(sub_pos), s_flag)
	elif (sub_pos == '-') & (que_pos != '-'):
		sub_homo_col = '-'*taxa_number
		que_homo_col = extract_column_flagged(que_aln, int(que_pos), q_flag)
	else:
		pass
	return que_homo_col, sub_homo_col


	
def unmerged(tup_lst):
	"""
	((a, b), (c, d)) ---> ((a, b), (c, d))
	"""
	list_1 = [i[0] for i in tup_lst]
	list_2 = [i[1] for i in tup_lst]

	return [list_1, list_2]




def de_duplicator(tup_lst):
	# it duplicates sequences with same flag in a tuple
	# Here use map() to speed up
	int_dict = {}
	for t in tup_lst:
		if t[0] not in int_dict:
			int_dict[t[0]]=t[1]
		else:
			continue
	return [int_dict[k] for k in int_dict]

def ColSampler(ref_sub, matrix):
	
	# ref_sub, matrix = multi_arg for nproc later
	"""
	Sampler samples columns from single alignments iteratively using homologous site
	positions in corresponding alignment and build final column-concatenated alignment.
	"""
	# single_aln_lst = [a for a in subprocess.getoutput('ls {}/*'.format(pars['alns_folder'])).split('\n')]

	aln_dict = {}	
	for aln in single_aln_lst:
		aln_dict[aln.split('/')[1]] = reorder_aln(AlignIO.read(aln, 'fasta'))
	
	"""
	Generate a dictionary, aln file name is the key and taxa-reordered alignment is the value.
	It is for sampling column from alignment file.
	"""



	sub_matrix=find_coordinates_between_fragements(matrix[ref_sub])
	
#$$$$# Here you can think about breaking all fragments into smaller fragments (500bp) labeled 
#$$$# with genome name and position so as to nporc

	# test = defaultdict(list)
	# ref_ori_dict = defaultdict(list)

	homo_site_list = []
	"""
	Find coordinates of subject seq and query seq for one sub-matrix
	"""
	for frg in sub_matrix:
		pos_pos = tuple(zip(frg[0], frg[1]))
		site_site = tuple(zip(frg[2], frg[3]))
		flag_flag = tuple(zip(frg[4], frg[5]))
		for n in range(len(pos_pos)):

			q_ref = ref_sub.split('@')[0] # extract query seq name
			s_ref = ref_sub.split('@')[1] # extract subject seq name

			que_flag = q_ref + '@' + str(pos_pos[n][0])+'@'+ site_site[n][0]+'@'+ flag_flag[n][0]# make a unique flag for sub seq [genome_name@position]
			sub_flag = s_ref +"@"+ str(pos_pos[n][1])+'@'+ site_site[n][1] +'@'+ flag_flag[n][1]# make a unique flag for query seq
			homo_site_list.append((que_flag, sub_flag))
	return homo_site_list


def nproc_sampler(processor, ref_sub_lst, matrix):

	matrix_lst = [matrix] * len(ref_sub_lst)
	return multi_map(processor, ColSampler, zip(ref_sub_lst, matrix_lst))


def check_max_uniq(site_lst):

	votes = {'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0}
	for s in site_lst:
		votes[s] += 1

	return max(votes, key=votes.get)


def consensus_col_builder(lst_homo_cols, refs_num = 2): # ? Check here, might be some problems.
	lst_homo_cols = [c.split("$")[-1] for c in lst_homo_cols]
	consensus_col = []
	if len(lst_homo_cols) >= refs_num: 
		for taxa in range(len(lst_homo_cols[0])):
			site_lst = []
			for col in range(len(lst_homo_cols)):
				site_lst.append(lst_homo_cols[col][taxa])
			# consensus_col.append(site_checker(site_lst))
			consensus_col.append(check_max_uniq(site_lst))
	
		return ''.join(consensus_col)
	else:
		pass

def aln_builder(aln, all_cols_list): # aln is the template alignment file	
	rec_lst = []
	reordered = sorted(record.id for record in aln)
	rotated = list(zip(*reversed(all_cols_list)))
	for taxa_idx in range(len(reordered)):
		seq = ''.join(rotated[taxa_idx])
		rec_lst.append(SeqRecord(Seq(seq, generic_dna), id = reordered[taxa_idx], description = ''))		
	return rec_lst
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# Last stage is to write a class which specifically handle
# the generator of all homologous column pairs.
# It should have features:
# 1. Create a dictionary list in which each query site is key
# and the list of subject sites is the value.
# 2. Generate different alignment files with minimum homologous col members (voting rule or specific cutoff)
# 3. Calculate size of shared sites by different Refs   
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

class column_toolkit(object):
	"""
	Object column_toolkit is to handle all homologous sites
	"""
	def __init__(self, homo_sites):
		self.homo_sites = homo_sites
	def merge_col_voting(self, number_refs = 2):
		merged_col_lst = [consensus_col_builder(self.homo_sites[c], number_refs) for c in self.homo_sites]
		res = list(filter(None, merged_col_lst))
		return res

if __name__== '__main__':


	pars = read_args(sys.argv)
	proc_log = open('processing_log.txt', 'w')

	single_aln_lst = [a for a in subprocess.getoutput('ls {}/*'.format(pars['alns_folder'])).split('\n')]
	aln_dict = {}	
	for aln in single_aln_lst:
		aln_dict[aln.split('/')[1]] = reorder_aln(AlignIO.read(aln, 'fasta'))	

	refs_name_headers = refs_collector(pars['alns_folder']) # refenrece genome names in a list
	ref_sub_groups = [i[0]+'@'+i[1] for i in itertools.combinations(refs_name_headers, 2)]
	proc_log.write("Write aligment RefSeq in each aligment into one fasta file ---> Refs.fna\n")
	make_blast_db('Refs.fna') # make blast database
	proc_log.write("Prepare blast database ---> [Refs.fna.nhr, Refs.fna.nin, Refs.fna.nsq]\n")
	blast_('Refs.fna', 'Refs.fna', 4) # call blast
	proc_log.write("Generate blast output ---> INTER_blast_opt_tmp.tab\n")
	QC_on_blastn('INTER_blast_opt_tmp.tab', pars['length'], pid = pars['identity']) # QC on blast result
	proc_log.write("Generate filtered blast output with length {} and identity {} ----> INTER_blastn_opt_tmp_cleaned.tab\n".\
		format(pars['length'], pars['identity']))
	aln_paths = [a for a in subprocess.getoutput('ls {}/*'.format(pars['alns_folder'])).split('\n')]

	matrice = {pair: partition_genomes('INTER_blastn_opt_tmp_cleaned.tab')[pair] for pair in ref_sub_groups}

	proc_log.write("-----------------------------------------\nToatl length homologous fragments\n\
	# Genome pairs\t Length\n")

	sites_containner = defaultdict(set)
	for pair in ref_sub_groups:
		sub_matric_frag_len = find_coordinates_between_fragements(matrice[pair])
		tot_len = sum([len(i[0]) for i in sub_matric_frag_len])
		proc_log.write("{}\t{}nt\n".format(pair.replace('@','X'), tot_len))		
		col_sampled = ColSampler(pair, matrice)
		for c in col_sampled:
			sites = col_finder(aln_dict, c[0], c[1])
			que_tag = c[0]
			que_col = c[0]+'$'+sites[0]
			sub_col = c[1]+'$'+sites[1]
			sites_containner[que_tag].add(que_col)
			sites_containner[que_tag].add(sub_col)

	col_tool_obj = column_toolkit(sites_containner)
	cols_merged = col_tool_obj.merge_col_voting(pars['homo_site_in_refs'])
	merged = aln_builder(AlignIO.read('{}/{}'.format(pars['alns_folder'], refs_name_headers[0]), 'fasta'),\
	 cols_merged)
	opt_merged_col_name = 'merged_aln_{}_refs.fna'.format(str(pars['homo_site_in_refs']))
	SeqIO.write(merged, opt_merged_col_name, 'fasta')

# 1. Add a feature of calculating mutation/site

#------------------------------------------------------------------------------
	#Distribute all columns into different files

	# for i in range(len(aln_paths)):
	# 	lst = [] # list of columns at one position  	
	# 	for j in all_cols:
	# 		lst.append(list(j[i]))

	# 	final = aln_builder(AlignIO.read('{}/{}'.format(pars['alns_folder'], refs_name_headers[0]), 'fasta'), lst)
	# 	SeqIO.write(final,'HomoAln_{}.fna'.format(str(i)), 'fasta')
#-------------------------------------------------------------------------------


# ------------------------------------------------------------------------------	

# Notes: 
#1. sampled coulmn redundancy can be solved by applying harder quality control on selection blastn results
#2. be careful, reference alignment file name should be consistent with reference genome name inside fasta
