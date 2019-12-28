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
from collections import ChainMap

"""
NAME: HomoAlignsGenerator.py

DESCRIPTION: HomoAlignsGenerator.py is a python program which samples alignments columns
			 based on coordinates of homolougous sites of references, and outputs homologous
			 columns into seperate files. Those sampled files are supposed to be homolougous 
			 alignments. 
"""

__author__ = "Kun D. Huang"
__version__ = "0.1"
__date__ = "17.06.2019"

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
	parser.add_argument('-c',
						'--coreness',
						help = 'The coreness of references',
						default = 1.0,
						type = float)

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
			new_header = l[1]+'@'+l[0]
			partitioned_dict[new_header].append(l)
		else:
			continue

	return partitioned_dict

def SubWise_Refs(Refs_PairWise):
	"""
	It generates a dictionary in which subject reference is the key
	and value is a list of query references. 
	"""
	RefSeqs_SubWise = defaultdict(list)
	for g in Refs_PairWise:
		RefSeqs_SubWise[g.split('@')[0]].append(g.split('@')[1])
	return RefSeqs_SubWise	


def sub_que_matcher(sseq, qseq, s_idx, q_idx):
	"""
	It takes subject sequence and query sequence and their
	alignment start positions, and return the postion of homologous
	sites on genome before alignment.
	"""
	s_pos_lst = []
	q_pos_lst = []
	if len(sseq.replace('-','')) != len(qseq.replace('-','')):
		for n in range(1,len(sseq)+1):
			idx = n - 1
			if sseq[idx] != '-':
				s_pos_lst.append(s_idx)
				s_idx += 1
				if qseq[idx] != '-':
					q_pos_lst.append(q_idx)
					q_idx += 1
				else:
					q_pos_lst.append('-')
	else:
		for n in range(len(sseq)):
			s_pos_lst.append(s_idx)
			s_idx += 1
			q_pos_lst.append(q_idx)
			q_idx += 1

	return s_pos_lst, q_pos_lst			



def homo_sites_finder(sstart, send, sseq, qstart, qend, qseq):
	"""
	It finds paired site postions of homologous sequences 
	"""
	s_poses = []
	q_poses = []

	if ((send - sstart) > 0) & ((qend - qstart) > 0):
		# subject and query sequences are both --> direction
		s_idx = sstart - 1
		q_idx = qstart - 1
		s_poses.extend(sub_que_matcher(sseq, qseq, s_idx, q_idx)[0])
		q_poses.extend(sub_que_matcher(sseq, qseq, s_idx, q_idx)[1])

	elif ((send - sstart) > 0) & ((qend - qstart) < 0):
		# subject follows --> direction, and query follows <-- direction
		s_idx = sstart - 1
		q_idx = qend - 1 
		s_poses.extend(sub_que_matcher(sseq, qseq, s_idx, q_idx)[0])
		q_poses.extend(sub_que_matcher(sseq, qseq, s_idx, q_idx)[1])

	
	elif ((send - sstart) < 0) & ((qend - qstart) > 0):
		# subject follows <-- direction, query follows --> direction
		s_idx = send - 1 
		q_idx = qstart - 1
		s_poses.extend(sub_que_matcher(sseq, qseq, s_idx, q_idx)[0][::-1])
		q_poses.extend(sub_que_matcher(sseq, qseq, s_idx, q_idx)[1])


	else:
		# subject follows <-- direction, query follows <-- direction
		s_idx = send - 1
		q_idx = qend - 1
		s_poses.extend(sub_que_matcher(sseq, qseq, s_idx, q_idx)[0][::-1])
		q_poses.extend(sub_que_matcher(sseq, qseq, s_idx, q_idx)[1][::-1])

	return s_poses, q_poses

def merge(list_1, list_2):
	"""
	([1, 2, 3],[a, b, c]) ----> ((1, a), (2, b), (3, c))
	"""

	merged_list = tuple(zip(list_1, list_2))
	return merged_list

def find_coordinates_between_fragements(sub_matrix):

	"""
	It finds original postions of sseq and qseq for each pair of reference genomes
	"""

	s_pos_lst_genome = []
	q_pos_lst_genome = []
	for hit in sub_matrix:
		hit_object =blastn_sort(hit)
		s_pos = homo_sites_finder(hit_object.sstart, hit_object.send, hit_object.sseq,\
		 hit_object.qstart, hit_object.qend, hit_object.qseq)[0]
		q_pos = homo_sites_finder(hit_object.sstart, hit_object.send, hit_object.sseq,\
		 hit_object.qstart, hit_object.qend, hit_object.qseq)[1]

		s_pos_lst_genome.append(s_pos)
		q_pos_lst_genome.append(q_pos)

	return merge(s_pos_lst_genome, q_pos_lst_genome)


		

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

def ColSampler(ref_sub):
	
	"""
	Sampler samples columns from single alignments iteratively using homologous site
	positions in corresponding alignment and build final column-concatenated alignment.
	"""

	matrix = partition_genomes('INTER_blastn_opt_tmp_cleaned.tab') # Partition the whole blast table into sub-matrix, based on ref seq name
	single_aln_lst = [a for a in subprocess.getoutput('ls {}/*'.format(pars['alns_folder'])).split('\n')]

	aln_dict = {}	
	for aln in single_aln_lst:
		aln_dict[aln.split('/')[1]] = reorder_aln(AlignIO.read(aln, 'fasta'))
	
	"""
	Generate a dictionary, aln file is the key and taxa-reordered alignment is the value.
	"""


	ref_ori_dict = defaultdict(list)

	test = defaultdict(list)

	sub_matrix=find_coordinates_between_fragements(matrix[ref_sub])
	"""
	Find coordinates of subject seq and query seq for one sub-matrix
	"""
	for frg in sub_matrix:
		pos_pos = merge(frg[0], frg[1])
		for n in range(len(pos_pos)):

			ref = ref_sub.split('@')[0] # extract subject seq name
			p_ref = ref_sub.split('@')[1] # extract query seq name

			ref_flag = ref + '@' + str(pos_pos[n][0]) # make a unique flag for sub seq [genome_name@position]
			sub_flag = p_ref +"@"+ str(pos_pos[n][1]) # make a unique flag for query seq
			
			"""
			e.g. ref_flag: GCA_003466205.fna@183073; sub_flag: GCA_002834165.fna@77800
			site 183073 one genome GCA_003466205 is the homologous site of 77800 on GCA_002834165
			"""
			
			s_unmerged_pos_pos = unmerged(pos_pos)[0]
			q_unmerged_pos_pos = unmerged(pos_pos)[1]

			# print(s_unmerged_pos_pos)  
			sub_ref_col = col_finder_2(aln_dict, ref, pos_pos[n][0], s_unmerged_pos_pos)
			q_ref_col=col_finder_2(aln_dict, p_ref, pos_pos[n][1], q_unmerged_pos_pos)
			# print(sub_ref_col, q_ref_col)
			test[ref_flag].append((ref_flag, sub_ref_col))
			test[ref_flag].append((sub_flag, q_ref_col))

	return test


def nproc_sampler(processor, ref_sub_lst):


	return multi_map(processor, ColSampler, ref_sub_lst)

# def site_checker(site_lst):
# 	"""
# 	It takes a list of homo sites. If all sites share same nucleotide
# 	if return the nucleotide otherwise it fills in gap
# 	"""
# 	if all(site == site_lst[0] for site in site_lst) == True:
# 		return site_lst[0]
# 	else:
# 		return '-'

def check_max_uniq(site_lst):

	votes = {'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0}
	for s in site_lst:
		votes[s] += 1

	return max(votes, key=votes.get)

def consensus_col_builder(lst_homo_cols): # ? Check here, might be some problems.
	consensus_col = []
	for taxa in range(len(lst_homo_cols[0])):
		site_lst = []
		for col in range(len(lst_homo_cols)):
			site_lst.append(lst_homo_cols[col][taxa])
		# consensus_col.append(site_checker(site_lst))
		consensus_col.append(check_max_uniq(site_lst))
	
	return ''.join(consensus_col)


if __name__== '__main__':

	pars = read_args(sys.argv)
	refs_name_headers = refs_collector(pars['alns_folder']) # refenrece genome names in a list

	make_blast_db('Refs.fna') # make blast database
	blast_('Refs.fna', 'Refs.fna', 4) # call blast
	QC_on_blastn('INTER_blast_opt_tmp.tab', pars['length'], pid = pars['identity']) # QC on blast result

	aln_paths = [a for a in subprocess.getoutput('ls {}/*'.format(pars['alns_folder'])).split('\n')]

	RefSeqs_SubWise =SubWise_Refs(partition_genomes('INTER_blastn_opt_tmp_cleaned.tab'))	
	"""
	|
	|
	V
	RefSeqs_SubWise is dictory which builds a link between subject sequence and query sequences.
	Key is the subject sequence and corresonding query sequences are in the list as value.
	{'GCA_003466205.fna': ['GCA_002834165.fna', 'GCA_002834225.fna', 'GCA_002834235.fna', 'GCA_003466165.fna'],
	 'GCA_002834235.fna': ['GCA_002834165.fna', 'GCA_002834225.fna', 'GCA_003466165.fna', 'GCA_003466205.fna'],
	 'GCA_003466165.fna': ['GCA_002834165.fna', 'GCA_002834225.fna', 'GCA_002834235.fna', 'GCA_003466205.fna'],
	 'GCA_002834225.fna': ['GCA_002834165.fna', 'GCA_002834235.fna', 'GCA_003466165.fna', 'GCA_003466205.fna'],
	 'GCA_002834165.fna': ['GCA_002834225.fna', 'GCA_002834235.fna', 'GCA_003466165.fna', 'GCA_003466205.fna']}
	"""

	ref_sub_groups = [] # The list which holds unique labels. Thus, each sub seq and que seq is paired. 

	for ref in RefSeqs_SubWise:
		for sub in RefSeqs_SubWise[ref]:
			uniq_label = ref+'@'+sub # uniq_label is composed of ref seq name and sub seq name, delimited by '@'.
			ref_sub_groups.append(uniq_label)
		
	A = nproc_sampler(pars['nproc'], ref_sub_groups)

	a = defaultdict(list) # It merges a list of dictionaries into one dictionary
 
	for to_merge in A:
		for single_work in to_merge:
			a[single_work].extend(to_merge[single_work])

	all_cols_for_merging = [] # supercore

	for i in a: # i=genome@position; e.g. GCA_002834235.fna@306351
		x = de_duplicator(a[i])
		if len(x)/len(aln_paths) >= pars['coreness']:
			all_cols_for_merging.append(x)
		
	

	consensus_col_lst = []

	for i in all_cols_for_merging:
			consensus_col = consensus_col_builder(i)
			consensus_col_lst.append(consensus_col)
			

	def aln_builder(aln, all_cols_list): # aln is the template alignment file
		
		rec_lst = []
		reordered = sorted(record.id for record in aln)

		rotated = list(zip(*reversed(all_cols_list)))
		for taxa_idx in range(len(reordered)):
			seq = ''.join(rotated[taxa_idx])
			rec_lst.append(SeqRecord(Seq(seq, generic_dna), id = reordered[taxa_idx], description = ''))
		
		return rec_lst
#------------------------------------------------------------------------------
	#Distribute all columns into different files

	# for i in range(len(aln_paths)):
	# 	lst = [] # list of columns at one position  	
	# 	for j in all_cols:
	# 		lst.append(list(j[i]))

	# 	final = aln_builder(AlignIO.read('{}/{}'.format(pars['alns_folder'], refs_name_headers[0]), 'fasta'), lst)
	# 	SeqIO.write(final,'HomoAln_{}.fna'.format(str(i)), 'fasta')
#-------------------------------------------------------------------------------

	merged = aln_builder(AlignIO.read('{}/{}'.format(pars['alns_folder'], refs_name_headers[0]), 'fasta'), consensus_col_lst)
	SeqIO.write(merged, 'columned_merged_aln.fna', 'fasta')
# ------------------------------------------------------------------------------	

# Notes: 
#1. sampled coulmn redundancy can be solved by applying harder quality control on selection blastn results
#2. be careful, reference alignment file name should be consistent with reference genome name inside fasta
