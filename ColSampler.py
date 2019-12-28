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
		                help = 'minimum length of homologous sequence alignment in blast when\
		                 searching for coordinates between references.',
		                 type = str,
		                 default = '500')
	parser.add_argument('-i',
						'--identity',
						help = 'minium identity of homologous sequences when searching for\
						 coordinates between references.',
						 type = str,
						 default = '95.0')
	return vars(parser.parse_args())



def refs_collector(alns_folder):
	"""
	It takes a folder which contains single alignment and writes all reference
	genomes into one fasta file for blast procedure, meanwhile it returns a list
	of reference headers. 
	"""
	rec_lst = []
	ref_header_lst = []
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
	It partitions hits of genomes in pairwise style, ecluding self comparison. It 
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
	alignment start position, and return the postion of homologous
	sites on genome before alignment.

	$$$$ This part could be paralleled.
	"""
	s_pos_lst = []
	q_pos_lst = []
	for n in range(1,len(sseq)+1):
		idx = n - 1
		if list(sseq)[idx] != '-':
			s_pos_lst.append(s_idx)
			s_idx += 1
			if list(qseq)[idx] != '-':
				q_pos_lst.append(q_idx)
				q_idx += 1
			else:
				q_pos_lst.append('-')
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



def check_max_uniq(site_lst):
	"""
	This process could be paralleled 
	"""

	votes = {'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0}
	for s in site_lst:
		votes[s] += 1
	return max(votes, key=votes.get)


def consensus_col_builder(tuple_of_homo_cols):
	# It calculates consensus column using majority rule
	consensus_col = []
	for taxa in range(len(tuple_of_homo_cols[0])):
		site_lst = []
		for col in range(len(tuple_of_homo_cols)):
			site_lst.append(tuple_of_homo_cols[col][taxa])
		consensus_col.append(check_max_uniq(site_lst))

	return ''.join(consensus_col)		


def reorder_aln(aln):
	# it reorders headers in alignment object
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
	
	Ref_aln = Refs_aln_dict[ref_name]
	tax_number = len(Ref_aln[:, 2])
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

def ColSampler(single_aln_lst, RefSeqs_SubWise, matrix):
	# Sampler samples columns from single alignments iteratively using homologous site
	# positions in corresponding alignment and build final column-concatenated alignment.

	aln_dict = {}
	for aln in single_aln_lst:
		aln_dict[aln.split('/')[1]]=reorder_aln(AlignIO.read(aln, 'fasta'))
	
	ref_ori_dict = defaultdict(list)

	test = defaultdict(list)
	sample_sites_recorder = [] # It records all sampled sites to avoid redundancy
	for ref in RefSeqs_SubWise:
		for p_ref in RefSeqs_SubWise[ref]:
			key_idx = ref+'@'+p_ref
			sub_matrix=find_coordinates_between_fragements(matrix[key_idx])
			for frg in sub_matrix:
				pos_pos = merge(frg[0], frg[1])
				# print(ref, p_ref, pos_pos)
				for n in range(len(pos_pos)):
					ref_flag = ref+ '@' + str(pos_pos[n][0])
					sub_flag = p_ref +"@"+ str(pos_pos[n][1])
					if ref_flag not in sample_sites_recorder:
						sample_sites_recorder.append(ref_flag)
						s_unmerged_pos_pos = unmerged(pos_pos)[0]
						q_unmerged_pos_pos = unmerged(pos_pos)[1]  
						sub_ref_col = col_finder_2(aln_dict, ref, pos_pos[n][0], s_unmerged_pos_pos)
						q_ref_col=col_finder_2(aln_dict, p_ref, pos_pos[n][1], q_unmerged_pos_pos)
						test[ref_flag].append((ref_flag, sub_ref_col))
						test[ref_flag].append((sub_flag, q_ref_col))
						sample_sites_recorder.append(sub_flag)
						sample_sites_recorder.append(sub_flag)
					else:
						continue


	return test

if __name__== '__main__':

	pars = read_args(sys.argv)
	refs_name_headers = refs_collector(pars['alns_folder'])
	aln_paths = [a for a in subprocess.getoutput('ls {}/*'.format(pars['alns_folder'])).split('\n')]
	make_blast_db('Refs.fna') # make blast database
	blast_('Refs.fna', 'Refs.fna', 4) # call blast
	QC_on_blastn('INTER_blast_opt_tmp.tab', pars['length'], pid = pars['identity']) # QC on blast result
	RefSeqs_SubWise =SubWise_Refs(partition_genomes('INTER_blastn_opt_tmp_cleaned.tab'))
	matrix = partition_genomes('INTER_blastn_opt_tmp_cleaned.tab')

	all_cols = []
	a = ColSampler(aln_paths, RefSeqs_SubWise, matrix)
	
	for i in a:
		x = de_duplicator(a[i])
		all_cols.append(list(consensus_col_builder(x)))

	def aln_builder(aln, all_cols_list):
		rec_lst = []
		reordered = sorted([record.id for record in aln])

		rotated = list(zip(*reversed(all_cols_list)))
		for taxa_idx in range(len(reordered)):
			seq = ''.join(rotated[taxa_idx])
			rec_lst.append(SeqRecord(Seq(seq, generic_dna), id = reordered[taxa_idx], description = ''))
		
		return rec_lst

	aln_builder(AlignIO.read('{}/{}'.format(pars['alns_folder'], refs_name_headers[0]), 'fasta'), all_cols)	
	final = aln_builder(AlignIO.read('{}/{}'.format(pars['alns_folder'], refs_name_headers[0]), 'fasta'), all_cols)

	SeqIO.write(final,'ColumnConcatenated_aln.fna', 'fasta')

	

# Next work: put homologous columns in one list and build consensus columns and concatenate them
# together, be careful when you re-assign taxa names to the sequences. Don`t make mistake.

# Notes: 
#1. sampled coulmn redundancy can be solved by applying harder quality control on selection blastn results
#2. be careful, reference alignment file name should be consistent with reference genome name inside fasta





