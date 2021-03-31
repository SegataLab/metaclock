#!/usr/bin/env python
    
import subprocess
import argparse
import os
import sys
import shutil
import logging
from collections import defaultdict
from collections import Counter
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from itertools import groupby
from operator import itemgetter
from .metaclock_mac import build_db_file
from .metaclock_mac import blast_genomes
from .metaclock_mac import QC_on_blastn
from .metaclock_mac import blastn_sort
from .metaclock_mac import complement_seq
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align import MultipleSeqAlignment
from .utils import utils
from collections import ChainMap
import itertools
from .utils import gff3_parser
from operator import itemgetter
from datetime import datetime
from .utils import prokka_annotation
import multiprocessing as mp
from logging.config import fileConfig



"""
NAME: metaclock_combiner.py

DESCRIPTION: metaclock_combiner.py is a python program which samples columns in the alignment
             based on coordinates of homolougous sites of references, and outputs homologous
             columns into seperate files. Those sampled files are then merged into one alignment.

             New feature: it allows to choose columns from CDS and non-CDS regions.  
"""

__author__ = "Kun D. Huang"
__version__ = "0.1"
__date__ = "21.01.2021"


log_config = "/".join(os.path.abspath(__file__).split('/')[:-1]) + '/metaclock_configs/logging_config.ini'
fileConfig(log_config)
logger = logging.getLogger()

def read_args(args):

    parser = argparse.ArgumentParser()
    parser.add_argument('alns_folder',
                        nargs = "?",
                        metavar = "alignments_folder",
                        help = 'Specify the folder which contains alignments built using\
                         single reference with metaclock_mac. Note: alignment file name should be consistent\
                         with reference genome name inside fasta',
                        type = str)

    parser.add_argument('-l',
                        '--length',
                        help = 'The minimum length of homologous sequence alignment to keep as a hit. default: [500]',
                        type = str,
                        default = '500')

    parser.add_argument('-i',
                        '--identity',
                        help = 'The minium identity of homologous sequences to keep as a hit. default: [95.0]',
                         type = str,
                         default = '95.0')

    parser.add_argument('-p',
                        '--nproc',
                        help = 'Number of processors to use',
                        default = 1,
                        type = int)

    parser.add_argument('-hf',
                        '--homo_site_in_refs',
                        help = 'The number of references site must be in to be homologous. default: [The total number of references used]',
                        default = None,
                        type = int)

    parser.add_argument('-d',
                        '--output_folder',
                        help = 'Specify a name for output folder. default: [outputs] in current directory',
                        default = 'outputs',
                        type = str)

    parser.add_argument('-int',
                        '--intermediate_folder',
                        help = 'Specify a name for intermediate folder. default: [intermediates] in current directory',
                        default = 'intermediates',
                        type = str)

    parser.add_argument('-vt',
                        '--voting_threshold',
                        help = 'Specify a threshold for majority rule in merging columns. default: [0.5]',
                        default = 0.51,
                        type = float)
    parser.add_argument('-cds',
                        '--CDS_region',
                        help = 'Consider sites only within coding sequences. Note: make sure prokka is installed in the system path.',
                        action = 'store_true')

    return vars(parser.parse_args())


class RangeDict(dict):
    """
    It convert the exact range which should be key to the range of number as key
    """
    def __getitem__(self, item):
        if type(item) != range: # or xrange in python2
            for key in self:
                if item in key:
                    return self[key]
        else:
            return super().__getitem__(item)

def refs_collector(alns_folder, opt_dir):
    
    """
    It takes a folder which contains single alignments and writes all reference
    genomes into one fasta file for blast procedure, meanwhile it returns a list
    of reference headers. 
    """
    rec_lst = [] # a list of SeqRecords, each for one reference genome
    ref_header_lst = [] # a list of reference names
    all_alns = subprocess.getoutput('ls {}/*'.format(alns_folder)).split('\n')

    for aln in all_alns:
        RefSeq_name = aln.split('/')[-1]
        aln_dict = SeqIO.to_dict(SeqIO.parse(aln,'fasta'))
        rec_lst.append(SeqRecord(Seq(str(aln_dict[RefSeq_name].seq)), id = RefSeq_name, description = ''))
        ref_header_lst.append(RefSeq_name)

    opt_file = opt_dir + '/RefSeqs.fna'
    SeqIO.write(rec_lst, opt_file, 'fasta')

    return {'RefSeq_header_list': ref_header_lst, 'RefSeqs': opt_file}



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


def merge_cds(list_1, list_2, list_3, list_4, list_5, list_6, list_7, list_8):

    merged_list = tuple(zip(list_1, list_2, list_3, list_4, list_5, list_6, list_7, list_8))
    return merged_list

def merge(list_1, list_2, list_3, list_4, list_5, list_6):

    merged_list = tuple(zip(list_1, list_2, list_3, list_4, list_5, list_6))
    return merged_list

def multi_map(w,f,l):
    """
    multiprocessor is a function taking an integer(number of processor), a defined function,
    and a list of works as arguments. 
    """
    # This is to multiprocess command line on shell
    pool = mp.Pool(processes = w)
    return pool.map(f,l)

def site_coding_profile_document(gff_file):
    """
    This function takes gff_file as inpput
    and document coding and non-coding info for each site of RefSeq
    in the form of dict in which site is key and coding/non-coding (1/0) is value.
    Zero-based.
    """
    RefSeq_length = int(open(gff_file).readlines()[1].rstrip().split(' ')[-1])
    # Extract whole genome length from second line of gff file
    parsed_gff3 = gff3_parser.Parser(gff_file).element()
    dict_ = {}
    col_coding_profile = {}

    for i in parsed_gff3:
        if int(i.start) - int(i.end) < 0:
            dict_[range(int(i.start) - 1, int(i.end))] = i.type
        else:
            dict_[range(int(i.end) - 1, int(i.start))] = i.type
 
    c_rangedict = RangeDict(dict_)

    return c_rangedict

def sites_indexer(ori_list, col_coding_profile):

#   """
#   This function index coding site and filetering six lists using same index
#   """
    
    return ["C" if col_coding_profile[i] != None else 'N' for i in ori_list]

def find_coordinates_between_fragements(ref_sub, sub_matrix, gff_folder = None):

    """
    It finds original postions of sseq and qseq for each pair of reference genomes
    """
    if gff_folder:

        q_coding_profile = site_coding_profile_document(gff_folder + '/' + ref_sub.split('@')[0]+'.gff')
        s_coding_profile = site_coding_profile_document(gff_folder + '/' + ref_sub.split('@')[1]+'.gff')

        s_pos_lst_genome = []
        q_pos_lst_genome = []
        s_site_lst_blast = []
        q_site_lst_blast = []
        s_flag_lst_blast = []
        q_flag_lst_blast = []
        s_flag_lst_coding = []
        q_flag_lst_coding = []

        for hit in sub_matrix:
            hit_object =blastn_sort(hit)
            hit_senario = [hit_object.sseq, hit_object.sstart, hit_object.send, hit_object.qseq\
            , hit_object.qstart, hit_object.qend]

            q_pos = utils.homo_site_mapper(hit_senario)[0]
            s_pos = utils.homo_site_mapper(hit_senario)[1]
            q_site = utils.homo_site_mapper(hit_senario)[2]
            s_site = utils.homo_site_mapper(hit_senario)[3]
            q_flag = utils.homo_site_mapper(hit_senario)[4]
            s_flag = utils.homo_site_mapper(hit_senario)[5]
            q_cd_flag = utils.sites_indexer(q_pos, q_coding_profile)
            s_cd_flag = sites_indexer(s_pos, s_coding_profile)



            s_pos_lst_genome.append(s_pos)
            q_pos_lst_genome.append(q_pos)
            s_site_lst_blast.append(s_site)
            q_site_lst_blast.append(q_site)
            s_flag_lst_blast.append(s_flag)
            q_flag_lst_blast.append(q_flag)
            s_flag_lst_coding.append(s_cd_flag)
            q_flag_lst_coding.append(q_cd_flag)




        return merge_cds(q_pos_lst_genome, s_pos_lst_genome, q_site_lst_blast,\
         s_site_lst_blast, q_flag_lst_blast, s_flag_lst_blast,\
          q_flag_lst_coding, s_flag_lst_coding)

    else:

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
            q_pos = utils.homo_site_mapper(hit_senario)[0]
            s_pos = utils.homo_site_mapper(hit_senario)[1]
            q_site = utils.homo_site_mapper(hit_senario)[2]
            s_site = utils.homo_site_mapper(hit_senario)[3]
            q_flag = utils.homo_site_mapper(hit_senario)[4]
            s_flag = utils.homo_site_mapper(hit_senario)[5]
        
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



def ColSampler(ref_sub, matrix, gff_folder = None):
    
    # ref_sub, matrix = multi_arg for nproc later
    """
    Sampler samples columns from single alignments iteratively using homologous site
    positions in corresponding alignment and build final column-concatenated alignment.
    """
    # single_aln_lst = [a for a in subprocess.getoutput('ls {}/*'.format(pars['alns_folder'])).split('\n')]

    # aln_dict = {}   
    # for aln in single_aln_lst:
    #     aln_dict[aln.split('/')[1]] = reorder_aln(AlignIO.read(aln, 'fasta'))
    
    """
    Generate a dictionary, aln file name is the key and taxa-reordered alignment is the value.
    It is for sampling column from alignment file.
    """


    if gff_folder:
        sub_matrix=find_coordinates_between_fragements(ref_sub, matrix[ref_sub], gff_folder)
        
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
            qcd_scd =  tuple(zip(frg[6], frg[7]))
            for n in range(len(pos_pos)):
                if qcd_scd[n][0] == 'C' and qcd_scd[n][1] == 'C':
                    q_ref = ref_sub.split('@')[0] # extract query seq name
                    s_ref = ref_sub.split('@')[1] # extract subject seq name
                    que_flag = q_ref + '@' + str(pos_pos[n][0])+'@'+ site_site[n][0]+'@'+ flag_flag[n][0]# make a unique flag for sub seq [genome_name@position]
                    sub_flag = s_ref +"@"+ str(pos_pos[n][1])+'@'+ site_site[n][1] +'@'+ flag_flag[n][1]# make a unique flag for query seq
                    homo_site_list.append((que_flag, sub_flag))
            else:
                continue
        return homo_site_list
    else:
        sub_matrix=find_coordinates_between_fragements(ref_sub, matrix[ref_sub])
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


def sampler_cds(args):

    seq_pair, matrice, aln_dict, gff_folder = args
    seq_pair_1 = seq_pair.split('@')[0]
    seq_pair_2 = seq_pair.split('@')[1]
    logger.info('Sampling homologous sites between {} and {}....'.format(seq_pair_1, seq_pair_2))
    sites_containner_sub = []
    col_sampled = ColSampler(seq_pair, matrice, gff_folder)
    no_col_sampled = str(len(col_sampled))
    logger.info('We found {} homologous sites between {} and {}'.format(no_col_sampled, seq_pair_1, seq_pair_2))
    logger.info('Start locating site coordinates on the genome for {} and {}....'.format(seq_pair_1, seq_pair_2))
    
    for c in col_sampled:
        sites = col_finder(aln_dict, c[0], c[1])
        que_tag = c[0]
        que_col = c[0]+'$'+sites[0]
        sub_col = c[1]+'$'+sites[1]
        sites_containner_sub.append((que_tag, que_col))
        sites_containner_sub.append((que_tag, sub_col))

    return sites_containner_sub

def sampler(args):


    seq_pair, matrice, aln_dict = args
    seq_pair_1 = seq_pair.split('@')[0]
    seq_pair_2 = seq_pair.split('@')[1]

    logger.info('Sampling homologous sites between {} and {}....'.format(seq_pair_1, seq_pair_2))

    sites_containner_sub = []
    col_sampled = ColSampler(seq_pair, matrice)
    no_col_sampled = str(len(col_sampled))
    logger.info('We found {} homologous sites between {} and {}'.format(no_col_sampled, seq_pair_1, seq_pair_2))
    logger.info('Start locating site coordinates on the genome for {} and {}....'.format(seq_pair_1, seq_pair_2))
    for c in col_sampled:
        sites = col_finder(aln_dict, c[0], c[1])
        que_tag = c[0]
        que_col = c[0]+'$'+sites[0]
        sub_col = c[1]+'$'+sites[1]
        sites_containner_sub.append((que_tag, que_col))
        sites_containner_sub.append((que_tag, sub_col))
    return sites_containner_sub


def nproc_sampler(processor, ref_sub_lst, matrix, aln_dict, gff_folder = None):

    logger.info("The whole data has been partitioned and processed using {} processors.".format(processor))

    if gff_folder:
        matrix_lst = [matrix] * len(ref_sub_lst)
        aln_dict_lst = [aln_dict] * len(ref_sub_lst)
        gff_folder_lst = [gff_folder] * len(ref_sub_lst)
        
        return multi_map(processor, sampler_cds, zip(ref_sub_lst, matrix_lst, aln_dict_lst, gff_folder_lst))
    
    else:

        matrix_lst = [matrix] * len(ref_sub_lst)
        aln_dict_lst = [aln_dict] * len(ref_sub_lst)

        return multi_map(processor, sampler, zip(ref_sub_lst, matrix_lst, aln_dict_lst))




def check_max_uniq(site_lst, majority_threshold):

    votes = {'A': 0, 'T': 0, 'G': 0, 'C': 0, '-': 0}
    for s in site_lst:
        if s in votes:
            votes[s] += 1
        else:
            consensus_ = '-'

    consensus_ = max(votes, key=votes.get)
    dominance_allele_ratio = votes[consensus_]/sum(votes[i] for i in votes)

    if dominance_allele_ratio >= majority_threshold:

        consensus = consensus_
    
    else:
        consensus = '-' 

    return consensus

def consensus_col_builder(lst_homo_cols, majority_threshold, refs_num = 2):
    lst_homo_cols = [c.split("$")[-1] for c in lst_homo_cols]
    consensus_col = []
    if len(lst_homo_cols) >= refs_num: 
        for taxa in range(len(lst_homo_cols[0])):
            site_lst = []
            for col in range(len(lst_homo_cols)):
                site_lst.append(lst_homo_cols[col][taxa])
            # consensus_col.append(site_checker(site_lst))
            consensus_col.append(check_max_uniq(site_lst, majority_threshold))
    
        return ''.join(consensus_col)
    else:
        pass

def aln_builder(aln, all_cols_list): # aln is the template alignment file   
    rec_lst = []
    reordered = sorted(record.id for record in aln)
    rotated = list(zip(*reversed(all_cols_list)))
    for taxa_idx in range(len(reordered)):
        seq = ''.join(rotated[taxa_idx])
        rec_lst.append(SeqRecord(Seq(seq), id = reordered[taxa_idx], description = ''))        
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

def creat_uniqe_columns(super_lst):

    super_dict = defaultdict(set)
    for one_process in super_lst:
        for one_col in one_process:
            col_label = one_col[0]
            col_seq = one_col[1]
            super_dict[col_label].add(col_seq)
    return super_dict

class column_toolkit(object):

    """
    Object column_toolkit is to handle all homologous sites
    """
    def __init__(self, homo_sites):
        self.homo_sites = homo_sites
    def merge_col_voting(self, majority_threshold = 0.51, number_refs = 2):
        merged_col_lst = [consensus_col_builder(self.homo_sites[c], majority_threshold, number_refs) for c in self.homo_sites]
        res = list(filter(None, merged_col_lst))
        
        homo_sites = str(len(res))
        logger.info('{} sites are shared by {} genomes.'.format(homo_sites, number_refs))
        logger.info('Consensus site was decided by voting threshold {}'.format(majority_threshold))
        return res

def abspath_finder(file):
    # Convert wahtever path into abs path
    return os.path.abspath(file)

def create_folder(name):
    folder_name = abspath_finder(name) # determine the abs path
    if os.path.exists(folder_name):
        logger.info('{} exists: pass!'.format(folder_name))
        pass
    else:
        os.makedirs(folder_name)
        logger.info('Creating folder {}'.format(folder_name))
    return folder_name

def time_now():

    now = datetime.now()
    current_time = now.strftime("[%H:%M:%S]")

    return current_time


def blastn_nproc(nproc):

    if nproc >= 8:
        return 7
    else:
        return nproc


def main():
    
    pars = read_args(sys.argv)
    alns = abspath_finder(pars['alns_folder'])
    logger.info('Check input alignments: {}....'.format(alns))
    opt_dir = create_folder(pars['output_folder'])
    inter = create_folder(pars['intermediate_folder'])
    
    single_aln_lst = subprocess.getoutput('ls {}/*'.format(alns)).split('\n')
    logger.info('Combining multiple alignments starts....')
    for n in single_aln_lst:
        logger.info(n)

    aln_dict = {}
    for aln in single_aln_lst:
        aln_dict[aln.split('/')[-1]] = reorder_aln(AlignIO.read(aln, 'fasta'))

    refs_collected = refs_collector(alns, inter)
    refs_header_names = refs_collected['RefSeq_header_list']
    refs_sequences = refs_collected['RefSeqs']

    total_number_refs = len(refs_header_names)

    ref_sub_groups = [i[0]+'@'+i[1] for i in itertools.combinations(refs_header_names, 2)]


    logger.info('Start building database for blastn....')
    db_file = build_db_file('contigs', refs_sequences, inter, 'samples') # Here argument 'samples' is just a place holder.
    logger.info('Database building is completed and database files prefix is: {}....'.format(db_file))

    logger.info('Start searching for homologous regions among reference genomes....')
    blastn_threads = blastn_nproc(pars['nproc'])
    raw_blastn_opt = blast_genomes(inter, db_file, refs_sequences, blastn_threads, True)
    logger.info('Search is completed and the raw results are in: {}'.format(raw_blastn_opt))

    logger.info('Start cleaning raw results from homologs searching....')
    cleaned_blastn = QC_on_blastn(raw_blastn_opt, pars['length'], pars['identity'], inter)
    logger.info('Cleaning is completed and the results are in {}....'.format(cleaned_blastn))

    logger.info('Start partitioning data for speeding up computtaion using multiple processors....')

    partioned_matrice = {pair: partition_genomes(cleaned_blastn)[pair] for pair in ref_sub_groups}


    if pars['CDS_region']:
        logger.info("Combining multiple aligbments using only protein coding sites....")
        logger.info('This process would take a bit long....')
        gff_dir = prokka_annotation.auto_annotation(single_aln_lst, inter)

        all_sites_nproc = nproc_sampler(pars['nproc'], ref_sub_groups, partioned_matrice, aln_dict, gff_dir)
    else:
        logger.info("Combining multiple aligbments using all sites....")
        logger.info("This process would take a bit long....")
        all_sites_nproc = nproc_sampler(pars['nproc'], ref_sub_groups, partioned_matrice, aln_dict)
    
    logger.info('Start cleaning up homologous sites found....')
    sites_containner = creat_uniqe_columns(all_sites_nproc)
    logger.info('Sites cleaning is completed!')

    logger.info('Start building consensus sites')
    col_tool_obj = column_toolkit(sites_containner)
    logger.info('Start building consensus alignment....')
    if pars['homo_site_in_refs']:
        logger.info('The homologous site is defined by that a site must be in {} references'.format(pars['homo_site_in_refs']))
        cols_merged = col_tool_obj.merge_col_voting(pars['voting_threshold'], pars['homo_site_in_refs'])
    else:
        logger.info('The homologous site is defined by that a site must be in {} references'.format(total_number_refs))
        cols_merged = col_tool_obj.merge_col_voting(pars['voting_threshold'], total_number_refs)

    logger.info('Consensus alignment reconstruction is completed!')
    merged = aln_builder(AlignIO.read('{}/{}'.format(alns, refs_header_names[0]), 'fasta'), cols_merged)
    opt_file = opt_dir + '/multiple_alignments_combination.fna'
    SeqIO.write(merged, opt_file, 'fasta')
    logger.info('Please check your consensus alignment in: {}'.format(opt_file))


if __name__== '__main__':

    main()




# Notes: 
#1. sampled coulmn redundancy can be solved by applying harder quality control on selection blastn results
#2. be careful, reference alignment file name should be consistent with reference genome name inside fasta
#3. Add multiple processing which I removed last time.