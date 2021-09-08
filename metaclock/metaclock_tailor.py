#!/usr/bin/env python

"""
NAME: Python program metaclock_landscape.py - mutational landscape
Main features:
             1. Partition whole genome alignment into genetic blocks (i.e. coding and noncoding child alignments).
             2. Assess each child alignment and generate child alignment table characterized by various indices.
"""

import argparse
import sys, random
import subprocess
import os
import logging
import shutil
import collections
from Bio import AlignIO
from .utils import gff3_parser
import multiprocessing as mp
from collections import defaultdict
import itertools
import statistics
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import _Matrix
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
from .utils import utils
from .utils import PreClock_LRM
import numpy as np
import matplotlib.pyplot as plt
from ete3 import Tree
import matplotlib.gridspec as gridspec
from operator import itemgetter
from logging.config import fileConfig

log_config = "/".join(os.path.abspath(__file__).split('/')[:-1]) + '/metaclock_configs/logging_config.ini'
fileConfig(log_config)
logger = logging.getLogger()

def add_landscape_options(subparsers):

    landscape_parser = subparsers.add_parser('landscape',
                                             help = 'Tailor the whole genome alignment based on the genome landscape.',
                                             description = 'Assess mutational landscape of whole genome alignment and select interesting regions.')

    landscape_parser.add_argument('wga',
                        nargs = '?',
                        metavar = 'whole_genome_alignment',
                        help = 'Input the whole genome alignment generated by metaclock_mac.',
                        type = str)

    landscape_parser.add_argument('-o',
                        '--opt_dir',
                        help = 'Sepecify a directory name for outputs. default: "opt_dir" in current working directory.',
                        type = str,
                        default = 'opt_dir')

    landscape_parser.add_argument('-gff3',
                        '--gff3_annotation',
                        help = 'Input the gff3 annotation file. If no gff3 file was provided automatic \
                        annotation would be performed. (Please use gff file from the same reference \
                        sequence in resconstructing whole genome alignment in order to maintain the consistency.)',
                        type = str,
                        default = None)

    landscape_parser.add_argument('-s',
                        '--select',
                        help = 'Flag this option to select your interesting regions. Note: other options need to be specificed.',
                        action = 'store_true')
 
    landscape_parser.add_argument('-nproc',
                        '--number_of_processors',
                        help = 'Specify the number of processor to be used. default: [1]',
                        type = int,
                        default = 1)

    landscape_parser.add_argument('-ns',
                        '--number_of_samples',
                        help = 'Specify the minimum number of samples with <10 percent missing information to select the child alignment. Default: [1]',
                        type = float,
                        default = 1)

    landscape_parser.add_argument('-l',
                        '--minimum_length',
                        help = 'Specify the minimum length of a child alignment to be selected. default: [100]',
                        type = int,
                        default = 100)

    landscape_parser.add_argument('-mi',
                        '--maximum_missing_information',
                        help = 'Specify the minimum missing information of a child alignment to be selected. default: [1.0]',
                        type = float,
                        default = 1.0)


    landscape_parser.add_argument('-snv_density_max',
                        '--maximum_snv_density',
                        help = 'Specify the maximum SNV density for selecting a child alignment in order to exclude hypervariable regions. default: [1.0]',
                        type = float,
                        default = 1)

    landscape_parser.add_argument('-snv_density_min',
                        '--minimum_snv_density',
                        help = 'Specify the minimum SNV density for selecting a child alignment in order to exclude hypervariable regions. default: [0.0]',
                        type = float,
                        default = 0.0)

    landscape_parser.add_argument('-Pdist_mean',
                        '--average_pairwise_distances',
                        help = 'Specify the maximum average pairwise distances for selecting a child alignment. default: [1.0]',
                        type = float,
                        default = 1.0)

    landscape_parser.add_argument('-Pdist_stdv',
                        '--pairwise_distances_stdv',
                        help = 'Specify the maximum standard deviation of pairwise distances. default: [1.0]',
                        type = float,
                        default = 1.0)
    
    landscape_parser.add_argument('-f',
                        '--feature',
                        help = 'Specify the types [CDS, tRNA, nCDS] to select. If more than one type comma should be used to delimit.\
                        e.g. CDS,tRNA,nCDS: Output all kinds of alignment.',
                        type = str,
                        default = 'CDS,tRNA,nCDS')
    
    landscape_parser.add_argument('-fast_TempEst',
                        '--fast_temporal_signal_test',
                        help = 'Inputting a mapping file containing sample names and sample dates will allow a fast estimation for temporal signal.\
                        Note: this option has to be used together with -raxml.',
                        type = str,
                        default = None)

    landscape_parser.add_argument('-raxml',
                        '--raxml_tree',
                        help = 'Reconstruct a maximum likelihood tree using raxmlHPC.',
                        action = 'store_true')

def add_basic_options(subparsers):

    basic_parser = subparsers.add_parser('basic',
                                             help = 'Tailor the whole genome alignment using the basic mode - missing information trimming.',
                                             description = 'Assess mutational landscape of whole genome alignment and select interested regions.')

    basic_parser.add_argument('wga',
                        nargs = '?',
                        metavar = 'whole_genome_alignment',
                        help = 'Input the whole genome alignment generated by metaclock_mac.',
                        type = str)

    basic_parser.add_argument('opt',
                    nargs = '?',
                    metavar = 'output_file',
                    help = 'Specify output name for the tailored genome alignment in fasta format.',
                    type = str)

    basic_parser.add_argument('-a',
                    "--automated_tailoring",
                    help = "This option allows automated tailoring using trimAl for minimizing missing information in the alignment",
                    action = 'store_true')

    basic_parser.add_argument('-s',
                    "--stats",
                    help = "Specify a file name for statistic report of tailored alignment.",
                    type = str,
                    default = 'stats.tsv')

    basic_parser.add_argument('-tt',
                    "--target_tailoring",
                    help = "Input a list of samples to keep for alignment tailoring. The limit for missing info in each column can\
                     be given and it is delimited by comma (e.g. shortlist.txt,0.1), or column tailoring will be performed automatically using trimAl (e.g. list.txt).",
                    type = str)
    basic_parser.add_argument('-raxml',
                    '--raxml_output',
                    help = 'Specify a folder name for storing ML phylogeny results. Default: None',
                    type = str,
                    default = None)
    basic_parser.add_argument('-nproc',
                    '--number_of_processors',
                    help = 'Specify the number of processor to be used. default: [1]. Note: This is for running ML phylogeny.',
                    type = int,
                    default = 1)


    

def tailoring_cols(Seq_dict, c):
    
    core_col = []
    Seq_lst = [key for key in Seq_dict]
    Seq_lst_lst = [list(Seq_dict[g].seq) for g in Seq_lst]

    for C in range(len(list(Seq_dict[Seq_lst[0]].seq))):
        col_sites = [i[C] for i in Seq_lst_lst]
        tot_sites_num = len(col_sites)
        if float(col_sites.count('-')/tot_sites_num) > c+0.00001:
            continue
        else:
            core_col.append(C)
    
    rec_list=[]
    for g in Seq_dict:
        seq_lst = list(Seq_dict[g].seq)
        core_seq = itemgetter(*core_col)(seq_lst)
        seq = "".join(core_seq)
        _id = g
        rec_list.append(SeqRecord(Seq(seq), id = _id, description = ''))

    return rec_list


class tailor(object):

    def __init__(self, aln):
        self.aln = aln

    def short_list_tailor(self, short_list, c):
        # c is the cutoff for column
        # If c is not specified Trimal automated trimming will be activated
        Seq_dict = SeqIO.to_dict(self.aln)
        selected_aln = []
        for i in Seq_dict:
            if i in short_list:
                selected_aln.append(SeqRecord(Seq(str(Seq_dict[i].seq)), i, description = ''))
            else:
                continue
        if c:
            selected_aln_dict = SeqIO.to_dict(selected_aln)
            return tailoring_cols(selected_aln_dict, c)
        else:
            # Here trimal gappyout is implemented. 
            # Make sure trimal is set in your system
            SeqIO.write(selected_aln, 'tmp_feed_to_trimal.fna', 'fasta')
            cmd = 'trimal -gappyout -in tmp_feed_to_trimal.fna -out TrimalGappyout_opt.fna'
            subprocess.call(cmd, shell = True)


        

def short_list_tailaring_pars_check(par):
    if ',' in par:
        par_split = [p for p in par.split(',') if len(p) > 0]
        if len(par_split) == 2:
            return (par_split[0], float(par_split[1])+0.0000000001)
        else:
            sys.exit('Please specify short list file and gap cutoff for each column. e.g. text.txt,0.1')
    else:
        return (par, None)

def parse_gff(wga_aln, gff_file):

    """
    This function pack all genetic content (nCDS and CDS) into a dictionary.
    key: genetic content label - [type]_[start]_[end] * zero-based and python slicing style
    value: a named tuple (type, function, child_aln)
    """
    # Uniform gff into zero_based and adjust gene coordinates to python style in slicing.
    parsed_gff3 = gff3_parser.Parser(gff_file).element()
    whole_genome_coordinates = list(range(len(wga_aln[1,:])))
    gc_base = {} # gc_base is a dictionary: genetic content label is the key and the namedtiple is the value
    genetic_contents = collections.namedtuple('genetic_contents', 'Type Function Child_aln')
    obj_list = list(parsed_gff3)
    first_nCDS = genetic_contents(Type = 'nCDS',\
     Function = 'NA', Child_aln = wga_aln[:, 0: int(obj_list[0].start) -1])
    first_nCDS_label = 'nCDS' + '_0_' + str(int(obj_list[0].start) -1) 
    gc_base[first_nCDS_label] = first_nCDS 

    for i in range(len(obj_list)-1):
        I = obj_list[i]
        II = obj_list[i+1]
        if I.type != 'repeat_region':

            cds_label = I.type + '_' + str(int(I.start)-1) + '_' + str(int(I.end))
            cds_slicing_idx = (int(I.start)-1, int(I.end))
            cds_type = I.type
            cds_product = I.attributes['product']
            ncds_label = 'nCDS_' + str(int(I.end)) + '_' + str(int(II.start) - 1)
            ncds_slicing_idx = (int(I.end), int(II.start) - 1)
            ncds_type = 'nCDS'
            ncds_product = 'NA'
            gc_base[cds_label]=genetic_contents(Type = cds_type,\
             Function = cds_product, Child_aln = wga_aln[:, cds_slicing_idx[0]: cds_slicing_idx[1]])
            gc_base[ncds_label]=genetic_contents(Type = ncds_type,\
             Function = ncds_product, Child_aln = wga_aln[:, ncds_slicing_idx[0]: ncds_slicing_idx[1]])
            
    last_cds_label = obj_list[-1].type + '_' + str(int(obj_list[-1].start)-1) + '_' + obj_list[-1].end
    last_cds_slicing_idx = (int(obj_list[-1].start)-1, int(obj_list[-1].end))
    last_cds_type = obj_list[-1].type
    last_cds_product = obj_list[-1].attributes['product']
    last_ncds_label='nCDS_'+ obj_list[-1].end + '_' + str(len(wga_aln[1,:]))
    last_ncds_slicing_idx = (int(obj_list[-1].end), len(wga_aln[1,:]))
    last_ncds_type = 'nCDS'
    last_ncds_product = 'NA'
    gc_base[last_cds_label] = genetic_contents(Type = last_cds_type,\
     Function = last_cds_product, Child_aln = wga_aln[:, last_cds_slicing_idx[0]: last_cds_slicing_idx[1]])
    gc_base[last_ncds_label] = genetic_contents(Type = last_ncds_type,\
     Function = last_ncds_product, Child_aln = wga_aln[:, last_ncds_slicing_idx[0]: last_ncds_slicing_idx[1]])

    return gc_base

        
class ChildAlnStats(object):

    def __init__(self, child_aln_obj):
        self.child_aln_obj = child_aln_obj # It handles AlignIO object
    def aln_len(self):
        child_aln_len = len(self.child_aln_obj[1,:])
        return child_aln_len
    def variable_sites(self):
        # It calculates alignment longer than 100bp
        multi_count = 0
        bi_count = 0
        snv = collections.namedtuple('snv', 'SNV_no SNV_density')
        if self.aln_len() >= 10:
            for j in range(len(self.child_aln_obj[1])):
                col = set([i for i in list(self.child_aln_obj[:,j]) if i != '-'])
                if len(col) == 2:
                    bi_count += 1
                elif len(col) > 2:
                    multi_count += 1
            all_count = bi_count + multi_count
            if all_count > 0:
                snv_density = all_count/self.aln_len()
                # bi_ratio = bi_count/all_count
            else:
                # bi_ratio = 0.0
                snv_density = 0.0
            return snv(SNV_no = all_count, SNV_density = snv_density)
        else:
            return snv(SNV_no = 'NA', SNV_density = 'NA')

    def avg_pw_dist(self):
        dist_value = collections.namedtuple('dist_value', 'Avg Stdv')
        if self.variable_sites().SNV_no != 'NA':
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(self.child_aln_obj)
            pw_dists = [y for x in dm.matrix for y in x[:-1]] # flatten a nested list into one list
            # and remove self comparison
            avg = statistics.mean(pw_dists)
            stdv = statistics.stdev(pw_dists)
            return dist_value(Avg = avg, Stdv = stdv)
        else:
            return dist_value(Avg = 'NA', Stdv = 'NA')

    def tm_pw_dist(self, time_differences):
        dist_value = collections.namedtuple('dist_value', 'coefficint b_mu_rate')
        if self.variable_sites().SNV_no != 'NA':
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(self.child_aln_obj)

            genetic_diff = []
            time_diff = []

            for t in time_differences:
                pair_1, pair_2 = t[0][0], t[0][1]
                t_difference = t[1]
                g_difference = dm[pair_1, pair_2]

                genetic_diff.append(g_difference)
                time_diff.append(t_difference)
            df_ = pd.DataFrame.from_dict({'g': genetic_diff, 't': time_diff})
            col_g = df_['g']
            col_t = df_['t']
            col_eff = col_t.corr(col_g)
            
            df_['mu_rate'] = df_['g']/df_['t']

            # tm_dist_avg = df_['mu_rate'].mean()
            tm_dist_avg = df_['g'].sum()/df_['t'].sum()

            return dist_value(coefficint = col_eff, b_mu_rate = tm_dist_avg)
        else:

            return dist_value(coefficint = 'NA', b_mu_rate = 'NA')



    def missing_value(self):
        if self.aln_len() >= 20:
            init = ''
            for i in self.child_aln_obj:
                init += i.seq
            gaps = init.count('-')
            tot_len = len(init)
            return gaps/tot_len
        else:
            return 'NA'
    def no_sample_minimum_samples(self, mv = 0.1):
        # Count the number of samples with minimum missing information
        # in each child alignment. Default: [0.05]
        if self.aln_len() >= 20:
            counter = 0
            for i in self.child_aln_obj:
                s_mv = i.seq.count('-')/len(i.seq)
                if s_mv <= mv:
                    counter += 1
            return str(counter)
        else:
            return 'NA'



def multi_dist(packed_args):
    # Namedtuple causes problem in multiprocessing
    # here use positional index

    label, aln = packed_args
    cAln_obj = ChildAlnStats(aln)   
    apw_dist = cAln_obj.avg_pw_dist()
    return (label, (apw_dist.Avg, apw_dist.Stdv))

def multi_tm_dist(packed_args):
    label, aln, td_lst = packed_args
    cAln_obj = ChildAlnStats(aln)
    tm_dist = cAln_obj.tm_pw_dist(td_lst)

    return (label, (tm_dist.coefficint, tm_dist.b_mu_rate))

def multi_variable_sites(packed_args):
    label, aln = packed_args
    cAln_obj = ChildAlnStats(aln)
    vs = cAln_obj.variable_sites()

    return (label, (vs.SNV_no, vs.SNV_density))

def dict_judge(label, dict_):
    value = list(dict_.values())[0]
    if label in dict_:
        return dict_[label]
    else:
        if isinstance(value, tuple):
            return ['NA']*len(value)
        else:
            return 'NA'

def conv_2_str(lst):

    converted_lst = [str(i) for i in lst]
    return converted_lst

def type_check(type_):
  base_set = set(['tRNA', 'CDS', 'nCDS'])

  if len(type_) == 3 and set(type_).issubset(base_set):
    feature = type_
  elif len(type_) == 2 and set(type_).issubset(base_set):
    type_.append('Holder')
    feature = type_
  elif len(type_) == 1 and set(type_).issubset(base_set):
    type_.append('Holder')
    type_.append('Holder')
    feature = type_
  else:
    sys.exit('Please specify the features [tRNA,CDS,nCDS] you want to select!')

  return feature



def aln_filter(a_tab, par_tuple):
    # Here also set min/max bound to make sure function works !
    # So here adds a parameter control to make sure inputs are within bounds


  type_, Length, SNV_density_max, SNV_density_min, missing_value, mini_sample, avg_dist, stdv_dist = par_tuple
  
  type_ = type_check(type_)
  
  len_range = (min(list(a_tab['Length'])), max(list(a_tab['Length'])))
  SNV_density_range = (min(list(a_tab['SNV_density'])), max(list(a_tab['SNV_density'])))
  missing_value_range = (min(list(a_tab['Missing_value'])), max(list(a_tab['Missing_value'])))
  mini_sample_range = (min(list(a_tab['#Samples(missing information < 10%)'])),\
   max(list(a_tab['#Samples(missing information < 10%)'])))
  avg_dict_range = (min(list(a_tab['Avg_genetic_distances'])), max(list(a_tab['Avg_genetic_distances'])))
  stdv_dict_range = (min(list(a_tab['Stdv_genetic_distances'])), max(list(a_tab['Stdv_genetic_distances'])))
  if any( [Length > len_range[1],\
           SNV_density_max < SNV_density_range[0],\
           missing_value < missing_value_range[0],\
           mini_sample > mini_sample_range[1],\
           avg_dist < avg_dict_range[0],\
           stdv_dist < stdv_dict_range[0]] ):
      sys.exit("Please be aware the range:\n\
             Length: {} - {}\n\
             SNV density: {} - {}\n\
             Missing value: {} - {}\n\
             #Sample(missing value < 5 percent): {} - {}\n\
             Average distance: {} - {}\n\
             Stdv of average distance: {} - {}".format(len_range[0], len_range[1],\
                SNV_density_range[0], SNV_density_range[1],\
                missing_value_range[0], missing_value_range[1], mini_sample_range[0],\
                mini_sample_range[1], avg_dict_range[0], avg_dict_range[1],\
                stdv_dict_range[0], stdv_dict_range[1]))
  
  else:
      a_tab = a_tab.loc[a_tab['Length'] >= Length]
      a_tab = a_tab.loc[(a_tab['Type'] == type_[0]) | (a_tab['Type'] == type_[1]) | (a_tab['Type'] == type_[2])]
      a_tab = a_tab.loc[(a_tab['SNV_density'] <= SNV_density_max) & (a_tab['SNV_density'] >= SNV_density_min)]
      a_tab = a_tab.loc[a_tab['Missing_value'] <= missing_value]
      a_tab = a_tab.loc[a_tab['#Samples(missing information < 10%)'] >= mini_sample]
      a_tab = a_tab.loc[a_tab['Avg_genetic_distances'] <= avg_dist]
      a_tab = a_tab.loc[a_tab['Stdv_genetic_distances'] <= stdv_dist]

      return a_tab  

def Assessment_tab(wga_aln, gff3, nproc):
    logger.info("Start parsing gff3 file.....")
    gContent_db = parse_gff(wga_aln, gff3)
    logger.info("Parsing is completed!")

    nCDS_length = [len(gContent_db[i].Child_aln[1,:]) for i in gContent_db if i.startswith('nCDS')]
    coding_percent = str(100*(1 - sum(nCDS_length)/len(wga_aln[1,:])))+'%'
    label_aln_pairs = [(i, gContent_db[i].Child_aln) for i in gContent_db]
    
    logger.info("Assess lengths of partitioned alignments....")
    len_dict = {a: ChildAlnStats(gContent_db[a].Child_aln).aln_len() for a in gContent_db}
    logger.info("Length assessment is completed!")
    
    logger.info("Assess snv density for each partitioned alignment....")
    vs_dict = utils.multi_proc_dict(nproc, multi_variable_sites, label_aln_pairs)
    logger.info("snv density assessment is completed!")
    
    logger.info("Assess missing information in each partitioned alignment....")
    mv_dict = {a: ChildAlnStats(gContent_db[a].Child_aln).missing_value() for a in gContent_db}
    logger.info("Missing information assessment is completed!")
    
    logger.info("Assess the number of samples in the partitioned alignment having missing information < 10%")
    no_samples_min_mv_dict = {a: ChildAlnStats(gContent_db[a].Child_aln).no_sample_minimum_samples(0.1)\
     for a in gContent_db}
    logger.info("#Samples assessment is completed!")

    logger.info("Assess average pairwise genetic distances for each partitioned alignment....")
    avg_pwd_dict = utils.multi_proc_dict(nproc, multi_dist, label_aln_pairs)
    logger.info("Genetic assessment is completed!")
        # {'label': (mean, stdv)}                

    dict_lst_BigMatrix = defaultdict(list)            

    for label in gContent_db:

        if len_dict[label] >= 100: 
            dict_lst_BigMatrix[label].append([label, gContent_db[label].Type, gContent_db[label].Function, len_dict[label],\
             dict_judge(label, vs_dict)[1], dict_judge(label, mv_dict), dict_judge(label, no_samples_min_mv_dict),\
                dict_judge(label, avg_pwd_dict)[0], dict_judge(label, avg_pwd_dict)[1]])            
    
    return dict_lst_BigMatrix

def create_folder(name):
    folder_name = os.path.abspath(name) # determine the abs path
    if os.path.exists(folder_name):
        logger.info('{} exists: pass!'.format(folder_name))
        pass
    else:
        os.makedirs(folder_name)
        logger.info('Creating folder {}'.format(folder_name))
    return folder_name

def run_raxml(ipt_aln, nproc, opt_dir):

    opt_name = ipt_aln.split('/')[-1].replace('.fna', '')
    cmd = 'raxmlHPC-PTHREADS-SSE3 -T {} -f a -# 100 -p 12345 -x 12345 -s {} -m GTRGAMMA -n {} -w {}'.format(nproc, ipt_aln, opt_name, opt_dir)
    logger.info('{}'.format(cmd))
    subprocess.call(cmd, shell = True)

    return opt_dir + '/RAxML_bipartitions.' + opt_name

def visual(tab, opt):
    df_ = pd.read_csv(tab, sep = '\t')
    v1 = df_['Length']
    v2 = df_['SNV_density']
    v3 = df_['Missing_value']
    v4 = df_['#Samples(missing information < 10%)'].dropna().astype(int)
    v5 = df_['Avg_genetic_distances']
    v6 = df_['Stdv_genetic_distances']




    v4_delimits = np.asarray(range(0, v4.max()))
    CDFs = np.asarray([utils.EvalCdf(v4, i) for i in v4_delimits])

        
    fig = plt.figure(figsize = (15, 15), constrained_layout=True)
    spec = gridspec.GridSpec(ncols=3, nrows=2, figure=fig)
    fig_ax1 = fig.add_subplot(spec[0, 0])
    fig_ax1.set_xlabel('Length (bp)', fontsize=16)
    fig_ax1.set_ylabel('Count', fontsize = 16)
    fig_ax1.tick_params(axis='both', which='major', labelsize=13)
    utils.hist_plot(fig_ax1, v1, {'bins': 200})
    
    fig_ax2 = fig.add_subplot(spec[0, 1])
    fig_ax2.set_xlabel('SNV density', fontsize=16)
    fig_ax2.set_ylabel('Count', fontsize = 16)
    fig_ax2.tick_params(axis='both', which='major', labelsize=13)
    utils.hist_plot(fig_ax2, v2, {'bins': 200})
    
    fig_ax3 = fig.add_subplot(spec[0, 2])
    fig_ax3.set_xlabel('Gap score', fontsize=16)
    fig_ax3.set_ylabel('Count', fontsize = 16)
    fig_ax3.tick_params(axis='both', which='major', labelsize=13)
    utils.hist_plot(fig_ax3, v3, {'bins': 200})
    
    fig_ax4 = fig.add_subplot(spec[1, 0])
    fig_ax4.set_xlabel('#Samples (Gap score < 10%)', fontsize = 16)
    fig_ax4.set_ylabel('CDF (Cumulative distribution frequency)', fontsize = 16)
    fig_ax4.tick_params(axis = 'both', which = 'major', labelsize = 13)
    utils.hist_plot(fig_ax4, v4, {'bins': v4.max(), 'histtype': 'step', 'cumulative': True})

    fig_ax5 = fig.add_subplot(spec[1, 1])
    fig_ax5.set_xlabel('Average pairwise genetic distances', fontsize = 16)
    fig_ax5.set_ylabel('Count', fontsize = 16)
    fig_ax5.tick_params(axis = 'both', which = 'major', labelsize = 13)
    utils.hist_plot(fig_ax5, v5, {'bins': 200})

    fig_ax6 = fig.add_subplot(spec[1, 2])
    fig_ax6.set_xlabel('Standard deviation of pairwise genetic distances', fontsize = 16)
    fig_ax6.set_ylabel('Count', fontsize = 16)
    fig_ax6.tick_params(axis = 'both', which = 'major', labelsize = 13)
    utils.hist_plot(fig_ax6, v6, {'bins': 200})




    fig.savefig(opt)


def main():

    parser = argparse.ArgumentParser('basic', 'landscape')
    subparsers = parser.add_subparsers(help = 'program mode', dest = 'mode')
    
    add_basic_options(subparsers)
    add_landscape_options(subparsers)
    

    args = parser.parse_args()

    if args.mode == 'basic':

        wga_aln = AlignIO.read(args.wga, 'fasta')
        aln_obj = tailor(wga_aln)

        if args.target_tailoring:
            logger.info("Tailoring is being peformed on target samples.")
            sl_args = short_list_tailaring_pars_check(args.target_tailoring)
            short_list = [i.rstrip() for i in open(sl_args[0]).readlines()]
            c = sl_args[1]

            if c:
                logger.info("The maximum gaps allowed in each column is: {}".format(str(c)))
                SeqIO.write(aln_obj.short_list_tailor(short_list, c), args.opt ,'fasta')
            else:
                logger.info("Column-wise gaps are being trimmed usng trimAl.")
                aln_obj.short_list_tailor(short_list, c)
                subprocess.call('mv TrimalGappyout_opt.fna {}'.format(args.opt), shell = True)
                subprocess.call('rm tmp_feed_to_trimal.fna', shell = True)
            if args.stats:
                utils.out_stats(args.opt, args.stats)
            
            if args.raxml_output:
                logger.info('Start inferring maximum likelihood tree....')
                logger.info('ML phylogeny is output in: {}'.format(args.raxml_output))
                opt_dir = create_folder(args.raxml_output)
                run_raxml(args.opt, args.number_of_processors, opt_dir)
                logger.info('RAxML tree reconstruction is completed and outputs are in: {}'.format(opt_dir))
            else:
                logger.info('ML phylogeny reconstruction was skipped!')

            logger.info("Thanks for using and welcome back!")

        elif args.automated_tailoring:
            logger.info("Automated tailoring using trimAl is being performed only on columns and no samples have been removed.")
            cmd = 'trimal -gappyout -in {} -out TrimalGappyout_opt.fna'.format(args.wga)
            subprocess.call(cmd, shell = True)
            subprocess.call('mv TrimalGappyout_opt.fna {}'.format(args.opt), shell = True)
            logger.info("Thanks for using and welcome back!")
            utils.out_stats(args.opt, opt_stats = args.stats)

            if args.raxml_output:
                logger.info('Start inferring maximum likelihood tree....')
                logger.info('ML phylogeny is output in: {}'.format(args.raxml_output))
                opt_dir = create_folder(args.raxml_output)
                run_raxml(args.opt, args.number_of_processors, opt_dir)
                logger.info('RAxML tree reconstruction is completed and outputs are in: {}'.format(opt_dir))
            else:
                logger.info('ML phylogeny was skipped!')

        else:
            sys.exit('Please choose one of these tailor strategies: automated tailoring or target tailoring.')

    elif args.mode == 'landscape':

        wga_aln = AlignIO.read(args.wga, 'fasta')
        opt_dir = create_folder(args.opt_dir)
        file_name = opt_dir + '/assessment.txt'
        opt_file_name = open(file_name, 'w')

        opt_file_name.write('Child_alignment\t' + 'Type\t'+ 'Function\t' + 'Length\t' + 'SNV_density\t'\
          + 'Missing_value\t' + '#Samples(missing information < 10%)\t'+'Avg_genetic_distances\t'\
          + 'Stdv_genetic_distances' + '\n')

        if args.gff3_annotation:
            logger.info("The given annotation profile: {}".format(args.gff3_annotation))

            dict_lst_BigMatrix = Assessment_tab(wga_aln, args.gff3_annotation, args.number_of_processors)
        
        else:
            logger.info("No gff3 file is given, so automated annotation starts.....")

            opt_RefSeq = random.choice([i for i in wga_aln if i.seq.count('-') == 0])
            logger.info('The most complete genome {} is chosen from the alignment for annotation.'.format(opt_RefSeq.id))
            opt_RefSeq_abs = os.path.abspath(args.opt_dir + '/'+ opt_RefSeq.id)
            SeqIO.write(opt_RefSeq, opt_RefSeq_abs, 'fasta')
            logger.info('The genome sequence for annotation is in: {}'.format(opt_RefSeq_abs))
            
            prokka_opt_dir = os.path.abspath(args.opt_dir + '/{}'.format(opt_RefSeq.id) + '.prokka')

            prokka_file_prefix = opt_RefSeq.id

            cmd = 'prokka {} --outdir {} --prefix {}'.format(opt_RefSeq_abs, prokka_opt_dir, prokka_file_prefix)
            logger.info('prokka annotation starts:')
            logger.info(cmd)
            subprocess.call(cmd, shell = True)
            gff3_annotation_file = os.path.abspath(prokka_opt_dir + '/' + prokka_file_prefix + '.gff')
            logger.info('The resulted annotation file: {}'.format(gff3_annotation_file))
            
            dict_lst_BigMatrix = Assessment_tab(wga_aln, gff3_annotation_file, args.number_of_processors)
            subprocess.call('rm {}'.format(opt_RefSeq_abs), shell = True)
            subprocess.call('rm -rf {}'.format(prokka_opt_dir), shell = True)
        
        logger.info('Start outputting assessed results....')
        for aln in dict_lst_BigMatrix:
            write_in_lines = "\t".join(conv_2_str(dict_lst_BigMatrix[aln][0]))
            opt_file_name.write(write_in_lines + '\n')
        logger.info('Assessment is completed and output is in: {}'.format(file_name))
        
        opt_file_name.close()

        opt_stats = opt_dir + '/assess_stats.png'

        logger.info('Start visualizing assessment....')
        visual(file_name, opt_stats)
        logger.info('Visualization is completed and output is in: {}'.format(opt_stats))

        if args.select:
            logger.info("Start selecting prefered regions....")
            tab_df = pd.read_csv(file_name, sep = '\t')
            feature_ = args.feature.split(',')
            pars = (feature_, args.minimum_length, args.maximum_snv_density, args.minimum_snv_density,\
                    args.maximum_missing_information, args.number_of_samples, args.average_pairwise_distances,\
                    args.pairwise_distances_stdv)


            selected_alns = aln_filter(tab_df, pars)
            logger.info('Selecting is completed.')

            concatenate_selected_alns = os.path.abspath(args.opt_dir + '/selected_region_concatenation.fna')
            
            logger.info('Concatenation starts....')
            init_aln = wga_aln[:, 0:1]
            for i in selected_alns.index:
                aln_header = selected_alns.loc[i,:]["Child_alignment"]
                slicing_index = (int(aln_header.split('_')[1]), int(aln_header.split('_')[2]))
                sub_aln = wga_aln[:,slicing_index[0]: slicing_index[1]]
                init_aln += sub_aln
            init_aln = init_aln[:, 1:]
            SeqIO.write(init_aln, concatenate_selected_alns, 'fasta')
            logger.info('Concatenation is finished and the output is in: {}'.format(concatenate_selected_alns))

            if args.raxml_tree:
                logger.info('Start inferring maximum likelihood tree....')
                opt_dir = create_folder(args.opt_dir + '/raxml')
                best_tree = run_raxml(concatenate_selected_alns, args.number_of_processors, opt_dir)
                logger.info('RAxML tree reconstruction is completed and outputs are in: {}'.format(opt_dir))

                if args.fast_temporal_signal_test:
                    logger.info('Start fast estimation of temporal signal.....')
                    logger.info('Use raxml tree from: {}'.format(best_tree))
                    opt_png = concatenate_selected_alns.replace('.fna', '') + '_TempEst.png'
                    mp_file = args.fast_temporal_signal_test
                    PreClock_LRM.temp_est(best_tree, mp_file, opt_png)



if __name__ == '__main__':

    main()
