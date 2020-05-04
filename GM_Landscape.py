#!/usr/bin/env python

"""
NAME: Python program GM_Landscape.py - Genomic Mutational Landscape
Main features:
             1. Partition whole genome alignment into genetic blocks (i.e. coding and noncoding child alignments).
             2. Assess each child alignment and generate child alignment table characterized by various indices.
             3. Calculate selective power for each child alignment (dN/dS, Pi, TajimaD).
             4. Select good quality child alignment and run raxml tree individually.
             5. Calculate tree distance and cluster similar trees together to amplify evolutionary signal on whole genome.  
"""

import argparse
import sys
import subprocess
import os
import shutil
import collections
from Bio import AlignIO
from gff3_parser import Parser
import multiprocessing as mp
from collections import defaultdict
import itertools
import statistics
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import IUPAC, Gapped
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import _Matrix
from Bio import SeqIO
import pandas as pd
from utils import multi_proc_dict
import numpy as np
import matplotlib.pyplot as plt
from visual_utils import hist_plot
from BuildGenomeAln import build_raxml
from ete3 import Tree



def add_assessing_cmd_options(subparsers):
    # Add command options for Assessment feature
    assess_parser = subparsers.add_parser('Assessment',
                                           help = 'Assessing all child alignments.', 
                                           description = 'This feature screens whole genome alignment\
                                           genetic-element by genetic-element. It generates a table\
                                            in which each row depicts one child alignment.')
    assess_parser.add_argument('wga',
                                nargs = '?',
                                metavar = 'whole_genome_alignment',
                                help = 'Input the whole genome alignment generated in preceding steps.',
                                type = str)
    assess_parser.add_argument('opt_dir',
                                nargs = '?',
                                metavar = 'output_directory',
                                help = 'Sepecify a directory name for storing outputs',
                                type = str)
    assess_parser.add_argument('-gff3',
                                '--gff3_annotation',
                                help = 'Input the gff3 annotation file. If no gff3 file was provided automatic \
                                annotation would be performed. (Please use the sequence built in whole genome alignment \
                                for annotation to maintain consistence)',
                                type = str,
                                default = None)
    assess_parser.add_argument('-ml',
                                '--minimum_length',
                                help = 'Specify the minimum length of child alignment to assess. default: [100]',
                                type = int,
                                default = 100)
    assess_parser.add_argument('-min_mv',
                               '--minimum_missing_information',
                               help = 'Specify the minimum missing information of each sample in each child alignment\
                                to count the number of samples. Default: [0.05]',
                               type = float,
                               default = 0.05)
    assess_parser.add_argument('-tf',
                            '--time_file',
                            help = 'Input the tip times file. 1st column is name and 2nd column is time, delimited by tab.',
                            type = str,
                            default = None)

    assess_parser.add_argument('-nproc',
                           '--number_of_processors',
                           help = 'Specify the number of processor to be used. default: [1]',
                           type = int,
                           default = 1)


def add_select_cmd_options(subparsers):
	# Add command options for selecting child alignment
    select_parser = subparsers.add_parser('Select',
	                                       help = 'Select child alignments.',
	                                       description = 'This feature generates and selects interesting child alignments\
	                                       based on various indices, such as length, missing value, dN/dS etc.')

    select_parser.add_argument('wga',
	                            nargs = '?',
	                            metavar = 'whole_genome_alignment',
	                            help = 'Input the whole genome alignment generated in preceding steps.',
                                type = str)
    select_parser.add_argument('a',
                               nargs = '?',
                               metavar = 'Assessment',
                               help = 'Input the assessment file generated in Assessment phase.',
                               type = str)
    select_parser.add_argument('-o',
                               '--output_file',
                               help = 'Specify the file name to store output table. Default: [standard output]',
                               type = str,
                               default = None)
    select_parser.add_argument('-l',
                               '--length',
                               help = 'Specify the minimum length of child alignment. Default: [20bp]',
                               type = int,
                               default = 20)

    select_parser.add_argument('-mv',
                               '--missing_value',
                               help = 'Specify the maximum missing value in the child alignment. Default: [1]',
                               type = float,
                               default = 1.0)
    select_parser.add_argument('-ms_cmv',
                               '--minimum_no_sample_controled_mv',
                               help = 'Specify the minimum number of samples [whose missing value < 5 percent]in the child alignment. Default: [0]',
                               type = float,
                               default = 1)
    select_parser.add_argument('-vs',
                               '--variable_sites',
                               help = 'Specify the minimum number of variable sites. default: [0]',
                               type = int,
                               default = 0)

    select_parser.add_argument('-vs_density',
                               '--variable_sites_density',
                               help = 'Specify the maximum variable sites density. default: [1]',
                               type = float,
                               default = 1)

    select_parser.add_argument('-a_Pdist',
                               '--average_pairwise_distance',
                               help = 'Specify the maximum average pairwise distance. default: [0]',
                               type = float,
                               default = 1.0)

    select_parser.add_argument('-a_Pdist_stdv',
                               '--average_pairwise_distance_stdv',
                               help = 'Specify the maximum Stdv of average pairwise distance. default: [0]',
                               type = float,
                               default = 1.0)
    select_parser.add_argument('-min_cor',
                               '--minimum_corrlation',
                               help = 'Specify the minimum correlation. default: [0]',
                               type = float,
                               default = 0.0)
    select_parser.add_argument('-tm_avg_dist',
                               '--time_measured_avg_distance',
                               help = 'Specify the max time-measured avg genetic distance. default: [1.0]',
                               type= float,
                               default= 1.0)
    select_parser.add_argument('-f',
                                 '--feature',
                                 help = 'Specify the types to output [CDS, rRNA, tRNA, non_CDS].\
                                  If more than one type comma should be used to delimit.\
                                 e.g. CDS,rRNA,tRNA,non_CDS: Output all kinds of alignment.',
                                 type = str,
                                 default = 'CDS,rRNA,tRNA,non_CDS')

    select_parser.add_argument('-ia',
                               '--individual_alignments',
                               help = 'If this option is chosen selected individual alignments will be output\
                               in  the folder called individual_alignments.',
                               action = 'store_true')

    select_parser.add_argument('-c',
                               '--concatenate',
                               help = 'Specify an output file name where to store concatenated selected\
                               individual alignments.',
                               type = str,
                               default = None)

    select_parser.add_argument('-raxml_t',
                               '--raxml_threads',
                               help = 'Specify the thread number to run raxml tree for concatenated alignment \
                               if ML-tree is expected.',
                               type = str,
                               default = 0)
    select_parser.add_argument('-fast_TempEst',
                               '--fast_temporal_signal_est',
                               help = 'Inputting a mapping file will allow a fast estimation for temporal signal.\
                               Note: this option has to be selected together with -raxml_t.',
                               type = str,
                               default = None)

def add_visual_cmd_options(subparsers):
	# Add command options for visualizing stats
    visual_parser = subparsers.add_parser('Visual',
	                                       help = 'Visualizing assessment',
	                                       description = 'This featue visualizes assessment.')
    visual_parser.add_argument('at',
	                            nargs = '?',
	                            metavar = 'assessment_table',
	                            help = 'Input the assessment table generated in preceding steps.',
                                type = str)
    visual_parser.add_argument('-a',
                               '--all',
                               help = 'Automatic visualization of all assessed aspects.',
                               action = 'store_true')



def parse_gff(wga_aln, gff_file):

    """
    This function pack all genetic content (nCDS and CDS) into a dictionary.
    key: genetic content label - [type]_[start]_[end] * zero-based and python slicing style
    value: a named tuple (type, function, child_aln)
    """
    # Uniform gff into zero_based and adjust gene coordinates to python style in slicing.
    parsed_gff3 = Parser(gff_file).element()
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
        if self.aln_len() >= 10:
            init = ''
            for i in self.child_aln_obj:
                init += i.seq
            gaps = init.count('-')
            tot_len = len(init)
            return gaps/tot_len
        else:
            return 'NA'
    def no_sample_minimum_samples(self, mv = 0.05):
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
def aln_filter(a_tab, par_tuple):
	# Here also set min/max bound to make sure function works !
	# So here adds a parameter control to make sure inputs are within bounds


	Length, SNV, SNV_density, missing_value, mini_sample, avg_dist, stdv_dist, cor, tm_dist = par_tuple

	len_range = (min(list(a_tab['Length'])), max(list(a_tab['Length'])))
	SNV_range = (min(list(a_tab['#SNV'])), max(list(a_tab['#SNV'])))
	SNV_density_range = (min(list(a_tab['SNV_density'])), max(list(a_tab['SNV_density'])))
	missing_value_range = (min(list(a_tab['Missing_value'])), max(list(a_tab['Missing_value'])))
	mini_sample_range = (min(list(a_tab['#Samples(missing information < 5%)'])),\
	 max(list(a_tab['#Samples(missing information < 5%)'])))
	avg_dict_range = (min(list(a_tab['Avg_genetic_distance'])), max(list(a_tab['Avg_genetic_distance'])))
	stdv_dict_range = (min(list(a_tab['Stdv_genetic_distance'])), max(list(a_tab['Stdv_genetic_distance'])))
	
	if any( [Length > len_range[1],\
	         SNV > SNV_range[1],\
	         SNV_density < SNV_density_range[0],\
	         missing_value < missing_value_range[0],\
	         mini_sample > mini_sample_range[1],\
	         avg_dist < avg_dict_range[0],\
	         stdv_dist < stdv_dict_range[0]] ):

	    sys.exit("Please be aware the range:\n\
	    	     Length: {} - {}\n\
	    	     #SNV: {} - {}\n\
	    	     SNV density: {} - {}\n\
	    	     Missing value: {} - {}\n\
	    	     #Sample(missing value < 5 percent): {} - {}\n\
	    	     Average distance: {} - {}\n\
	    	     Stdv of average distance: {} - {}".format(len_range[0], len_range[1],\
	    	     	SNV_range[0], SNV_range[1], SNV_density_range[0], SNV_density_range[1],\
	    	     	missing_value_range[0], missing_value_range[1], mini_sample_range[0],\
	    	     	mini_sample_range[1], avg_dict_range[0], avg_dict_range[1],\
	    	     	stdv_dict_range[0], stdv_dict_range[1]))
	else:

	    a_tab = a_tab.loc[a_tab['Length'] >= Length]
	    a_tab = a_tab.loc[a_tab['#SNV'] >= SNV]
	    a_tab = a_tab.loc[a_tab['SNV_density'] <= SNV_density]
	    a_tab = a_tab.loc[a_tab['Missing_value'] <= missing_value]
	    a_tab = a_tab.loc[a_tab['#Samples(missing information < 5%)'] >= mini_sample]
	    a_tab = a_tab.loc[a_tab['Avg_genetic_distance'] <= avg_dist]
	    a_tab = a_tab.loc[a_tab['Stdv_genetic_distance'] <= stdv_dist]
	    a_tab = a_tab.loc[(a_tab['correlation'] >= cor) | (a_tab['correlation'].isnull())]
	    a_tab = a_tab.loc[a_tab['avg_time_measured_dist'] <= tm_dist]

	    return a_tab  





def main():

    parser = argparse.ArgumentParser('Assessment', 'Select', 'Visual')
    subparsers = parser.add_subparsers(help = 'program mode', dest = 'mode')
    add_assessing_cmd_options(subparsers)
    add_select_cmd_options(subparsers)
    add_visual_cmd_options(subparsers)

    args = parser.parse_args()


    if args.mode == 'Assessment':
        wga_aln = AlignIO.read(args.wga, 'fasta')

        opt_dir = os.getcwd()+'/{}'.format(args.opt_dir)
        if os.path.exists(opt_dir):
            shutil.rmtree(opt_dir)
        os.makedirs(opt_dir)
        opt_file_name = open(opt_dir + '/Assessment.txt', 'w')
        opt_file_name.write('Child_alignment\t' + 'Function\t' + 'Length\t' + '#SNV\t' + 'SNV_density\t'\
        	+ 'Missing_value\t' + '#Samples(missing information < 5%)\t'+'Avg_genetic_distance\t'\
        	+ 'Stdv_genetic_distance' + '\t' + 'correlation'+'\t'+'avg_time_measured_dist'+'\n')
        
        ipt_mp_dict = {i.rstrip().split('\t')[0]: float(i.rstrip().split('\t')[1])\
        for i in open(args.time_file).readlines()}
        combined = itertools.combinations(ipt_mp_dict.keys(), 2)
        time_diff_lst = [(i, abs(ipt_mp_dict[i[0]] - ipt_mp_dict[i[1]])+1) \
        for i in combined]


        def Assessment_tab(gff3):

            gContent_db = parse_gff(wga_aln, gff3)
            nCDS_length = [len(gContent_db[i].Child_aln[1,:]) for i in gContent_db if i.startswith('nCDS')]
            coding_percent = str(100*(1 - sum(nCDS_length)/len(wga_aln[1,:])))+'%'

            label_aln_pairs = [(i, gContent_db[i].Child_aln) for i in gContent_db]
            len_dict = {a: ChildAlnStats(gContent_db[a].Child_aln).aln_len() for a in gContent_db}
            vs_dict = multi_proc_dict(args.number_of_processors, multi_variable_sites, label_aln_pairs)
            mv_dict = {a: ChildAlnStats(gContent_db[a].Child_aln).missing_value() for a in gContent_db}  
            no_samples_min_mv_dict = {a: ChildAlnStats(gContent_db[a].Child_aln).no_sample_minimum_samples(args.minimum_missing_information)\
             for a in gContent_db}

            avg_pwd_dict = multi_proc_dict(args.number_of_processors, multi_dist, label_aln_pairs)
                # {'label': (mean, stdv)}
            	
            label_aln_td_lst = [(i, gContent_db[i].Child_aln, time_diff_lst) for i in gContent_db]
            avg_tm_pwd_dict = multi_proc_dict(args.number_of_processors, multi_tm_dist, label_aln_td_lst)





            dict_lst_BigMatrix = defaultdict(list)

            
            for label in gContent_db:
                if len_dict[label] >= args.minimum_length: 
                    dict_lst_BigMatrix[label].append([label, gContent_db[label].Function, len_dict[label],\
                     dict_judge(label, vs_dict)[0], dict_judge(label, vs_dict)[1], \
                        dict_judge(label, mv_dict), dict_judge(label, no_samples_min_mv_dict),\
                        dict_judge(label, avg_pwd_dict)[0], dict_judge(label, avg_pwd_dict)[1],\
                        dict_judge(label, avg_tm_pwd_dict)[0], dict_judge(label, avg_tm_pwd_dict)[1]])
            
            return dict_lst_BigMatrix

        if args.gff3_annotation:

            dict_lst_BigMatrix = Assessment_tab(args.gff3_annotation)
      
        else:

            opt_RefSeq = [i for i in wga_aln\
             if (not i.id.startswith('a__')) and (not i.id.startswith('m__'))][0]
            SeqIO.write(opt_RefSeq, opt_RefSeq.id, 'fasta')
            RefSeq_file = os.getcwd()+'/{}'.format(opt_RefSeq.id)
            prokka_opt_dir = os.getcwd()+'/{}'.format(opt_RefSeq.id) + '.prokka'
            prokka_file_prefix = opt_RefSeq.id
            cmd = 'prokka {} --outdir {} --prefix {}'.format(RefSeq_file, prokka_opt_dir, prokka_file_prefix)
            subprocess.call(cmd, shell = True)
            gff3_annotation_file = opt_RefSeq.id + '.prokka/' + prokka_file_prefix + '.gff'

            dict_lst_BigMatrix = Assessment_tab(gff3_annotation_file)
        
        for aln in dict_lst_BigMatrix:
            write_in_lines = "\t".join(conv_2_str(dict_lst_BigMatrix[aln][0]))
            opt_file_name.write(write_in_lines + '\n')
        opt_file_name.close()

       

    elif args.mode == 'Select':
        wga_aln = AlignIO.read(args.wga, 'fasta')
        tab_df = pd.read_csv(args.a, sep = '\t')
        pars = (args.length, args.variable_sites, args.variable_sites_density,\
        	args.missing_value, args.minimum_no_sample_controled_mv, args.average_pairwise_distance,\
        	args.average_pairwise_distance_stdv, args.minimum_corrlation, args.time_measured_avg_distance)
        selected_alns = aln_filter(tab_df, pars)

        if args.output_file:
        	opt_file_path = os.getcwd()+'/{}'.format(args.output_file)
        	selected_alns.to_csv(opt_file_path, index = False, header = True, sep = '\t')
        else:
        	sys.stdout.write("\t".join(selected_alns.columns)+'\n')
        	for i in selected_alns.index:
        		selected_row = selected_alns.loc[i,:]
        		
        		sys.stdout.write(selected_row["Child_alignment"]+'\t'\
        			+ str(selected_row['Function'])+'\t'\
        			+ selected_row['Length'].astype(str)+'\t'\
        			+ selected_row['#SNV'].astype(str)+'\t'\
        			+ selected_row['SNV_density'].astype(str)+'\t'\
        			+ selected_row['Missing_value'].astype(str)+'\t'\
        			+ selected_row['#Samples(missing information < 5%)'].astype(str)+'\t'\
        			+ selected_row["Avg_genetic_distance"].astype(str)+'\t'\
        			+ selected_row["Stdv_genetic_distance"].astype(str)+'\t'\
        			+ selected_row["correlation"].astype(str)+'\t'\
        			+ selected_row["avg_time_measured_dist"].astype(str) +'\n')

        if args.individual_alignments:
        	opt_dir = os.getcwd() + '/individual_alignments'
        	if os.path.exists(opt_dir):
        	    shutil.rmtree(opt_dir)
        	os.makedirs(opt_dir)
        	for i in selected_alns.index:
        		aln_header = selected_alns.loc[i,:]["Child_alignment"]
        		slicing_index = (int(aln_header.split('_')[1]), int(aln_header.split('_')[2]))
        		sub_aln = wga_aln[:,slicing_index[0]: slicing_index[1]]
        		SeqIO.write(sub_aln, opt_dir+'/'+aln_header+'.fna', 'fasta')   	
        else:
        	pass

        if args.concatenate:
        	opt_file = os.getcwd()+'/{}'.format(args.concatenate)
        	init_aln = wga_aln[:, 0:1]
        	for i in selected_alns.index:
        		aln_header = selected_alns.loc[i,:]["Child_alignment"]
        		slicing_index = (int(aln_header.split('_')[1]), int(aln_header.split('_')[2]))
        		sub_aln = wga_aln[:,slicing_index[0]: slicing_index[1]]
        		init_aln += sub_aln
        	init_aln = init_aln[:, 1:]
        	SeqIO.write(init_aln, opt_file, 'fasta')
        	build_raxml(args.concatenate, args.raxml_threads)
        	if args.raxml_threads != 0 and args.fast_temporal_signal_est:
        		tre = subprocess.getoutput('ls raxml_tree/RAxML_bipartitions.*')
        		opt_png = args.concatenate.replace('.fna', '') + '_TempEst.png'
        		cmd = 'PreClock_LRM.py {} {} {}'.format(tre, args.fast_temporal_signal_est, opt_png)
        		subprocess.call(cmd, shell = True)

        		

        else:
        	pass

    elif args.mode == 'Visual':
        print('In the development ......')
        if args.all:
        	df_ = pd.read_csv(args.at, sep = '\t')
        	v1 = df_[df_['avg_time_measured_dist'] > 0]['avg_time_measured_dist']

        	v2 = df_['correlation'].dropna()
        	fig, (ax1, ax2) = plt.subplots(2,1)
        	hist_plot(ax1, v1, {'bins': 200})
        	ax1.text(0.003, 30, 'Max: {} mutation/site/year\
        		\nMin: {} mutation/site/year'.format(str(v1.max()), str(v1.min())))
        	hist_plot(ax2, v2, {'bins': 50})
        	fig.savefig('test.png')



        else:
            pass




    else:
        sys.exit('Please choose modes procvided in the menu !')
if __name__ == '__main__':
    main()





