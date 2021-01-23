#!/usr/bin/env python

import sys
import math
import subprocess
import argparse
import os
import itertools
import shutil
from collections import defaultdict
from functools import partial
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import pandas as pd
from operator import itemgetter
from .utils import utils


def read_args(args):

    parser = argparse.ArgumentParser()

    parser.add_argument('ga_file',
                        nargs = "?",
                        metavar = "genome_alignment_reconstruction",
                        help = "Input a genome alignment for tailoring. [fasta] format",
                        type = str)

    parser.add_argument('tailored_ga_file',
                        nargs = "?",
                        metavar = "output",
                        help = "Specify output name for the tailored genome alignment in fasta format",
                        type = str)

    parser.add_argument('-a',
                        "--automated_tailoring",
                        help = "This option allows automated tailoring using trimAl for minimizing missing information in the alignment",
                        action = 'store_true')

    parser.add_argument('-s',
                        "--stats",
                        help = "Specify a file name for statistic report of tailored alignment.",
                        type = str)

    parser.add_argument('-tt',
                        "--target_tailoring",
                        help = "Input a list of samples to keep for alignment tailoring. The limit for missing info in each column can\
                         be given and it is delimited by comma (e.g. shortlist.txt,0.1), or column tailoring will be performed automatically using trimAl (e.g. list.txt).",
                        type = str)

    return vars(parser.parse_args())


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

def main()
    pars = read_args(sys.argv)
    wga_aln = AlignIO.read(pars['ga_file'], 'fasta')
    aln_obj = tailor(wga_aln)


    if pars['target_tailoring']:
        sys.stdout.write("Tailoring is being peformed on target samples.\n")

        sl_args = short_list_tailaring_pars_check(pars['target_tailoring'])
        short_list = [i.rstrip() for i in open(sl_args[0]).readlines()]
        c = sl_args[1]

        if c:
            sys.stdout.write("The maximum gaps allowed in each column is: {}\n".format(str(c)))
            SeqIO.write(aln_obj.short_list_tailor(short_list, c), pars['tailored_ga_file'] ,'fasta')
        else:
            sys.stdout.write("Column-wise gaps are being trimmed usng trimAl.")
            aln_obj.short_list_tailor(short_list, c)
            subprocess.call('mv TrimalGappyout_opt.fna {}'.format(pars['tailored_ga_file']), shell = True)
            subprocess.call('rm tmp_feed_to_trimal.fna', shell = True)
        if pars['stats']:
            utils.out_stats(pars['tailored_ga_file'], opt_stats = pars['stats'])

        sys.stdout.write("Thanks for using and welcome back!\n")
    
    elif pars['automated_tailoring']:
        sys.stdout.write("Automated tailoring using trimAl is being performed only on columns and no samples have been removed.\n")
        cmd = 'trimal -gappyout -in {} -out TrimalGappyout_opt.fna'.format(pars['ga_file'])
        subprocess.call(cmd, shell = True)
        subprocess.call('mv TrimalGappyout_opt.fna {}'.format(pars['tailored_ga_file']), shell = True)
        sys.stdout.write("Thanks for using and welcome back!\n")
        utils.out_stats(pars['tailored_ga_file'], opt_stats = pars['stats'])


    else:
        sys.exit('Please choose one of these tailor strategies: automated tailoring or target tailoring.')            

if __name__ == '__main__':
    main()


    








