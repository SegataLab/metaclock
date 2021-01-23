#!/usr/bin/env python

import os, sys
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, IUPAC, Gapped
from Bio.Seq import Seq
import itertools
import collections
import argparse
import multiprocessing as mp
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import seaborn as sns
import matplotlib.pyplot as plt

"""
This script is to calculate SNP rates between taxa in the alignment igonring sites containing gap.
"""

def read_args(args):

    parser = argparse.ArgumentParser()

    parser.add_argument('ipt_msa',
                        nargs = "?",
                        metavar = "input_msa",
                        help = "Input a multiple sequence alignment in fasta form.",
                        type = str)

    parser.add_argument('-nproc',
                        '--number_processors',
                        help = 'Specify the number of processors to use. default: [1]',
                        type = int,
                        default = 1)

    return vars(parser.parse_args())






def calc_transitions(args):

    seq1_header, seq2_header, seq1, seq2 = args

    SNP = 0
    no_gap_len = 0
    gap_len = 0
    for i in range(len(seq1)):
        
        SEQ1 = seq1[i].upper()
        SEQ2 = seq2[i].upper()

        if (SEQ1 != '-') and (SEQ2 != '-'):
            no_gap_len += 1
            if SEQ1 != SEQ2:
                SNP += 1
            else:
                SNP += 0
        else:
            gap_len += 1
    
    if no_gap_len != 0:
        SNP_rate = float(SNP/no_gap_len)
    else:
        sys.stdout.write('{} and {} are have zero sites identical, hence we assign rate of 1.\n'.format(seq1_header, seq2_header))
        SNP_rate = 1

    gap_ratio = float(gap_len/len(seq1))        


    return seq1_header, seq2_header, SNP_rate, gap_ratio


def nproc_calc_tran(seq_pairs, seq_dict, nproc = 1):

    packed_args = [(seq[0], seq[1], str(seq_dict[seq[0]].seq), str(seq_dict[seq[1]].seq)) for seq in seq_pairs]

    pool = mp.Pool(processes = nproc)
    results = pool.map(calc_transitions, packed_args)

    return results


def draw_SNV_rates(snv_rates, opt_dir = os.getcwd()):
    df_ = pd.DataFrame(snv_rates, columns = ['Taxa_x', 'Taxa_y', 'SNP_rate', 'gap_ratio'])
    df_["SNP_rate"] = pd.to_numeric(df_["SNP_rate"])
    df_combine = df_[df_.Taxa_x != df_.Taxa_y]
    df_combine.to_csv('pairwise_snv_rates.tsv', index = False, sep = '\t')
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()

    rates = df_.pivot("Taxa_x", "Taxa_y", "SNP_rate")
    x = df_combine['SNP_rate']
    fig1 = sns.clustermap(rates)
    ax = fig1.ax_heatmap
    ax.set_ylabel("")
    ax.set_xlabel("")
    fig1.savefig(opt_dir + '/snv_rates_heatmap.png', bbox = 'tight')

    sns.distplot(x, kde = True, bins = 100, norm_hist=True, ax = ax2)
    fig2.savefig(opt_dir + '/snv_rates_distribution.png', bbox = 'tight')



    
def plot_snv_rates_main(ipt_msa_file, nproc, output_dir = os.getcwd()):
    seqs_dict = SeqIO.to_dict(SeqIO.parse(open(ipt_msa_file), 'fasta'))
    seq_pairs =  itertools.product(seqs_dict.keys(), seqs_dict.keys())
    snv_rates = nproc_calc_tran(seq_pairs, seqs_dict, nproc)

    draw_SNV_rates(snv_rates, opt_dir = output_dir)

def main():
    pars = read_args(sys.argv)
    plot_snv_rates_main(pars['ipt_msa'], pars['number_processors'])

if __name__ == "__main__":
    main()
    







