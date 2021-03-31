#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
import itertools
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pandas as pd
import seaborn as sns
from .utils import utils


def read_args(args):

    parser = argparse.ArgumentParser()

    parser.add_argument('ipt_aln',
                        nargs = "?",
                        metavar = "input_alignment",
                        help = "Input a MSA alignment in fasta",
                        type = str)

    parser.add_argument('opt_name',
                        nargs = "?",
                        metavar = "output_name",
                        help = "Specify the output name for resulted figure",
                        type = str)
    parser.add_argument('-w',
                        '--window_size',
                        nargs = '?',
                        help = "Sliding window size. default:[200bp]",
                        type = int,
                        default = 200)
    parser.add_argument('-s',
                        '--step_size',
                        nargs = '?',
                        help = 'Moving steps of sliding window. default:[50bp]',
                        type = int,
                        default = 50)
    parser.add_argument('-l_dist',
                        '--label_distance',
                        nargs = '?',
                        help = 'The distance betwen text label and the plotting. default: [0.15]',
                        type = float,
                        default = 0.15)
    parser.add_argument('-opt_text',
                        '--output_text',
                        nargs = '?',
                        help = 'Output a statistics text.',
                        type = str,
                        default = None)


    return vars(parser.parse_args())




def GenomeBar(genome_dict, opt_png, f_name = None, min_v = 0, max_v = 1, l_dist = 0.1):

    fig = plt.figure(figsize=(15,20))   
    x_scale = np.linspace(0,1, 100)


    axprops = dict(xticks=[], yticks=[]) 
    barprops = dict(aspect='auto', cmap=plt.cm.Blues, interpolation='nearest',vmin = min_v, vmax = max_v)

    bot = 0.8
    plt.axes([0.85, 0.85, 0.1, 0.000001])
    plt.yticks([])
    plt.xticks([0, 0.25, 0.5, 0.75, 1])

    ax_scale = fig.add_axes([0.85, 0.85, 0.1, 0.01], **axprops)
    ax_scale.imshow(x_scale.reshape((1, -1)), **barprops)
    ax_scale.text(0.5, 1.2, 'Reconstruction score', ha = 'center', transform=ax_scale.transAxes)

    for g in genome_dict:
        header = g
        x = np.asarray(genome_dict[g])
        ax = fig.add_axes([0.25, bot, 0.7, 0.01], **axprops)
        ax.imshow(x.reshape(1,-1), **barprops)
        ax.text(-l_dist, l_dist, header, ha = 'left', transform=ax.transAxes)
        bot -= 0.01

    return plt.savefig(opt_png)

def color_scale(sites_lst=range(20)):
    x = np.asarray(sites_lst)
    axprops = dict(xticks=[], yticks=[]) 
    barprops = dict(aspect='auto', cmap=plt.cm.seismic, interpolation='nearest',vmin = 0, vmax = 20)
    fig=plt.figure(figsize=(8,3))
    plt.axes([0.05,0.1,0.1,0.00001])
    plt.yticks([])
    plt.xticks(np.arange(0, 25, step=5))

    ax2 = fig.add_axes([0.05, 0.15, 0.1, 0.07], **axprops) # Draw a rectangle 
    ax2.imshow(x.reshape((1,-1)), **barprops)

    
    return plt.savefig('color_scale.png')
 

def win_counting_mv(g_seq, win_size, step_size):
    mv_lst = []
    for i in range(0, len(g_seq), step_size):
        stretch_mv = g_seq[i: i+win_size].count('-')
        mv_ratio = 1- (stretch_mv/win_size)
        mv_lst.append(mv_ratio)
    return mv_lst 

        
def window_sliding_mv(fna, win_size, step_size):
    # This function calculates missing value along genome in sliding window manner
    # and return a dictionary in which header is key and the list of windowed values is value
    Seq_dict = SeqIO.to_dict(SeqIO.parse(fna, 'fasta'))
    mv_dict = {}
    for i in Seq_dict:
        seq = Seq_dict[i].seq
        mv_dict[i] = win_counting_mv(seq, win_size, step_size)
    return mv_dict

def out_stats(fna_file, output_file):
    seqio_dict = SeqIO.to_dict(SeqIO.parse(fna_file, "fasta"))

    opt = open(output_file,"w")
    opt.write("SeqHeader\t"+ "GapRatio\t" + "AlignmentLength\t" + "SequenceLength\t" +"GCcontent\n")
    for seq in seqio_dict:
        Seq = str(seqio_dict[seq].seq)
        S = utils.stats(Seq)
        line = str(seq) + "\t" + str(S.calc_GapRatio())+ "\t" + str(S.calc_genome()[0])+ "\t" + str(S.calc_genome()[1]) + "\t" + str(S.calc_CG()) + '\n'
        opt.write(line)
    opt.close()

def main():
    pars = read_args(sys.argv)
    windowed_mv = window_sliding_mv(pars['ipt_aln'], pars['window_size'], pars['step_size'])
    GenomeBar(windowed_mv, pars['opt_name'], l_dist = pars['label_distance'])
    if pars['output_text']:
        out_stats(pars['ipt_aln'], pars['output_text'])


if __name__ == '__main__':
    main()


