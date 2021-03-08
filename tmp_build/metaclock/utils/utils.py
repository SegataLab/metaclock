#!/usr/bin/env python
from __future__ import division
import re 
import os
import multiprocessing as mp
import numpy as np
import bz2
import sys
from Bio import SeqIO
from itertools import chain
import argparse
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from operator import itemgetter 
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

    

def pattern_filtering(p,l):
    """
    pattern_filtering is a function taking a defined pattern and a list as arguments.
    It filter out elements in a list which does not match the given pattern 
    """
    new_l = []
    for f in l:
        if re.match(p,f):
            new_l.append(f)
    return new_l        

def multi_map(w,f,l):
    """
    multiprocessor is a function taking an integer(number of processor), a defined function,
    and a list of works as arguments. 
    """
    # This is to multiprocess command line on shell
    pool = mp.Pool(processes = w)
    return pool.starmap(f,l)

def multi_proc_dict(w, f, dict_lst):
    
    pool = mp.Pool(processes = w)
    results = pool.map(f, dict_lst)

    return {i[0]: i[1] for i in results}



def unique(list_):
    """
    It generates a list which only contains unique elements
    """
    x = np.array(list_)
    return np.unique(x)

def find(s):
    """
    It finds all indexes of gaps in a sequence
    """

    return [i for i, ltr in enumerate(s) if ltr == '-']

def openr( fn, mode = 'r'):
    if fn is None:
        return sys.stdin
    return bz2.BZ2File(fn) if fn.endswith(".bz2") else open(fn, mode)


def openw( fn ):
    if fn is None:
        return sys.stdout
    return bz2.BZ2File(fn, 'w') if fn.endswith(".bz2") else open(fn, "w")

def is_numbers(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

class stats(object):

    def __init__(self, sequence):
        self.sequence = sequence

    def calc_GapRatio(self):
        length=len(self.sequence)
        gaps = self.sequence.count('N') + self.sequence.count('-')
        r = round(gaps/length,2)
        return r

    def calc_genome(self):
        genome1 = self.sequence.replace("N","")
        genome2 = genome1.replace("-","")
        return len(self.sequence),len(genome2)
    
    def calc_CG(self):
        GC = self.sequence.count("C") +self.sequence.count("c")+self.sequence.count("g")+ self.sequence.count("G")
        if GC != 0:
                return round(GC/len(self.sequence.replace("-","")),2)
        else:
                return "0"

def out_stats(fna_file, opt_stats= 'stats.tsv'):
    seqio_dict = SeqIO.to_dict(SeqIO.parse(fna_file, "fasta"))
    opt_file = opt_stats
    opt = open(opt_file,"w")
    opt.write("SeqHeader\t"+ "GapRatio\t" + "AlignLength\t" + "SeqLength\t" +"GCcontent\n")
    for seq in seqio_dict:
        Seq = str(seqio_dict[seq].seq)
        S = stats(Seq)
        line = str(seq) + "\t" + str(S.calc_GapRatio())+ "\t" + str(S.calc_genome()[0])+ "\t" + str(S.calc_genome()[1]) + "\t" + str(S.calc_CG()) + '\n'
        opt.write(line)
    opt.close()

def homo_site_mapper(senario):
    """
    It takes subject sequence and query sequence and their
    alignment start positions, and map the postions of homologous
    sites back to genome position. Zero-based
    """ 
    sseq, s_start, s_end, qseq, q_start, q_end = senario
    s_pos_lst = []
    q_pos_lst = []
    s_idx = s_start - 1
    q_idx = q_start - 1
    s_site_lst = []
    q_site_lst = []
    s_com_rev_flag = []
    q_com_rev_flag = []
    if ((s_end - s_start) > 0) & ((q_end - q_start) > 0):
        for idx in range(len(sseq)):
            s_site_lst.append(sseq[idx])
            q_site_lst.append(qseq[idx])
            s_com_rev_flag.append('NO')
            q_com_rev_flag.append('NO')
            if (sseq[idx] != '-') & (qseq[idx] != '-'):
                s_pos_lst.append(s_idx)
                q_pos_lst.append(q_idx)
                s_idx += 1
                q_idx += 1
            elif (sseq[idx] != '-') & (qseq[idx] == '-'):
                s_pos_lst.append(s_idx)
                q_pos_lst.append('-')               
                s_idx += 1
            elif (sseq[idx] == '-') & (qseq[idx] != '-'):
                s_pos_lst.append('-')
                q_pos_lst.append(q_idx)
                q_idx += 1
            else:
                s_pos_lst.append('-')
                q_pos_lst.append('-')
    
    elif ((s_end - s_start) > 0) & ((q_end - q_start) < 0):
        for idx in range(len(sseq)):
            s_site_lst.append(sseq[idx])
            q_site_lst.append(qseq[idx])
            s_com_rev_flag.append('NO')
            q_com_rev_flag.append('YES')            
            if (sseq[idx] != '-') & (qseq[idx] != '-'):
                s_pos_lst.append(s_idx)
                q_pos_lst.append(q_idx)
                s_idx += 1
                q_idx -= 1
            elif (sseq[idx] != '-') & (qseq[idx] == '-'):
                s_pos_lst.append(s_idx)
                q_pos_lst.append('-')
                s_idx += 1
            elif (sseq[idx] == '-') & (qseq[idx] != '-'):
                s_pos_lst.append('-')
                q_pos_lst.append(q_idx)
                q_idx -= 1
            else:
                s_pos_lst.append('-')
                q_pos_lst.append('-')

    elif ((s_end - s_start) < 0) & ((q_end - q_start) > 0):
        for idx in range(len(sseq)):
            s_site_lst.append(sseq[idx])
            q_site_lst.append(qseq[idx])
            s_com_rev_flag.append('YES')
            q_com_rev_flag.append('NO')
            if (sseq[idx] != '-') & (qseq[idx] != '-'):
                s_pos_lst.append(s_idx)
                q_pos_lst.append(q_idx)
                s_idx -= 1
                q_idx += 1
            elif (sseq[idx] != '-') & (qseq[idx] == '-'):
                s_pos_lst.append(s_idx)
                q_pos_lst.append('-')
                s_idx -= 1
            elif (sseq[idx] == '-') & (qseq[idx] != '-'):
                s_pos_lst.append('-')
                q_pos_lst.append(q_idx)
                q_idx += 1
            else:
                s_pos_lst.append('-')
                q_pos_lst.append('-')

    else:
        for idx in range(len(sseq)):
            s_site_lst.append(sseq[idx])
            q_site_lst.append(qseq[idx])
            s_com_rev_flag.append('YES')
            q_com_rev_flag.append('YES')            
            if (sseq[idx] != '-') & (qseq[idx] != '-'):
                s_pos_lst.append(s_idx)
                q_pos_lst.append(q_idx)
                s_idx -= 1
                q_idx -= 1
            elif (sseq[idx] != '-') & (qseq[idx] == '-'):
                s_pos_lst.append(s_idx)
                q_pos_lst.append('-')
                s_idx -= 1
            elif (sseq[idx] == '-') & (qseq[idx] != '-'):
                s_pos_lst.append('-')
                q_pos_lst.append(q_idx)
                q_idx -= 1
            else:
                s_pos_lst.append('-')
                q_pos_lst.append('-')           
                    
    return q_pos_lst, s_pos_lst, q_site_lst, s_site_lst, q_com_rev_flag, s_com_rev_flag

def ancient_sample_tailor(aln_dict, num_a = 1):
    core_col = []
    a_seqs = [list(aln_dict[a].seq) for a in aln_dict if a.startswith('a__')]   
    for c in range(len(a_seqs[0])):
        col_sites = [i[c] for i in a_seqs]
        aDNA_nucleotide = len([i for i in col_sites if i != '-'])
        if aDNA_nucleotide >= num_a:
            core_col.append(c)
        else:
            continue

    record_lst = []
    for g in aln_dict:
        seq_lst = list(aln_dict[g].seq)
        core_seq = itemgetter(*core_col)(seq_lst)
        seq="".join(core_seq)
        _id = g
        record_lst.append(SeqRecord(Seq(seq), id = _id, description = ''))
    return record_lst



def distribution(value_list):
    sns.set(style = 'white', color_codes = True)
    ax = sns.distplot(value_list, kde=False, bins = 100)
    ax.set(xlabel = 'Missing information (%)', ylabel = 'Number of columns')
    
    return plt
def Barplot(df):
    sns.set(style = 'white', color_codes = True)
    ax = sns.barplot(x = 'Cutoffs', y = 'N. columns', hue = 'Samples', data = df)

    return plt

def EvalCdf(sample, x):
    count = 0.0
    for value in sample:
        if value <= x:
            count += 1
    prob = count / len(sample)
    return prob

def hist_plot(ax, data, param_dict):

    out = ax.hist(data, **param_dict)

    return out


def draw_damage_pattern(G2A_files, C2T_files, opt_dir):
    substitution_type, postions, freq = [], [], []

    for f in G2A_files:
        opt_f = [l.rstrip().split('\t') for l in open(f).readlines()]
        for p in opt_f[1:]:
            substitution_type.append('G>A')
            postions.append(int(p[0]))
            freq.append(float(p[1]))

    for f in C2T_files:
        opt_f = [l.rstrip().split('\t') for l in open(f).readlines()]
        for p in opt_f[1:]:
            substitution_type.append('C>T')
            postions.append(int(p[0]))
            freq.append(float(p[1]))

    d = {'Substitution type': substitution_type, 'Positions': postions, 'Frequency': freq}
    df_ = pd.DataFrame(data=d)
    print(df_)
    
    fig, ax = plt.subplots()
    sns.relplot(x="Positions", y="Frequency", kind = "line", data = df_, col = 'Substitution type', ax = ax)
    plt.savefig(opt_dir + '/damage_pattern.png', dpi=300, bbox_inches = 'tight')








