#!/usr/bin/env python

import subprocess
import sys
import shutil
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from datetime import datetime

"""
This script is to automatically prepare gffs for Multiple RefSeqs merging mode.
"""
def time_now():

    now = datetime.now()
    current_time = now.strftime("[%H:%M:%S]")

    return current_time

def label_length_check(label):
    if len(label) >= 20:
        short_lab = 'Seq_00001'
    else:
        short_lab = label

    return short_lab

def extract_refseq(alns, inter):
    opt_aln_lst = []
    for aln in alns:
        rec_lst = []
        label = aln.split('/')[-1]
        aln_dict = SeqIO.to_dict(SeqIO.parse(aln,'fasta'))
        aln_id = label_length_check(aln_dict[label].id)
        aln_seq = aln_dict[label].seq
        aln_opt = inter + '/' + label
        opt_aln_lst.append(aln_opt)
        rec_lst.append(SeqRecord(Seq(str(aln_seq),\
         generic_dna), id = aln_id, description = ''))
        SeqIO.write(rec_lst, aln_opt, 'fasta')
        sys.stdout.write('{} Reference genome has been extracted for annotation: {}\n'.format(time_now(), aln_opt))

    return opt_aln_lst

def run_prokka(seqs):
    all_refseqs = []
    all_outdirs = []

    for refseq in seqs:
        all_refseqs.append(refseq)
        outdir = refseq + '.prokka'
        all_outdirs.append(outdir)
        prefix = refseq.split('/')[-1]

        cmd = "prokka {} --outdir {} --prefix {} {}".format(refseq, outdir, prefix, refseq)
        subprocess.call(cmd, shell = True)

    return {'prokka_dirs': all_outdirs, 'refseqs': all_refseqs} 


def create_folder(name):
    folder_name = abspath_finder(name) # determine the abs path
    if os.path.exists(folder_name):
        print('{} exists: pass!'.format(folder_name))
        pass
    else:
        os.makedirs(folder_name)
        print('Creating folder {}'.format(folder_name))
    return folder_name

def abspath_finder(file):
    # Convert wahtever path into abs path
    return os.path.abspath(file)



def auto_annotation(all_alns_fna, inter):
    gffs_folder = create_folder(inter + '/gffs')
    sys.stdout.write('{} Start extracting reference genome from alignments.\n'.format(time_now()))
    all_refs = extract_refseq(all_alns_fna, inter)

    prokka_opt = run_prokka(all_refs)['prokka_dirs']
    refseqs = run_prokka(all_refs)['refseqs']

    for prokka_folder in prokka_opt:
        cmd = 'cp {}/*gff {}/.'.format(prokka_folder, gffs_folder)
        cmd_2 = 'rm -rf {}'.format(prokka_folder)
        subprocess.call(cmd, shell = True)
        subprocess.call(cmd_2, shell = True)

    for refseq in refseqs:
        cmd_3 = 'rm -r {}'.format(refseq)
        subprocess.call(cmd_2, shell = True)
    
    sys.stdout.write('{} All annotation files are in: {}\n'.format(time_now(), gffs_folder))
    return gffs_folder


