#!/usr/bin/env python

"""
NAME: Python program BuildGenomeAln.py is to build genome alignment which includes
ancient samples, NCBI references, or modern metagenome samples, using either bl
ast or mapping approach.  
"""

__author__ = "Kun D. Huang"
__version__ = "0.5"
__date__ = "25.02.2020"

import sys
import math
import timeit
import subprocess
import argparse
import os
import itertools
import pysam
import shutil
from collections import defaultdict
from functools import partial
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from utils import multi_map
from utils import unique
from utils import find
from utils import out_stats
from utils import ancient_sample_tailor
from utils import distribution
from utils import Barplot
from AlignStats import AlignStats
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def add_blast_approach_cmd_options(subparsers):
    """
    Add command options for blast approach.
    """
    blast_parser = subparsers.add_parser(
        'contigs_based',help = 'Use contigs_based approach',
        description = 'Use contigs_based approach.')

    blast_parser.add_argument(
        'ancient_sample',
        nargs = '?',
        metavar = 'AncientSamples',
        help = 'input folder names contains ancient microbiome reads \
            separated by comma or input a text file each line indicates \
            a folder name. If you input "None" no ancient sample integrated.',
        type = str)

    blast_parser.add_argument(
        'contigs',
        nargs = '?',
        metavar = 'contigs_isorefs_mags',
        help = 'specify the folder name which contains all query genomes \
            (either RefSeqs or MAGs), in FASTA form, used for blast.',
        type = str)

    blast_parser.add_argument(
        'ref_fna',
        nargs = '?',
        metavar = 'RefSeq_fasta',
        help = 'input a fasta file of reference sequence.',
        type = str)

    blast_parser.add_argument(
        '-trim',
        '--trim_reads_end',
        help = 'Trim the reads before computing the consensus. A value \
            of 10:10 means that the first and last 10 positions of \
            each read will be ignored. Default: None',
        type = str,
        default = None)

    blast_parser.add_argument(
        '-minqual',
        '--minimum_quality',
        help = 'specify the minimum quality of aligned reads to keep. \
            default: [30]',
        type = str,
        default = '30')

    blast_parser.add_argument(
        '-minlen',
        '--minimum_length',
        help = 'specify the minimum length of aligned reads to keep. \
            default: [30nt]',
        default = '30')

    blast_parser.add_argument(
        '-maxsnps',
        '--max_snps',
        help = 'maximum edit distance on the alignment for a read to pass. \
            default: [0.03]',
        type = str,
        default = '0.03')

    blast_parser.add_argument(
        '-mincov',
        '--minimum_coverage',
        help = 'minimum position coverage to keep for reconstruction. \
            default: [3]',
        default = '3',
        type = str)
    blast_parser.add_argument(
        '-mindomall',
        '--minimum_dominant_allele',
        help = 'Cutoff for degree of allele dominance for a position to be considered \
        polymorphic when building consensus sequence. default: [0.8]',
        default = 0.8,
        type = float)

    blast_parser.add_argument(
        '-m',
        '--search_mode',
        help = 'specify search mode in bowtie2. e.g. -a or -k,3. \
        Look for bowtie2 manual for details. default: [*-k,3*]. \
        Attention: * sign should close parameter so as to \
            escape special character, for example, *-a*, *-k,5*.',
        type = str,
        default = '*-k,3*')    

    blast_parser.add_argument(
        '-aln_len_blast',
        '--alignment_length',
        help = 'minimum alignment length of blast hits. default: [500nt]',
        type = str,
        default = '500')

    blast_parser.add_argument(
        '-pid_blast',
        '--identity',
        help = 'minimum identity of alignment of blast hits. default: [95.0]',
        type = str,
        default = '95.0')
    
    blast_parser.add_argument(
        '-blast_t',
        '--blast_threads',
        help = 'threading number of blastn. default: [1]',
        type = str,
        default = '1')

    blast_parser.add_argument(
        '-t',
        '--threads',
        help = 'threading number of bowtie2. default: [1]',
        type = str,
        default = '1')

    blast_parser.add_argument(
        '-p',
        '--processor',
        help = 'multiple processors to be used. Suggest specify \
            equivalent number of ancient samples for the maximum \
            computation trade-off. default: [1]',
        type = int,
        default = 1)
    
    blast_parser.add_argument(
        '-tr',
        '--trimmed_reads',
        help = 'generate a file contains reads with possible \
            degenerated sites being trimmed',
        action = 'store_true')

    blast_parser.add_argument(
        '-raxml_t',
        '--raxml_threads',
        help = 'Number of threads to use in building raxml tree. \
            (raxml tree will be built if specified.)',
        type = str,
        default = 0)


def add_mapping_approach_cmd_options(subparsers):

    """
    Add command options for mapping approach.
    """

    mapping_parser = subparsers.add_parser('reads_based',
                                            help = 'Use reads_based approach.',
                                            description = 'Use reads_based approach.')

    mapping_parser.add_argument('ref_fna',
                                nargs = '?',
                                metavar = 'RefSeq_fasta',
                                help = 'input a fasta file of reference sequence',
                                type = str)

    mapping_parser.add_argument('ancient_sample',
                                nargs = '?',
                                metavar = 'AncientSamples',
                                help = 'input folder names containing metagenomic reads\
                                separated by comma or input a text file each line indicates\
                                a folder name',
                                type = str,
                                default = None)

    mapping_parser.add_argument('modern_sample',
                                nargs = '?',
                                metavar = 'ModernSamples',
                                help = 'input folder names containing metagenomic reads\
                                separated by comma or input a text file each line indicates\
                                a folder name',
                                type = str)

    mapping_parser.add_argument('-trim',
                                '--trim_reads_end',
                                help = 'Trim the reads before computing the consensus. A value\
                                of 10:10 means that the first and last 10 positions of\
                                each read will be ignored. Default: None',
                                type = str,
                                default = None)

    mapping_parser.add_argument('-a_minqual',
                                '--a_minimum_quality',
                                help = 'specify the minimum quality of aligned reads to keep. default: [30](for ancient samples)',
                                type = str,
                                default = '30')

    mapping_parser.add_argument('-a_minlen',
                                '--a_minimum_length',
                                help = 'specify the minimum length of aligned reads to keep. default: [30](for ancient samples)',
                                default = '30')

    mapping_parser.add_argument('-a_maxsnps',
                                '--a_max_snps',
                                help = 'maximum edit distance on the alignment for a read to pass. default: [0.03](for ancient samples)',
                                type = str,
                                default = '0.03')

    mapping_parser.add_argument('-a_mincov',
                                '--a_minimum_coverage',
                                help = 'minimum position coverage to keep for reconstruction. default: [3](for ancient samples)',
                                default = '3',
                                type = str)    

    mapping_parser.add_argument('-a_bt_t',
                                '--a_threads',
                                help = 'threading number of bowtie2 set for ancient samples. default: [1](for ancient samples)',
                                default = '1',
                                type = str)

    mapping_parser.add_argument('-a_m',
                                '--a_search_mode',
                                help =  'specify search mode in bowtie2. e.g. -a or -k,3.\
                                Look for bowtie2 manual for details. default: [*-k,3*]. (for ancient samples)\
                                Attention: * sign should close parameter so as to escape special character\
                                , for example, *-a*, *-k,5*.',
                                type = str,
                                default = '-k,3')
    mapping_parser.add_argument('-m_minqual',
                                '--m_minimum_quality',
                                help = 'specify the minimum quality of aligned reads to keep. default: [30](for modern samples)',
                                type = str,
                                default = '30')

    mapping_parser.add_argument('-m_minlen',
                                '--m_minimum_length',
                                help = 'specify the minimum length of aligned reads to keep. default: [50](for modern samples)',
                                default = '50')

    mapping_parser.add_argument('-m_maxsnps',
                                '--m_max_snps',
                                help = 'maximum edit distance on the alignment for a read to pass. default: [0.03](for modern samples)',
                                type = str,
                                default = '0.03')

    mapping_parser.add_argument('-m_mincov',
                                '--m_minimum_coverage',
                                help = 'minimum position coverage to keep for reconstruction. default: [5](for modern samples)',
                                default = '5',
                                type = str)    

    mapping_parser.add_argument('-m_bt_t',
                                '--m_threads',
                                help = 'threading number of bowtie2 set for modern samples. default: [1](for modern samples)',
                                default = '1',
                                type = str)

    mapping_parser.add_argument('-m_m',
                                '--m_search_mode',
                                help = 'specify search mode in bowtie2. e.g. -a or -k,3.\
                                Look for bowtie2 manual for details. default: [*-k,3*]. (for modern samples)\
                                Attention: * sign should close parameter so as to escape special character\
                                , for example, *-a*, *-k,5*.',
                                type = str,
                                default = '-k,3')    

    mapping_parser.add_argument('-a_p',
                                '--a_processor',
                                help = 'multiple processors used for ancient samples. *the production of\
                                 threads and processor should be less than available CPUs. default: [1]',
                                default = 1,
                                type = int)

    mapping_parser.add_argument('-m_p',
                                '--m_processor',
                                help = 'multiple processors used for modern samples. *the production of\
                                 threads and processor should be less than available CPUs. default: [1]',
                                default = 1,
                                type = int)

    mapping_parser.add_argument('-tr',
                                '--trimmed_reads',
                                help = 'generate a file contains reads with possible degenerated sites being trimmed',
                                action = 'store_true')

    mapping_parser.add_argument('-raxml_t',
                                '--raxml_threads',
                                help = 'Number of threads to use in building raxml tree. (raxml tree will be built\
                                if specified)',
                                type = str,
                                default = 0)  


# def add_alignment_assessing_cmd_option(subparsers):
#   """
#   Add command options for alignment assessing. 
#   """

#   alignment_assessing_parser = subparsers.add_parser('alignment_assessing',
#                                           help = 'Assessment of whole genome alignment.',
#                                           description = 'Features of assessing alignment composition.')
#   alignment_assessing_parser.add_argument('-opt_dir',
#                                           '--output_directory',
#                                           help = 'Specify the output directory. Current working directory as default.',
#                                           type = str,
#                                           default = None)
#   alignment_assessing_parser.add_argument('genome_alignment',
#                                           nargs = '?',
#                                           help = 'Input the genome alignment in fasta file.',
#                                           metavar = 'GenomeAlignment',
#                                           type = str)
#   alignment_assessing_parser.add_argument('-cov',
#                                           '--coverage',
#                                           help = 'Assessing coverage of RefSeq in each sample.',
#                                           action = 'store_true')
#   alignment_assessing_parser.add_argument('-misinfo_col_dist',
#                                           '--missinginfo_column_distribution',
#                                           help = 'Plot the distribution of columns with missing information.',
#                                           action = 'store_true')
#   alignment_assessing_parser.add_argument('-vs',
#                                           '--variant_sites',
#                                           help = 'Report the absolute number and proportion of variant sites (biallelic and multiallelic).',
#                                           action = 'store_true')
#   alignment_assessing_parser.add_argument('-vs_ws',
#                                           '--variant_sites_window_sliding',
#                                           help = 'Specify the window size to slice whole alignment to calculate variant sites (biallelic and multiallelic).',
#                                           type = int,
#                                           default = 0)
#   alignment_assessing_parser.add_argument('-opt_hvf',
#                                           '--output_hyperviriable_regione_free_file',
#                                           help = 'specify the threshold to remove hypervariable regions.',
#                                           type = int,
#                                           default = 0)

  
def make_blast_db(ref_fna, opt_dir = os.getcwd()):

    """
    It creates database of reference genome for blast
    Input: reference genome in fasta format
    Program: call 'makeblastdb'
    Output: blast database files
    """

    cmd = 'makeblastdb -in {} -dbtype nucl -out {}'.format(opt_dir + '/' + ref_fna, opt_dir + '/' + ref_fna)
    subprocess.call(cmd, shell = True)

def generate_query_set_and_mf_file(contigs_folder):
    
    """
    Input: 'contigs' -> a folder of fasta files, each represents a genome
    Program: 
    1) It merges all individual fasta files into one fasta file, called 'INTER_QuerySet.fna',as an intermediate.
    2) It creates a mapping file delimited by 'comma', 1st column is genome name (i.e. file name) and 2nd column is 
    name of the contig which belongs to the genome.
    Output: 1) A merged fasta file, INTER_QuerySet.fna, for downstream analysis as input
            2) A mapping file, 1st column is genome name and 2nd column is contig name 
    """

    opt_file = open('INTER_genome_contigs_mp.csv', 'w')
    subprocess.call('cat {}/* > INTER_QuerySet.fna'.format(contigs_folder), shell = True)
    all_fna = subprocess.getoutput('ls {}/*'.format(contigs_folder)).split('\n')

    for f in all_fna:
        opened_fna = SeqIO.to_dict(SeqIO.parse(open(f), 'fasta'))
        for h in opened_fna:
            row = f.split('/')[1] + ',' + h + '\n'
            opt_file.write(row)
    opt_file.close()


def blast_(ref_fna, query, b_t, opt_dir = os.getcwd()):

    """
    It blasts all individual genome sequences which are merged in 'INTER_QuerySet.fna'
    against the reference genome which is preprocessed as blast database format

    Input: 1) reference genome name, 2) 'INTER_QuertSet.fna', 3) number of threads
    Program: Call 'blastn' with fixed command line
    Output: 'INTER_blast_opt_tmp.tab' a big table which contains all hits searched using blast
    """

    outfmt = '6 qaccver saccver pident length mismatch gapopen qstart qend\
 sstart send evalue bitscore qseq sseq'
    cmd = "blastn -db {} -query {} -outfmt '{}' -num_threads {}\
 -word_size 9 -out {}/INTER_blast_opt_tmp.tab".format(opt_dir+'/'+ref_fna, opt_dir+'/'+query, outfmt, b_t, opt_dir)
    subprocess.call(cmd, shell = True)

def QC_on_blastn(blast_tab, length, pid, opt_dir = os.getcwd()):

    """
    It applies QC to blast output
    Input: 1)'INTER_blast_opt_tmp.tab', 2) length of hits, 3) identity percentage of hits
    Program: call script 'bo6_screen.py'
    Output: filtered blast output, 'INTER_blastn_opt_tmp_cleaned.tab' 
    """

    cmd = 'cut -f 1-14 {} | bo6_screen.py --length {} --pid {} > {}/INTER_blastn_opt_tmp_cleaned.tab'.format(blast_tab, length, pid, opt_dir)
    subprocess.call(cmd, shell = True)

def bowtie2_build(ref_fna):

    cmd = 'bowtie2-build {} {}'.format(ref_fna, ref_fna)
    subprocess.call(cmd, shell = True)

def detect_reads_suffix(reads_foler):

    reads_file = subprocess.getoutput('ls {}/*fastq*'.format(reads_foler)).split('\n')[0]
    return reads_file.split('.')[-1]  

def mapping_(multi_args):

    ref, s, t, k = multi_args
    k = k.replace("*", "")
    suffix = detect_reads_suffix(s)
    if suffix == "bz2":
        cmd = 'bzcat {}/*fastq.bz2 | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | \
samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), s+'.bam')
        subprocess.call(cmd, shell = True)
    elif suffix == "gz":
        cmd = 'zcat {}/*fastq.gz | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | \
samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), s+'.bam')
        subprocess.call(cmd, shell = True)
    elif suffix == "fastq":
        cmd = 'cat {}/*fastq | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | \
samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), s+'.bam')
        subprocess.call(cmd, shell = True)
    else:
        sys.exit("Reads have to be in the form of .bz2, .gz or .fastq !")

def mapping_a(multi_args):

    ref, s, t, k = multi_args
    k = k.replace("*", "")
    suffix = detect_reads_suffix(s)
    labeled_ancient_opt_bam = '/'.join(s.split('/')[:-1])+'/a__'+ s.split('/')[-1] + '.bam'
    if suffix == "bz2":
        cmd = 'bzcat {}/*fastq.bz2 | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | \
samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), labeled_ancient_opt_bam)
        subprocess.call(cmd, shell = True)
    elif suffix == "gz":
        cmd = 'zcat {}/*fastq.gz | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | \
samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), labeled_ancient_opt_bam)
        subprocess.call(cmd, shell = True)
    elif suffix == "fastq":
        cmd = 'cat {}/*fastq | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | \
samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), labeled_ancient_opt_bam)
        subprocess.call(cmd, shell = True)
    else:
        sys.exit("Reads have to be in the form of .bz2, .gz or .fastq!")


def mapping_m(multi_args):

    ref, s, t, k = multi_args
    k = k.replace("*", "")
    suffix = detect_reads_suffix(s)
    labeled_modern_opt_bam = '/'.join(s.split('/')[:-1])+'/m__'+ s.split('/')[-1] + '.bam'
    if suffix == "bz2":
        cmd = 'bzcat {}/*.fastq.bz2 | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | \
samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), labeled_modern_opt_bam)
        subprocess.call(cmd, shell = True)
    elif suffix == "fastq.gz":
        cmd = 'zcat {}/*fastq.gz | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | \
samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), labeled_modern_opt_bam)
        subprocess.call(cmd, shell = True)
    elif suffix == "fastq":
        cmd = 'cat {}/*fastq | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | \
samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), labeled_modern_opt_bam)
        subprocess.call(cmd, shell = True)
    else:
        sys.exit("Reads have to be in the form of .bz2, .gz or .fastq!")

def create_mapping(args):

    if ',' in args.ancient_sample:
        sample_list = [os.getcwd() + '/' +  i for i in args.ancient_sample.split(',')]
    elif '.txt' in args.ancient_sample:
        sample_list = [os.getcwd() + '/' +  i.rstrip() for i in open(args.ancient_sample).readlines()]
    else:
        sample_list = [os.getcwd() + '/' + args.ancient_sample]

    threads = [args.threads] * len(sample_list)
    ref_lst = [args.ref_fna] * len(sample_list)
    mode = [args.search_mode] * len(sample_list)

    multi_map(args.processor, mapping_, zip(ref_lst, sample_list, threads, mode))   

def create_mapping_a(args):

    '''
    Mapping procedure specific for ancient samples
    '''
    
    if ',' in args.ancient_sample:
        sample_list = [os.getcwd() + '/' + i for i in args.ancient_sample.split(',')]
    elif '.txt' in args.ancient_sample:
        sample_list = [os.getcwd() + '/' + i.rstrip() for i in open(args.ancient_sample).readlines()]
    else:
        sample_list = [os.getcwd() + '/' + args.ancient_sample]
    
    threads = [args.a_threads] * len(sample_list)
    ref_lst = [args.ref_fna] * len(sample_list)
    mode = [args.a_search_mode] * len(sample_list)

    multi_map(args.a_processor, mapping_a, zip(ref_lst, sample_list, threads, mode))    

def create_mapping_m(args):

    if ',' in args.modern_sample:
        sample_list = [os.getcwd() + '/' + i for i in args.modern_sample.split(',')]
    elif '.txt' in args.modern_sample:
        sample_list = [os.getcwd() + '/' +  i.rstrip() for i in open(args.modern_sample).readlines()]
    else:
        sample_list = [os.getcwd() + '/' + args.modern_sample]

    threads = [args.m_threads] * len(sample_list)
    ref_lst = [args.ref_fna] * len(sample_list)
    mode = [args.m_search_mode] * len(sample_list)

    multi_map(args.m_processor, mapping_m, zip(ref_lst, sample_list, threads, mode))    


def bam_filter(multi_args):

    bam, mq, ml, m_snp = multi_args
    filter_path = os.path.dirname(__file__)+"/cmseq/filter.py"
    cmd = "samtools view -h -F 0x4 {} | python3 {} --minqual {} --minlen {} --maxsnps {} > {}".format(bam, filter_path,\
     mq, ml, m_snp, "filtered_INTER_"+bam)
    subprocess.call(cmd, shell =True)


def parallel_filter(args):

    bams = subprocess.getoutput("ls *.bam").split('\n')
    mq_lst = [args.minimum_quality] * len(bams)
    ml_lst = [args.minimum_length] * len(bams)
    m_snp = [args.max_snps] * len(bams)

    multi_map(args.processor, bam_filter, zip(bams, mq_lst, ml_lst, m_snp))

def parallel_filter_a(args):

    bams = subprocess.getoutput("ls a__*.bam").split('\n')
    mq_lst = [args.a_minimum_quality] * len(bams)
    ml_lst = [args.a_minimum_length] * len(bams)
    m_snp = [args.a_max_snps] * len(bams)

    multi_map(args.a_processor, bam_filter, zip(bams, mq_lst, ml_lst, m_snp))

def parallel_filter_m(args):

    bams = subprocess.getoutput("ls m__*.bam").split('\n')
    mq_lst = [args.m_minimum_quality] * len(bams)
    ml_lst = [args.m_minimum_length] * len(bams)
    m_snp = [args.m_max_snps] * len(bams)

    multi_map(args.m_processor, bam_filter, zip(bams, mq_lst, ml_lst, m_snp))


def build_consensus(multi_args):

    f_bam, mc, trim_end, domall = multi_args
    consensus_path = os.path.dirname(__file__)+"/cmseq/consensus.py"

    if trim_end != None:

        cmd = "python3 {} {} --sortindex --mincov {} --trim {} --dominant_frq_thrsh {} > {}".format(consensus_path, f_bam, mc, trim_end, domall, f_bam.split(".")[0]+'_consensus.fna')
        subprocess.call(cmd, shell = True)
    else:
        cmd = "python3 {} {} --sortindex --mincov {} --dominant_frq_thrsh {} > {}".format(consensus_path, f_bam, mc, domall, f_bam.split(".")[0]+'_consensus.fna')
        subprocess.call(cmd, shell = True)

def parallel_consensus(args):

    f_bams = subprocess.getoutput("ls filtered_INTER_*bam").split('\n')
    mc_lst = [args.minimum_coverage] * len(f_bams)
    trim_lst = [args.trim_reads_end] * len(f_bams)
    domall_lst = [args.minimum_dominant_allele] * len(f_bams)
    multi_map(args.processor, build_consensus, zip(f_bams, mc_lst, trim_lst, domall_lst))


def parallel_consensus_a(args):

    f_bams = subprocess.getoutput("ls filtered_INTER_a__*bam").split('\n')
    mc_lst = [args.a_minimum_coverage] * len(f_bams)
    trim_lst = [args.trim_reads_end] * len(f_bams) 
    multi_map(args.a_processor, build_consensus, zip(f_bams, mc_lst, trim_lst))

def parallel_consensus_m(args):

    f_bams = subprocess.getoutput("ls filtered_INTER_m__*bam").split('\n')
    mc_lst = [args.m_minimum_coverage] * len(f_bams)
    trim_lst = [args.trim_reads_end] * len(f_bams) 
    multi_map(args.m_processor, build_consensus, zip(f_bams, mc_lst, trim_lst))


class blastn_sort(object):
    """
    object blastn_sort is to organize information for each hit of blast in form 6 table
    """

    def __init__(self, hit):

        self.hit = hit
        self.qseqid = hit[0]
        self.sseqid = hit[1]
        self.pident = hit[2]
        self.length = hit[3]
        self.mismatch = hit[4]
        self.gapopen = hit[5]
        self.qstart = int(hit[6])
        self.qend = int(hit[7])
        self.sstart = int(hit[8])
        self.send = int(hit[9])
        self.evalue = hit[10]
        self.bitscore = float(hit[11])
        self.qseq = hit[12]
        self.sseq = hit[13]

def get_homo_query(homo_sub, homo_query):

    """
    This function takes sub_query alignment as arguments. To make sure reconst
    ructed homologous sequence from query has same length as reference sequence
    used for blast and mapping, it cuts gappy postion in subject sequence and s
    ame postion in query sequence as well. It returns modified query sequence.I
    f subject sequence is reversed, it will reverse query sequence.
    """

    gappy_idx = find(homo_sub) # find gappy positions in subject sequence
    query_lst = list(homo_query)
    for i in gappy_idx:
        query_lst[i] = '$'
    return "".join(query_lst).replace("$","")

def sort_contig(contig_coor_dir_bit_seq):
    
    """
    It sorts all hits of all contigs from one genome, based on start position.
    If hits share same start position, then consider bitscore.
    """
    sorted_dict = {}
    for c in contig_coor_dir_bit_seq:
        lst = sorted(contig_coor_dir_bit_seq[c], key = lambda x : (x[0], x[-1]))
        sorted_dict[c] = lst
    return sorted_dict

def bricklayer(contig, range_lst):

    """
    bricklayer is greedy algorithm which minimizes gappy sites.
    1) Create a contig backbone filled in with '-' using real ref contig length
    2) loop through the list of ranges
    3) replace '-' in backnone with homo query sequence form blas, using coord
    inates on subject
    4) Creat a nest of conditional states
        1> If next range is completely located within 'layed' region, continue
        2> If next range partially overlapped with 'layed' region, expand inte
        rected region
        3> replace '-' with range  
    """
    contig_lst = ['-']*len(contig) # Create a backnone using ref contig length in list
    range_start = range_lst[0][0] - 1
    range_end = range_lst[0][1]
    contig_lst[range_start : range_end] = range_lst[0][-1]
    layed_part = [[range_start, range_end]]

    for r in range(1, len(range_lst)):
        if (range_lst[r][0]-1 >= layed_part[-1][1]):
            contig_lst[range_lst[r][0] - 1 : range_lst[r][1]] = range_lst[r][-1]
            layed_part.append([range_lst[r][0] - 1, range_lst[r][1]])

        elif (range_lst[r][0] - 1 < layed_part[-1][1]) & (range_lst[r][1] <= layed_part[-1][1]):
            continue

        elif (range_lst[r][0] - 1 < layed_part[-1][1]) & (range_lst[r][1] >= layed_part[-1][1]):
            distance = int(range_lst[r][1] - layed_part[-1][1])
            contig_lst[layed_part[-1][1] : range_lst[r][1]] = range_lst[r][-1][-distance:]
            layed_part.append([layed_part[-1][1], range_lst[r][1]])

    return ''.join(contig_lst)  

def complement_seq(seq):
    mydna = Seq(seq, generic_dna)
    return str(mydna.complement())

def reconstruct(genome_contig_dir, query_genome_name, ref_dict):
    contig_coor_dir = defaultdict(list)
    for f in genome_contig_dir[query_genome_name]:
        contig_coor_dir[f[1]].append([f[2], f[3], f[4]])

    record_dict = {}
    for ctig in contig_coor_dir:
        x = sort_contig(contig_coor_dir)[ctig]
        seq = str(ref_dict[ctig].seq)
        recon = bricklayer(seq, x)
        record_dict[str(ctig)] = recon

    return record_dict


def concatenate_contigs(blast_tab, mp_file, ref_fna):
    ipt_mtx = [l.rstrip('\n').split('\t') for l in open(blast_tab).readlines()]
    mp_dict = {l.rstrip().split(',')[1] : l.rstrip().split(',')[0] for l in open(mp_file).readlines()}
    query_genome = unique([mp_dict[i[0]] for i in ipt_mtx])
    genome_contig_dir = defaultdict(list)

    for hit in ipt_mtx:
        hit_object = blastn_sort(hit)
        new_homo_query = get_homo_query(hit_object.sseq, hit_object.qseq)

        if hit_object.send - hit_object.sstart > 0:
            genome_contig_dir[mp_dict[hit_object.qseqid]].append((hit_object.qseqid, hit_object.sseqid,
                hit_object.sstart, hit_object.send, new_homo_query, hit_object.bitscore))
        else:
            genome_contig_dir[mp_dict[hit_object.qseqid]].append((hit_object.qseqid, hit_object.sseqid,
                hit_object.send, hit_object.sstart, complement_seq(new_homo_query[::-1]),hit_object.bitscore))

    ref_dict = SeqIO.to_dict(SeqIO.parse(open(ref_fna), 'fasta'))
    ref_len_dir = {c : len(ref_dict[c].seq) for c in ref_dict}
    genome_name_lst = list(set([mp_dict[q_g] for q_g in mp_dict]))
    genome_dict_ = {g : reconstruct(genome_contig_dir, g, ref_dict) for g in genome_name_lst}

    return genome_dict_

def reorder_contigs(ref_fna_dict, recon_genome_dict):
    """
    Works on single genome reordering !
    """
    ref_ctig_coordinate=sorted([(i, len(ref_fna_dict[i].seq)) for i in ref_fna_dict], key=lambda x: x[0])

    init=""
    for c in ref_ctig_coordinate:
        if c[0] in recon_genome_dict:
            init += recon_genome_dict[c[0]]
        else:
            init += len(ref_fna_dict[c[0]]) * "-"
    return init

def write_one_file(genomes_dict_, ref_fna):
    ref_ctig_dict = SeqIO.to_dict(SeqIO.parse(open(ref_fna), "fasta"))
    recon_genomes_list = []
    for g in genomes_dict_:
        # print(g)
        recon_seq = reorder_contigs(ref_ctig_dict, genomes_dict_[g])
        recon_genomes_list.append(SeqRecord(Seq(recon_seq, generic_dna), id = 'm__'+ g, description = '')) 

    return recon_genomes_list       

def output_trimmed_reads(trim_pos, bam_file):

    in_samfile = pysam.AlignmentFile(bam_file, 'rb')
    bam_file_opt = 'trimmed_'+'_'.join(bam_file.split('.')[0].split('_')[2:])+'.fastq'
    out_fastq = open(bam_file_opt, 'w')
    left, right = trim_pos.split(':')
    for read in in_samfile.fetch():
        if read.is_reverse:
            read_name = read.query_name
            read_seq = str(Seq(read.query_alignment_sequence[int(left): -int(right)], generic_dna)\
                .reverse_complement())
            read_qual = read.qual[int(left): -int(right)][::-1]

            out_fastq.write('@'+read_name+'\n')
            out_fastq.write(read_seq+'\n')
            out_fastq.write("+\n")
            out_fastq.write(read_qual+'\n')
        else:
            read_name = read.query_name
            read_seq = str(Seq(read.query_alignment_sequence[int(left): -int(right)], generic_dna))
            read_qual = read.qual[int(left): -int(right)]

            out_fastq.write('@'+read_name+'\n')
            out_fastq.write(read_seq+'\n')
            out_fastq.write("+\n")
            out_fastq.write(read_qual+'\n')
    out_fastq.close()           


def generate_par_report(args):
	"""
	This collects arguments of main parameters and outputs them in 'Parameters_setting.txt' 
	"""

	rep_opt = open('Parameters_setting.txt', 'w')
	if args.mode == 'contigs_based':
		blast_CPU = int(args.blast_threads)
		bowtie2_CPU = int(args.threads)*int(args.processor)
		raxml_CPU = int(args.raxml_threads)
		peak_cpu = max([blast_CPU, bowtie2_CPU, raxml_CPU])
		rep_opt.write('Alignment is generated in mode of [{}]\n'.format(args.mode))
		if args.ancient_sample == 'None':
			rep_opt.write('No ancient sample integrated!\n')
		else:
			pass
		if args.minimum_quality is not None and args.ancient_sample != 'None':
			rep_opt.write('Minimum quality of aligned reads kept [{}]\n'.format(args.minimum_quality))
		else:
			pass
		if args.minimum_length is not None and args.ancient_sample != 'None':
			rep_opt.write('Minimum length of aligned reads kept [{}nt]\n'.format(args.minimum_length))
		else:
			pass
		if args.minimum_coverage is not None and args.ancient_sample != 'None':		
			rep_opt.write('Minimum coverage of postions kept [{}]\n'.format(args.minimum_coverage))
		else:
			pass
		if args.max_snps is not None and args.ancient_sample != 'None':
			rep_opt.write('Maxium edit distance on the alignment for a read to pass [{}]\n'.format(args.max_snps))
		else:
			pass
		if args.identity is not None:
			rep_opt.write('Minimum alignment identity in blastn [{}]\n'.format(args.identity))
		else:
			pass
		if args.trim_reads_end is not None and args.ancient_sample != 'None':
			rep_opt.write('Length of reads end trimmed: {}\n'.format(args.trim_reads_end))
		else:
			pass
		if args.trimmed_reads and args.ancient_sample != 'None':
			rep_opt.write('Reads trimmed for damaged sites and used for reconstruction is output in fastq!\n')
		else:
			pass
		rep_opt.write('Minimum alignemnt length in blastn [{}]\n'.format(args.alignment_length))
		rep_opt.write('Peak CPU used [{}]\n'.format(str(peak_cpu)))
		rep_opt.close()

	if args.mode == 'reads_based':
		aDNA_CPU = int(args.a_threads)*int(args.a_processor)
		mDNA_CPU = int(args.m_threads)*int(args.m_processor)
		raxml_CPU = int(args.raxml_threads)
		peak_cpu = max([aDNA_CPU, mDNA_CPU, raxml_CPU])
		rep_opt.write('Alignment is generated in mode of [{}]\n'.format(args.mode))
		if args.ancient_sample == 'None':
			rep_opt.write('No ancient sample integrated!\n')
		else:
			pass
		if args.a_minimum_quality is not None and args.ancient_sample != 'None':
			rep_opt.write('Minimum quality of aligned reads (ancient samples) kept [{}]\n'.format(args.a_minimum_quality))
		else:
			pass
		if args.a_minimum_length is not None and args.ancient_sample != 'None':
			rep_opt.write('Minimum length of aligned reads (ancient samples) kept [{}nt]\n'.format(args.a_minimum_length))
		else:
			pass
		if args.a_minimum_coverage is not None and args.ancient_sample != 'None':		
			rep_opt.write('Minimum coverage (ancient samples) of postions kept [{}]\n'.format(args.a_minimum_coverage))
		else:
			pass
		if args.a_max_snps is not None and args.ancient_sample != 'None':
			rep_opt.write('Maxium edit distance on the alignment for a read (ancient samples) to pass [{}]\n'.format(args.a_max_snps))
		else:
			pass
		if args.m_minimum_length is not None:
			rep_opt.write('Minimum length of aligned reads (modern samples) kept [{}nt]\n'.format(args.m_minimum_length))
		else:
			pass
		if args.a_minimum_coverage is not None:		
			rep_opt.write('Minimum coverage (modern samples) of postions kept [{}]\n'.format(args.m_minimum_coverage))
		else:
			pass
		if args.a_max_snps is not None:
			rep_opt.write('Maxium edit distance on the alignment for a read (modern samples) to pass [{}]\n'.format(args.m_max_snps))
		else:
			pass
		if args.trim_reads_end is not None and args.ancient_sample != 'None':
			rep_opt.write('Length of reads end trimmed: {}\n'.format(args.trim_reads_end))
		else:
			pass
		if args.trimmed_reads and args.ancient_sample != 'None':
			rep_opt.write('Reads trimmed for damaged sites and used for reconstruction is output in fastq!\n')
		else:
			pass
		rep_opt.write('Peak CPU utilization: [{}]\n'.format(str(peak_cpu)))
		rep_opt.close()


def variant_sites(aln):
    aln_opt = AlignIO.read(aln, 'fasta')
    aln_obj = AlignStats(aln_opt)
    return aln_obj.variant_sites_ratio()    

def build_raxml(aln, raxml_t):
    if raxml_t != 0:
        raxml_dir = os.getcwd()+'/raxml_tree'
        if os.path.exists(raxml_dir):
            shutil.rmtree(raxml_dir)
            os.makedirs(raxml_dir)

        subprocess.call('cp {} {}/.'.format(aln, raxml_dir), shell = True)      
        opt_name = aln.replace('.fna', '')
        cmd = 'raxmlHPC-PTHREADS-SSE3 -T {} -f a -# 100 -p 12345 -x 12345 -s {} -m GTRGAMMA -n {} -w {}'.format(raxml_t, aln, opt_name, raxml_dir)
        # sys.stdout.write(cmd)
        subprocess.call(cmd, shell = True)
    else:
        return 'No tree is built!'

def main():

    parser = argparse.ArgumentParser('contigs_based', 'reads_based', 'alignment_assessing')
    # Holding all info necessary to parse the command line into python data type
    subparsers = parser.add_subparsers(help = 'program mode', dest = 'mode')
    # Creating sub-commands

    add_blast_approach_cmd_options(subparsers) # add commands to blast-approach functionality
    add_mapping_approach_cmd_options(subparsers) # add commands to mapping-approach functionality
    # add_alignment_assessing_cmd_option(subparsers)

    args = parser.parse_args() # Added arguments are called and stored in args

    if args.mode == 'contigs_based':

        generate_par_report(args)
        generate_query_set_and_mf_file(args.contigs)
        make_blast_db(args.ref_fna)
        blast_(args.ref_fna, 'INTER_QuerySet.fna', args.blast_threads)
        QC_on_blastn('INTER_blast_opt_tmp.tab', args.alignment_length, args.identity)
        bowtie2_build(args.ref_fna)
        if args.ancient_sample != 'None':    
            create_mapping(args)
            parallel_filter(args)
            parallel_consensus(args)
        else:
            pass
        ref_ctig_dict = SeqIO.to_dict(SeqIO.parse(open(args.ref_fna), "fasta"))
        genomes_contigs_dict = concatenate_contigs("INTER_blastn_opt_tmp_cleaned.tab", 'INTER_genome_contigs_mp.csv', args.ref_fna)
        rec_genomes_lst = write_one_file(genomes_contigs_dict, args.ref_fna)
        ref_ctig_dict_ = {g: str(ref_ctig_dict[g].seq) for g in ref_ctig_dict}
        ref_seq_record = [SeqRecord(Seq(reorder_contigs(ref_ctig_dict, ref_ctig_dict_), generic_dna),\
        id = args.ref_fna, description = '')]
        if args.ancient_sample != 'None':
            agenomes = subprocess.getoutput("ls *consensus.fna").split("\n")
            agenomes_dict = {'a__'+g.replace('_consensus.fna','').replace('filtered_INTER_', ''): {str(c).replace('_consensus', ''): str(SeqIO.to_dict(SeqIO.parse(open(g),\
                'fasta'))[c].seq).replace('N', '-')\
            for c in SeqIO.to_dict(SeqIO.parse(open(g), 'fasta'))}\
            for g in agenomes}
            agenomes_lst = [SeqRecord(Seq(reorder_contigs(ref_ctig_dict, agenomes_dict[ag]), generic_dna),\
              id = ag, description = '')\
            for ag in agenomes_dict]
        else:
            pass 
        rec_genomes_lst.extend(ref_seq_record)
        if args.ancient_sample != 'None':
            rec_genomes_lst.extend(agenomes_lst)
        else:
            pass

        SeqIO.write(rec_genomes_lst, '{}-GenomeAln_ContigsBased.fna'.format(args.ref_fna.split('.')[0]), 'fasta') 
        if args.trim_reads_end != None and args.trimmed_reads != False and args.ancient_sample != 'None':
            filtered_bam = subprocess.getoutput('ls filtered_INTER*.bam.sorted').split('\n')
            for bam in filtered_bam:
                output_trimmed_reads(args.trim_reads_end, bam)
        else:
            pass
        dir_name = 'intermediates'
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        else:
            print("Folder {} already exists !".format(dir_name))
        subprocess.call('mv *INTER* ./{}'.format(dir_name), shell = True)
        subprocess.call('mv *bt2 ./{}'.format(dir_name), shell = True)
        subprocess.call('mv *nin ./{}'.format(dir_name), shell = True)
        subprocess.call('mv *nsq ./{}'.format(dir_name), shell = True)
        subprocess.call('mv *nhr ./{}'.format(dir_name), shell = True)
        subprocess.call('mv *bam ./{}'.format(dir_name), shell = True)
        out_stats('{}-GenomeAln_ContigsBased.fna'.format(args.ref_fna.split('.')[0]))

        raw_genome_aln = subprocess.getoutput('ls *GenomeAln*').split('\n')[0]
        

        # build_raxml(fine_tune_aln, args.raxml_threads)

    elif args.mode == 'reads_based':

        generate_par_report(args)
        bowtie2_build(args.ref_fna)
        if args.ancient_sample != 'None':       
            create_mapping_a(args)
            parallel_filter_a(args)
            parallel_consensus_a(args)
        else:
            pass
        create_mapping_m(args)
        parallel_filter_m(args)
        parallel_consensus_m(args)


        ref_ctig_dict = SeqIO.to_dict(SeqIO.parse(open(args.ref_fna), "fasta"))
        ref_ctig_dict_ = {g: str(ref_ctig_dict[g].seq) for g in ref_ctig_dict}      
        ref_seq_record = [SeqRecord(Seq(reorder_contigs(ref_ctig_dict, ref_ctig_dict_), generic_dna),\
         id = args.ref_fna, description = '')]

        all_rec_genomes = subprocess.getoutput("ls *consensus.fna").split("\n")
        all_rec_genomes_dict = {g.replace('filtered_INTER_', ''): {str(c).replace('_consensus', ''): str(SeqIO.to_dict(SeqIO.parse(open(g),\
        'fasta'))[c].seq).replace('N', '-')\
        for c in SeqIO.to_dict(SeqIO.parse(open(g), 'fasta'))}\
        for g in all_rec_genomes}

        all_rec_genomes_lst = [SeqRecord(Seq(reorder_contigs(ref_ctig_dict, all_rec_genomes_dict[ag]), generic_dna),\
        id = ag, description = '')\
        for ag in all_rec_genomes_dict]
        all_rec_genomes_lst.extend(ref_seq_record)
        SeqIO.write(all_rec_genomes_lst, '{}-GenomeAln_ReadsBased.fna'.format(args.ref_fna.split('.')[0]), 'fasta')     

        if args.trim_reads_end != None and args.trimmed_reads != False and args.ancient_sample != 'None':
            filtered_bam = subprocess.getoutput('ls filtered_INTER_a__*.bam.sorted').split('\n')
            for bam in filtered_bam:
                output_trimmed_reads(args.trim_reads_end, bam)
        else:
            pass        
        dir_name = 'intermediates'
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        else:
            print("Folder {} already exists !".format(dir_name))
        subprocess.call('mv *INTER* ./{}'.format(dir_name), shell = True)
        subprocess.call('mv *bt2 ./{}'.format(dir_name), shell = True)
        subprocess.call('mv *bam ./{}'.format(dir_name), shell = True)
        out_stats('{}-GenomeAln_ReadsBased.fna'.format(args.ref_fna.split('.')[0]))
        

        # build_raxml(fine_tune_aln, args.raxml_threads)
#------------------------------------------------------------------------------------------------------------
    # elif args.mode == 'alignment_assessing':
    #     if args.output_directory == None:
    #         if args.missinginfo_column_distribution:
    #             misv_all = AlignStats(AlignIO.read(args.genome_alignment, 'fasta')).column_MissingValue_dist()
    #             distribution(misv_all).savefig('Columndistribution_MissingValue_AllTaxa.png', dpi = 150)
    #             misv_a_msa = MultipleSeqAlignment([i for i in AlignIO.read(args.genome_alignment, 'fasta') if i.id.startswith('a__')])
    #             misv_m_msa = MultipleSeqAlignment([i for i in AlignIO.read(args.genome_alignment, 'fasta') if not i.id.startswith('a__')])
    #             misv_a = AlignStats(misv_a_msa).column_MissingValue_dist()
    #             misv_m = AlignStats(misv_m_msa).column_MissingValue_dist()
    #             distribution(misv_a).savefig('Columndistribution_MissingValue_AncientTaxa.png', dpi = 150)
    #             distribution(misv_m).savefig('Columndistribution_MissingValue_ModernTaxa.png', dpi = 150)
    #             aln_obj = AlignStats(AlignIO.read(args.genome_alignment, 'fasta'))
    #             Barplot(aln_obj.column_MissingValue_bar()).savefig('NumberColumn_cutoffs.png', dpi = 150)
    #         else:
    #             pass
    #         if args.variant_sites:
    #             aln_obj = AlignStats(AlignIO.read(args.genome_alignment, 'fasta'))
    #             var_sites = aln_obj.variant_sites_ratio()
    #             opt_var = open('Variant_sites_report.txt', 'w')
    #             opt_var.write("Biallelic sites: {}\n".format(str(var_sites[2])))
    #             opt_var.write("Biallelic sites ratio: {}\n".format(str(var_sites[0])))
    #             opt_var.write("Multiallelic sites: {}\n".format(str(var_sites[3])))
    #             opt_var.write("Multiallelic sites ratio: {}\n".format(var_sites[1]))
    #             opt_var.close()
    #         else:
    #             pass
    #         if args.variant_sites_window_sliding != 0:
    #             window_size = args.variant_sites_window_sliding
    #             aln_ipt = AlignIO.read(args.genome_alignment, 'fasta')
    #             aln_ipt_len = len(aln_ipt[1, :]) 
    #             opt_var = open('window_sliced_variant_sites_report.txt', 'w')
    #             init_aln = aln_ipt[:,:1]
    #             for i in range(0, aln_ipt_len, window_size):
    #                 sliced_aln = aln_ipt[:, i:i+window_size]
    #                 sliced_aln_obj = AlignStats(sliced_aln).variant_sites_ratio()
    #                 opt_var.write(str(sliced_aln_obj[2])+'\t'+ str(sliced_aln_obj[3]) + '\n')
    #                 if sliced_aln_obj[2] == 0 or sliced_aln_obj[3] == 0:
    #                     init_aln += sliced_aln
    #             opt_var.close()
    #             SeqIO.write(init_aln[:, 2:], 'hyper_free_aln_0.fna', 'fasta')
    #             ########
                # Add an argument to threshold 
                ########
        #     else:
        #         pass
        # else:
        #     print("The feature of re-directing opts to a specified directory is still in the development !")    
#----------------------------------------------------------------------------------------------------------------
        
    else:
        sys.exit('Oops, choose either contigs_based, reads_based, or alignment_tailor. Please check help menu')



if __name__ == '__main__':

    start = timeit.default_timer() # Starting timer
    main() # Excuting the main script
    stop = timeit.default_timer() # Ending timer 
    time_callapsed = stop - start # 
    par_file = open("Parameters_setting.txt", "a+") # Reopen parameter setting file for appending
    par_file.write("Time used:{}(s)\n".format(str(time_callapsed))) # Appending time collapsed to parameter setting file
    par_file.close()








