#!/usr/bin/env python


import json
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

def main():
    # parse parameters
    parser = argparse.ArgumentParser(description = 'Reconstruct whole-genome-level MSA from large-scale datasets.')
    parser.add_argument('config_file', help='Input the configuration file.')
    parser.add_argument('-r',
                        '--reference',
                        help='Input reference genome sequence in the fasta format.',
                        type = str,
                        default = None)
    parser.add_argument('-a',
                        '--ancient_metagenomes',
                        help = 'An ancient-metagenome folder containing sub-folders, each sub-folder \
                        contains sequencing reads from one ancient sample. [Support suffix of .fastq/.fastq.gz/.fastq.bz2]',
                        type = str,
                        default = None)
    parser.add_argument('-m',
                        '--modern_metagenomes',
                        help = 'A modern-metagenome folder containing sub-folders, each sub-folder \
                        contains sequencing reads from one modern sample. [Support suffix of .fastq/.fastq.gz/.fastq.bz2]',
                        type = str,
                        default = None)
    parser.add_argument('-g',
                        '--genome_assemlies',
                        help = 'A genome-assemblies folder containing assembled genomes, each assembled \
                        genome is stored in a fasta file. [Support FASTA format]',
                        type = str,
                        default = None)
    parser.add_argument('-o',
                        '--output_dir',
                        help = 'A folder for storing outputs.',
                        default = 'Mac_output')

    args = parser.parse_args()
    

    try:
        with open(args.config_file, 'r') as config_file:
            configs_list = json.loads(config_file.read())
    except Exception as e:
        print(e) 

    Inters = create_folder('Intermediates') # Create a folder for storing intermediate files
    
    configs_list[0]['db_dir'] = Inters
    configs_list[1]['db_dir'] = Inters
    configs_list[2]['db_dir'] = Inters

    if args.reference:
        ref_genome = abspath_finder(args.reference) # get the abs path of refseq
        configs_list[0]["ref_genome"] = ref_genome
        configs_list[1]["ref_genome"] = ref_genome
        configs_list[2]["ref_genome"] = ref_genome
    else:
        pass
    if args.ancient_metagenomes:
        a_samples = obtain_samples(args.ancient_metagenomes)
        configs_list[0]['param_set']['sample_list'] = a_samples
    elif len(configs_list[0]['param_set']['sample_list']) != 0:
        a_samples = obtain_samples(configs_list[0]['param_set']['sample_list'])
        configs_list[0]['param_set']['sample_list'] = a_samples
    else:   
        a_samples = None
    # get all ancient metagenome samples in abs path
    if args.modern_metagenomes:
        m_samples = obtain_samples(args.modern_metagenomes)
        configs_list[1]['param_set']['sample_list'] = m_samples
    elif len(configs_list[1]['param_set']['sample_list']) != 0:
        m_samples = obtain_samples(configs_list[1]['param_set']['sample_list'])
        configs_list[1]['param_set']['sample_list'] = m_samples   
    else:
        m_samples = None
    # get all modern metagenome samples in abs path
    if args.genome_assemlies:
        genomes = abspath_finder(args.genome_assemlies)
        configs_list[2]['param_set']['sample_list'] = genomes
    else:
        genomes = None
    # get the folder of genome assemblies
    
    opt_dir = create_folder(args.output_dir)


    inter_results = []
    # deal with working directory
    for configs in configs_list:
        dest = workflow(configs)
        inter_results.extend(dest)


    # merge fiiles in inter_results
    output_file = opt_dir + '/Mac_genome_MSA.fna'
    merge_all(inter_results, ref_genome, output_file)

def obtain_samples(folder_path):

    """
    It takes abs path of an input folder,
    and return a list of abs paths of sub-folders
    """
    sample_list = []

    for i in subprocess.getoutput('ls {}/'.format(abspath_finder(folder_path))).split('\n'):
        single_sample_path = abspath_finder(folder_path + '/' + i)
        sample_list.append(single_sample_path)

    return sample_list 


def workflow(configs):
    mode = configs['mode']
    # build database
    db_dest = build_db_file(mode, configs['ref_genome'], configs['db_dir'], configs['param_set']['sample_list'])
    print(db_dest)
    # build mapping with database files
    reconstructed_genome = build_mapping(db_dest, **configs)
    
    return reconstructed_genome



def build_db_file(mode, ref_genome, db_dir, samples):
    """
    Args:
      ref_genome: input file path
      db_dir: destination directory path

    Return:
      destination files path
    """
    params = {
        'ref': ref_genome,
        'dest': db_dir + '/'+ ref_genome.split('/')[-1], 'Samples': samples,
    }

    if mode == 'reads' and params['Samples']:
        # build command
        cmd = "bowtie2-build -q {ref} {dest}".format(**params)
        run_cmd_in_shell(cmd)
    elif (mode == 'contigs') and params['Samples']:
        cmd = "makeblastdb -in {ref} -dbtype nucl -out {dest}".format(**params)
        run_cmd_in_shell(cmd)
    else:
    	print('database of {} mode was skipped!\n'.format(mode))
    	pass

        
    # run the command
    # TODO: exception handling


    return params['dest']


def build_mapping(db_dest, **kwargs):
    mode = kwargs['mode']
    param_set = kwargs['param_set']
    opt_all_files = []
    if mode == 'reads':
        age_type = kwargs['age_type']
        if (age_type == 1) and param_set['sample_list']:
            # Using ancient-specific param_set
            bam_files = bwt2_batch_mapping(param_set['sample_list'], db_dest, param_set['thread'], param_set['m_mode'], param_set['nproc']) # parameters specific to mapping ancient samples
            print('bam_files result: \n', bam_files)
            filtered_bams = batch_bam_filter(bam_files, param_set['min_q'], param_set['min_l'], param_set['max_snp_edist'], param_set['nproc'])
            print('filtered bam files result: \n', filtered_bams)
            opt_files = batch_consensus_builder(filtered_bams, param_set['min_c'], param_set['t_dist'], param_set['domi_ale_frq'], param_set['nproc'])
            print('reconstructed fasta files result: \n', opt_files)
            opt_all_files.extend(opt_files)


        elif (age_type == 2) and param_set['sample_list']:
            # Using modern-specific param_set
            bam_files = bwt2_batch_mapping(param_set['sample_list'], db_dest, param_set['thread'], param_set['m_mode'], param_set['nproc']) # parameters specific to mapping ancient samples
            print('bam_files result: \n', bam_files)
            filtered_bams = batch_bam_filter(bam_files, param_set['min_q'], param_set['min_l'], param_set['max_snp_edist'], param_set['nproc'])
            print('filtered bam files result: \n', filtered_bams)
            opt_files = batch_consensus_builder(filtered_bams, param_set['min_c'], None, param_set['domi_ale_frq'], param_set['nproc'])
            print('reconstructed fasta files result: \n', opt_files)
            opt_all_files.extend(opt_files)
 
    elif (mode == 'contigs') and param_set['sample_list']:
        configs = kwargs
        print('Generate mp file and concatenate ..........\n')
        ctigs_concatenated = generate_query_set_and_mp_file(configs['db_dir'], param_set['sample_list'])
        mp_file = ctigs_concatenated[0]
        query_set = ctigs_concatenated[1]
        db_dest = configs['db_dir'] + '/'+ configs['ref_genome'].split('/')[-1]
        print('Run blastn ........\n')
        raw_blastn_tab = blast_genomes(configs['db_dir'], db_dest, query_set, param_set['b_threads'])
        print('Run blastn qc ............\n')
        cleaned_blastn_tab = QC_on_blastn(raw_blastn_tab, param_set['b_len'], param_set['b_iden'], configs['db_dir'])
        print('Merge ..............\n')
        genomes_contigs_dict = concatenate_contigs(cleaned_blastn_tab, mp_file, configs['ref_genome'])
        print('Output single files ........\n')
        opt_files = output_single_files(genomes_contigs_dict, configs['ref_genome'], configs['db_dir'])
        opt_all_files.extend(opt_files)
    return opt_all_files


def single_mapping(output_dir, db_dest, sample, threads, m_mode):
    """
    Args: 
        ref_genome: the input file name of reference genome
        sample: the input folder containing metagenomic reads of a sample
        threads: tunable parameter for #threads in bowtie2
        m_mode: tunable parameter for searching mode in bowtie2

    """
    print("Mapping single sample {}".format(sample))
    m_mode = m_mode.replace("*", "") # Removing escape character for m_mode 
    suffix = detect_reads_suffix(sample)
    opt_raw_bam = output_dir + '/' + sample.split('/')[-1] + '.bam'

    if suffix == "bz2":
        cmd = 'bzcat {}/*fastq.bz2 | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | samtools view -bS - > {}'.format(sample, db_dest, threads, ' '.join(m_mode.split(',')), opt_raw_bam)
    elif suffix == "gz":
        cmd = 'zcat {}/*fastq.gz | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | samtools view -bS - > {}'.format(sample, db_dest, threads, ' '.join(m_mode.split(',')), opt_raw_bam)
    elif suffix == "fastq":
        cmd = 'cat {}/*fastq | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | samtools view -bS - > {}'.format(sample, db_dest, threads, ' '.join(m_mode.split(',')), opt_raw_bam)
    else:
        sys.exit("Reads have to be in the form of .bz2, .gz or .fastq!")

    run_cmd_in_shell(cmd)
    return opt_raw_bam

def bwt2_batch_mapping(sample_list, db_dest, threads, m_mode, processors):
    
  
    print("Batch mapping samples:\n {}".format("\n".join(sample_list)))
    
    proc_num = len(sample_list)
    output_dir = '/'.join(db_dest.split('/')[:-1])
    params = [[output_dir, db_dest, sample_list[i], threads, m_mode] for i in range(proc_num)]

    return multi_map(processors, single_mapping, params)


def run_cmd_in_shell(cmd):
    print('running command: ', cmd)
    subprocess.call(cmd, shell=True)


def detect_reads_suffix(reads_foler):
    """
    Args: reads folder
    Return: suffix of reads file inside the folder
    """
    reads_file = subprocess.getoutput('ls {}/*fastq*'.format(reads_foler)).split('\n')[0]
    return reads_file.split('.')[-1]

def single_bam_filter(output_dir, raw_bam, min_q, min_l, max_snp_edist):
    """
    Args: a given bam coupled with parameters for QC
        raw_bam: the raw bam file output from bwt2_batch_mapping()
        min_q: a tunable parameter for minimum quality
        min_l: a tunable parameter for minimum length
        max_snp_edist: a tunable parameter for maximum SNP edit distance 

    Return: a filtered bam [filtered_sample.bam] 
    """
    opt_filtered_bam = output_dir + '/' + 'filtered____'+ raw_bam.split('/')[-1]
    filter_path = os.path.dirname(__file__)+"/cmseq/filter.py"
    cmd = "samtools view -h -F 0x4 {} | python3 {} --minqual {} --minlen {} --maxsnps {} > {}".format(raw_bam, filter_path,\
     min_q, min_l, max_snp_edist, opt_filtered_bam)
    run_cmd_in_shell(cmd)
    return opt_filtered_bam

def batch_bam_filter(bams, min_q, min_l, max_snp_edist, processors):
    """
    Parallelize single_bam_filter()
    """
    print("Batch filter bams:\n {}".format("\n".join(bams)))

    proc_num = len(bams)
    output_dir = '/'.join(bams[0].split('/')[:-1])
    params = [[output_dir, bams[i], min_q, min_l, max_snp_edist] for i in range(proc_num)]

    
    return multi_map(processors, single_bam_filter, params)

def single_consensus_builder(output_dir, filtered_bam, min_c, t_dist, domi_ale_frq):
    """
    Args: a filtered bam coupled with parameters for QC
          filtered_bam: A filtered bam file output from batch_bam_filter()
          min_c: a tunable parameter for minimum coverage
          t_dist: a tunable parameter for trimming distance
          domi_ale_frq: a tunable parameter for dominant allele frequency 


    return: consensus sequences in a fasta file
    """
    opt_filtered_consensus = output_dir + '/' + 'consensus_'+ filtered_bam.split('/')[-1].replace('.bam', '.fna')
    consensus_path = os.path.dirname(__file__)+"/cmseq/consensus.py"

    if t_dist != None:

        cmd = "python3 {} {} --sortindex --mincov {} --trim {} --dominant_frq_thrsh {} > {}".format(consensus_path, filtered_bam, min_c, t_dist, domi_ale_frq, opt_filtered_consensus)
        run_cmd_in_shell(cmd)
    else:
        cmd = "python3 {} {} --sortindex --mincov {} --dominant_frq_thrsh {} > {}".format(consensus_path, filtered_bam, min_c, domi_ale_frq, opt_filtered_consensus)
        run_cmd_in_shell(cmd)

    return  opt_filtered_consensus 

def batch_consensus_builder(filtered_bams, min_c, t_dist, domi_ale_frq, processors):

    """
    Parallelize single_consensus_builder()
    """
    print("batch building consensus:\n {}".format("\n".join(filtered_bams)))
    proc_num = len(filtered_bams)
    output_dir = '/'.join(filtered_bams[0].split('/')[:-1])
    params = [[output_dir, filtered_bams[i], min_c, t_dist, domi_ale_frq] for i in range(proc_num)]


    return multi_map(processors, single_consensus_builder, params)


# # def output_trimmed_reads(trim_pos, bam_file):
# #     """
# #     This feature is optional.
# #     If chosen, trimmed reads extracted from filtered bams and re-direct to fastq files which
# #     are stored in outputs/trimmed_reads 
# #     """
# #     in_samfile = pysam.AlignmentFile(bam_file, 'rb')
# #     bam_file_opt = 'trimmed_'+'_'.join(bam_file.split('.')[0].split('_')[2:])+'.fastq'
# #     out_fastq = open(bam_file_opt, 'w')
# #     left, right = trim_pos.split(':')
# #     for read in in_samfile.fetch():
# #         if read.is_reverse:
# #             read_name = read.query_name
# #             read_seq = str(Seq(read.query_alignment_sequence[int(left): -int(right)], generic_dna)\
# #                 .reverse_complement())
# #             read_qual = read.qual[int(left): -int(right)][::-1]

# #             out_fastq.write('@'+read_name+'\n')
# #             out_fastq.write(read_seq+'\n')
# #             out_fastq.write("+\n")
# #             out_fastq.write(read_qual+'\n')
# #         else:
# #             read_name = read.query_name
# #             read_seq = str(Seq(read.query_alignment_sequence[int(left): -int(right)], generic_dna))
# #             read_qual = read.qual[int(left): -int(right)]

# #             out_fastq.write('@'+read_name+'\n')
# #             out_fastq.write(read_seq+'\n')
# #             out_fastq.write("+\n")
# #             out_fastq.write(read_qual+'\n')
# #     out_fastq.close()


# ###############################################################################################################
# # Above are all about bowtie2 process, and below are all about blastn process
# ###############################################################################################################

def generate_query_set_and_mp_file(output_dir, contigs_folder):
    
    """
    Input: 'contigs' -> a folder of fasta files, each represents a genome
    Program: 
    1) It merges all individual fasta files into one fasta file, called 'INTER_QuerySet.fna',as an intermediate.
    2) It creates a mapping file delimited by 'comma', 1st column is genome name (i.e. file name) and 2nd column is 
    name of the contig which belongs to the genome.
    Output: 1) A merged fasta file, INTER_QuerySet.fna, for downstream analysis as input
            2) A mapping file, 1st column is genome name and 2nd column is contig name 
    """

    mp_file = output_dir + '/genome_contigs_mp.csv'
    opened_mp = open(mp_file, 'w')

    subprocess.call('cat {}/* > {}/QuerySet.fna'.format(contigs_folder, output_dir), shell = True)
    all_fna = subprocess.getoutput('ls {}/*'.format(contigs_folder)).split('\n')

    for f in all_fna:
        opened_fna = SeqIO.to_dict(SeqIO.parse(open(f), 'fasta'))
        for h in opened_fna:
            row = f.split('/')[-1] + ',' + h + '\n'
            opened_mp.write(row)
    opened_mp.close()

    return mp_file, '{}/QuerySet.fna'.format(output_dir)

def blast_genomes(output_dir, db_dest, query_set, threads):

    outfmt = '6 qaccver saccver pident length mismatch gapopen qstart qend\
                    sstart send evalue bitscore qseq sseq'
    cmd = "blastn -db {} -query {} -outfmt '{}' -num_threads {}\
            -word_size 9 -out {}/raw_blastn_opt.tab".format(db_dest, query_set, outfmt, threads, output_dir)

    run_cmd_in_shell(cmd)
    return '{}/raw_blastn_opt.tab'.format(output_dir)

def QC_on_blastn(blast_tab, length, pid, output_dir):

    """
    It applies QC on blast output
    Input: 1)'INTER_blast_opt_tmp.tab', 2) length of hits, 3) identity percentage of hits
    Program: call script 'bo6_screen.py'
    Output: filtered blast output, 'INTER_blastn_opt_tmp_cleaned.tab' 
    """

    cmd = 'cut -f 1-14 {} | bo6_screen.py --length {} --pid {} > {}/cleaned_blastn_opt.tab'.format(blast_tab, length, pid, output_dir)
    run_cmd_in_shell(cmd)
    return '{}/cleaned_blastn_opt.tab'.format(output_dir)   
   


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

def reorder_contigs_2(ref_fna_dict, recon_genome_dict):
    """
    This version is to output single fasta files for each homology-guided reconstructed genomes
    Output a list of SeqRecord which contain all contigs.
    """
    ref_ctig_coordinate=sorted([(i, len(ref_fna_dict[i].seq)) for i in ref_fna_dict], key=lambda x: x[0])

    contigs_list = []
    for c in ref_ctig_coordinate:
        if c[0] in recon_genome_dict:
            ctig_record = SeqRecord(Seq(recon_genome_dict[c[0]], generic_dna), id = c[0] + '_consensus', description = '')
            contigs_list.append(ctig_record)
        else:
            ctig_seq = len(ref_fna_dict[c[0]]) * '-'
            ctig_record = SeqRecord(Seq(ctig_seq, generic_dna), id = c[0] + '_consensus', description = '')
            contigs_list.append(ctig_record)

    return contigs_list

def output_single_files(genomes_dict_, ref_fna, output_dir):
    ref_ctig_dict = SeqIO.to_dict(SeqIO.parse(open(ref_fna), "fasta"))
    files = []
    for g in genomes_dict_:
        recon_seq = reorder_contigs_2(ref_ctig_dict, genomes_dict_[g])
        opt_file_name = '{}/homolog____{}'.format(output_dir, g)
        files.append(opt_file_name)
        SeqIO.write(recon_seq, opt_file_name ,'fasta')

    return files


def merge_all(files_lst, ref_fna, output_file):
    sorted_headers = sorted([header for header in SeqIO.to_dict(SeqIO.parse(open(ref_fna), "fasta"))])
    seq_records = []
    refseq_to_dict = SeqIO.to_dict(SeqIO.parse(open(ref_fna), "fasta"))
    ref_seq_concatenated = ''.join(str(refseq_to_dict[header].seq) for header in sorted_headers)
    seq_records.append(SeqRecord(Seq(ref_seq_concatenated, generic_dna), id = ref_fna.split('/')[-1], description= ''))
    for file in files_lst:
        file_to_dict = SeqIO.to_dict(SeqIO.parse(open(file), "fasta"))
        seq_name = file.split('/')[-1].split('____')[-1]
        concatenated_seq = ''
        for header in sorted_headers:
            header = header + '_consensus'
            concatenated_seq += str(file_to_dict[header].seq)
        single_seq_record = SeqRecord(Seq(concatenated_seq, generic_dna), id = seq_name, description = '')
        seq_records.append(single_seq_record)
    SeqIO.write(seq_records, output_file, 'fasta')

    return output_file

def abspath_finder(file):
    # Convert wahtever path into abs path
    return os.path.abspath(file)

def create_folder(name):
    folder_name = abspath_finder(name) # determine the abs path
    if os.path.exists(folder_name):
        shutil.rmtree(folder_name)
    os.makedirs(folder_name)
    return folder_name

    
if __name__ == "__main__":
    main()



