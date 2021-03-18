#!/usr/bin/env python


import argparse
import itertools
import json
import logging
import math
import os
import pysam
import shutil
import subprocess
import sys
import timeit

from collections import defaultdict
from functools import partial
from logging.config import fileConfig

from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from .utils.data_types import AncientReadsType, ContigsType, ModernReadsType
from .utils import utils, AlignStats, SNP_rates

log_config = "/".join(os.path.abspath(__file__).split('/')[:-1]) + '/metaclock_configs/logging_config.ini'
fileConfig(log_config)
logger = logging.getLogger()


def main():
    # If the if_clean is True, we skip creation of DB files and bam files if
    # they exist.

    parser = argparse.ArgumentParser(description = 'Reconstruct whole-genome-level MSA from large-scale datasets.')
    parser.add_argument('config_file', help='Input the configuration file.')
    parser.add_argument('-r',
                        '--reference',
                        help='Input reference genome sequence in the fasta format.',
                        type = str,
                        default = None)
    parser.add_argument('-a_ipt',
                        '--ancient_metagenomes',
                        help = 'Input an ancient-metagenome folder containing sub-folders, each sub-folder \
                        contains sequencing reads from one ancient sample. [Support suffix of .fastq/.fastq.gz/.fastq.bz2]',
                        type = str,
                        default = None)
    parser.add_argument('-m_ipt',
                        '--modern_metagenomes',
                        help = 'Input a modern-metagenome folder containing sub-folders, each sub-folder \
                        contains sequencing reads from one modern sample. [Support suffix of .fastq/.fastq.gz/.fastq.bz2]',
                        type = str,
                        default = None)
    parser.add_argument('-g_ipt',
                        '--genome_assemlies',
                        help = 'Input a genome-assemblies folder containing assembled genomes, each assembled \
                        genome is stored in a fasta file. [Support FASTA format]',
                        type = str,
                        default = None)
    parser.add_argument('-o',
                        '--output_dir',
                        help = 'Specify an output folder for storing results. Default: [Mac_output] in the working directory.',
                        default = 'Mac_output')

    parser.add_argument('-int',
                        '--intermediate_dir',
                        help = 'Specify an intermediate folder for storing intermediate files. Default: [intermediates] in the working directory.',
                        default = None)
    parser.add_argument('-c',
                        '--clean',
                        help='Clean intermediate files and rerun from the beginning, \
                        otherwise rerun from the intermediate files',
                        action='store_true')
    parser.add_argument('-a',
                        '--authentication',
                        help='Autheticate the anicent origin of genomic information used in genome reconstruction.',
                        action='store_true')
    parser.add_argument('-snv_rate',
                        '--SNV_rate',
                        help='Estimate pairwise SNV rates between samples.',
                        action='store_true')

    args = parser.parse_args()



    if_clean = args.clean
    logger.info("Clean intermediate files: {}".format(if_clean))
    if_est_snv = args.SNV_rate
    logger.info("Estimate pairwise SNV rates: {}".format(if_est_snv))
    if_authenticate = args.authentication
    logger.info("Authenticate ancient origin of reads: {}".format(if_authenticate))

    #Loading config file
    try:
        with open(args.config_file, 'r') as config_file:
            configs_list = json.loads(config_file.read())
            logger.info("Loading configuration file: completed!")
    except Exception as e:
        logger.error(e)

    dummy_data = {
        'intermediate': '',
        'reference_genome': '',
        'samples': '',
    }
    if not configs_list.get('ancient_reads', None):
        configs_list['ancient_reads'] = dummy_data
    if not configs_list.get('modern_reads', None):
        configs_list['modern_reads'] = dummy_data
    if not configs_list.get('contigs', None):
        configs_list['contigs'] = dummy_data

    if args.intermediate_dir \
        or (not configs_list['ancient_reads']['intermediate'] \
            or not configs_list['modern_reads']['intermediate'] \
            or not configs_list['contigs']['intermediate']):
        if args.intermediate_dir:
            # Check if intermediate path is given by cmd, use this path if yes
            inters = create_folder(args.intermediate_dir)
        else:
            inters = create_folder('intermediates')

        # Add absolute path of intermediate folder to sample processing job
        configs_list['ancient_reads']['intermediate'] = inters
        configs_list['modern_reads']['intermediate'] = inters
        configs_list['contigs']['intermediate'] = inters
        logger.info('Intermediate files folder:\n{}'.format(inters))

    if args.reference:
        # If reference is given by command, overwrite config file
        # with updated data from cmd
        ref_genome = abspath_finder(args.reference) # get the abs path of refseq
        configs_list['ancient_reads']["reference_genome"] = ref_genome
        configs_list['modern_reads']["reference_genome"] = ref_genome
        configs_list['contigs']["reference_genome"] = ref_genome
        logger.info('Reading reference genome from:\n{}'.format(ref_genome))
    elif not configs_list['ancient_reads']['reference_genome'] \
        and not configs_list['modern_reads']['reference_genome'] \
        and not configs_list['contigs']['reference_genome']:
        sys.exit("Please input a reference genome in fasta format!")

    if args.ancient_metagenomes:
        #If path of ancient samples is given by cmd, it obtains
        # all samples' abs paths from parsed args and store them in a list,
        # and overwrites config with a list of sample paths.
        configs_list['ancient_reads']['samples'] = obtain_samples(
            args.ancient_metagenomes)
    elif len(configs_list['ancient_reads']['samples']) != 0:
        # If the path is given by config file, it obtains
        # all samples' abs paths from path written in config file
        # and store them in a list, and overwrites config
        # with a list of sample paths
        configs_list['ancient_reads']['samples'] = obtain_samples(
            configs_list['ancient_reads']['samples'])
    else:
        configs_list['ancient_reads']['samples'] = None

    if args.modern_metagenomes:
        # If path of modern samples is given by cmd, it obtains
        # all samples' abs paths from parsed args and store them in a list,
        # and overwrites config with a list of sample paths.
        configs_list['modern_reads']['samples'] = obtain_samples(args.modern_metagenomes)
    elif len(configs_list['modern_reads']['samples']) != 0:
        # If the path is given by config file, it obtains
        # all samples' abs paths from path written in
        # config file and store them in a list,
        # and overwrites config with a list of sample paths
        configs_list['modern_reads']['samples'] = obtain_samples(
            configs_list['modern_reads']['samples'])
    else:
        configs_list['modern_reads']['samples'] = None

    if args.genome_assemlies:
        configs_list['contigs']['samples'] = abspath_finder(args.genome_assemlies)
    elif len(configs_list['contigs']['samples']) != 0:
        configs_list['contigs']['samples'] = abspath_finder(
            configs_list['contigs']['samples'])
    else:
        configs_list['contigs']['samples'] = None

    opt_dir = create_folder(args.output_dir)
    logger.info('Outputs folder: {}'.format(opt_dir))

    inter_results = []
    str_to_class = {
        'ancient_reads': AncientReadsType,
        'contigs': ContigsType,
        'modern_reads': ModernReadsType,
    }
    for k,v in configs_list.items():
        if not v.get('parameter_set', None):
            continue
        configs = str_to_class[k](**v)
        dest = workflow(configs, if_clean, if_authenticate, opt_dir)
        inter_results.extend(dest)


    output_file = opt_dir + '/' + opt_dir.split('/')[-1] +'.fna'
    Mac_final = merge_all(inter_results, ref_genome, output_file)

    opt_stats = opt_dir + '/mac_stats.tsv'
    utils.out_stats(Mac_final, opt_stats = opt_stats)

    if if_est_snv:
        nproc = configs_list['ancient_reads']['parameter_set']['nproc']
        SNP_rates.plot_snv_rates_main(Mac_final, nproc, output_dir = opt_dir)
    logger.info('Genome alignment reconstuction is completed and welcome back !')


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


def workflow(configs, if_clean, if_authenticate, opt_dir):
    mode = configs.input_type
    # build database
    db_dest = []
    if not if_clean:
        if mode == 'reads':
            db_filenames = get_bowtie2_db_files(configs.intermediate, configs.reference_genome)
        elif mode == 'contigs':
            db_filenames = get_blastn_db_files(configs.intermediate, configs.reference_genome)

        skip = True
        for db_filename in db_filenames['db_files']:
            if not os.path.exists(db_filename):
                skip = False
                break
        if skip:
            # DB files exist, we skip db creation
            db_dest = db_filenames['db_prefix']
            logger.info('Using existing database files: {}'.format(db_dest))
    if not db_dest:
        logger.info('Creating new database files')
        db_dest = build_db_file(mode, configs.reference_genome, configs.intermediate, configs.samples)
    logger.info('Database files path:\n{}'.format(db_dest))
    # build mapping with database files
    reconstructed_genome = build_mapping(db_dest, if_clean, if_authenticate, opt_dir, configs)

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
        logger.info("Creating database for bowtie2: completed!")
    elif (mode == 'contigs') and params['Samples']:
        cmd = "makeblastdb -in {ref} -dbtype nucl -out {dest}".format(**params)
        run_cmd_in_shell(cmd)
        logger.info("Creating database for blastn: completed!")
    else:
        logger.info('Database of {} mode was skipped!'.format(mode))

    return params['dest']


def build_mapping(db_dest, if_clean, if_authenticate, opt_dir, configs):
    mode = configs.input_type
    param_set = configs.param_set
    opt_all_files = []
    if mode == 'reads':
        age_type = configs.age_type
        if (age_type == 1) and configs.samples:
            # Using ancient-specific param_set
            bam_files = bwt2_batch_mapping(
                configs.samples, db_dest, param_set['bowtie2_threads'],
                param_set['search_report_mode'], param_set['nproc'], if_clean)

            logger.info('Raw bam files: \n{}'.format('\n'.join(bam_files)))
            filtered_bams = batch_bam_filter(
                bam_files, param_set['minimum_mapping_quality'],
                param_set['minimum_mapping_length'],
                param_set['maximum_snp_edit_distance'], param_set['nproc'])
            logger.info(
                'Filtered bam files: \n{}'.format('\n'.join(filtered_bams)))
            opt_files = batch_consensus_builder(
                filtered_bams, param_set['minimum_coverage'],
                param_set['trim_distance'],
                param_set['dominant_allele_frequency'], param_set['nproc'])
            logger.info('Reconstructed fasta files: \n{}'.format('\n'.join(opt_files)))
            opt_all_files.extend(opt_files)
            if param_set['output_trimmed_reads'] == 1:
                for bam in filtered_bams:
                    sorted_bam = bam + '.sorted'
                    logger.info("Outputing trimmed reads from {}".format(sorted_bam))
                    output_trimmed_reads(param_set['trim_distance'], sorted_bam)
            else:
                pass

            if if_authenticate:
                damage_patterns = authenticate(
                    configs.intermediate, configs.reference_genome,
                    filtered_bams)
                G2A_files = damage_patterns[0]
                C2T_files = damage_patterns[1]
                utils.draw_damage_pattern(G2A_files, C2T_files, opt_dir)
            else:
                pass
        elif (age_type == 2) and configs.samples:
            # Using modern-specific param_set
            bam_files = bwt2_batch_mapping(
                configs.samples, db_dest, param_set['bowtie2_threads'],
                param_set['search_report_mode'], param_set['nproc'], if_clean)

            logger.info('Raw bam files: \n{}'.format('\n'.join(bam_files)))
            filtered_bams = batch_bam_filter(
                bam_files, param_set['minimum_mapping_quality'],
                param_set['minimum_mapping_length'],
                param_set['maximum_snp_edit_distance'], param_set['nproc'])
            logger.info('Filtered bam files: \n{}'.format('\n'.join(filtered_bams)))
            opt_files = batch_consensus_builder(
                filtered_bams, param_set['minimum_coverage'], None,
                param_set['dominant_allele_frequency'], param_set['nproc'])
            logger.info('Reconstructed fasta files: \n{}'.format('\n'.join(opt_files)))
            opt_all_files.extend(opt_files)
    elif (mode == 'contigs') and configs.samples:
        ctigs_concatenated = generate_query_set_and_mp_file(
            configs.intermediate, configs.samples)
        logger.info('Generating mp file and concatenating genomes: completed!\n')

        mp_file = ctigs_concatenated[0]
        query_set = ctigs_concatenated[1]
        db_dest = configs.intermediate + '/'+ configs.reference_genome.split('/')[-1]
        raw_blastn_tab = blast_genomes(
            configs.intermediate, db_dest, query_set,
            param_set['blastn_threads'], if_clean)
        logger.info('Genearating raw blastn results: completed!\n')

        cleaned_blastn_tab = QC_on_blastn(
            raw_blastn_tab, param_set['homolog_length'],
            param_set['homolog_identity'], configs.intermediate)
        logger.info('Applying QC on raw blastn results: completed!\n')

        genomes_contigs_dict = concatenate_contigs(
            cleaned_blastn_tab, mp_file, configs.reference_genome)
        logger.info('Reordering homologous fragments: completed!\n')

        opt_files = output_single_files(
            genomes_contigs_dict, configs.reference_genome, configs.intermediate)
        logger.info('Outputing single homologs files: completed!\n')

        opt_all_files.extend(opt_files)
    return opt_all_files


def single_mapping(output_dir, db_dest, sample, threads, m_mode, if_clean):
    """
    Args:
        ref_genome: the input file name of reference genome
        sample: the input folder containing metagenomic reads of a sample
        threads: tunable parameter for #threads in bowtie2
        m_mode: tunable parameter for searching mode in bowtie2

    """

    m_mode = m_mode.replace("*", "") # Removing escape character for m_mode
    suffix = detect_reads_suffix(sample)
    opt_raw_bam = get_single_bam_file(sample.split('/')[-1], output_dir)

    if os.path.exists(opt_raw_bam) and not if_clean:
        return opt_raw_bam

    if suffix == "bz2":
        cmd = 'bzcat {}/*fastq.bz2 | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | samtools view -bS - > {}'.format(sample, db_dest, threads, ' '.join(m_mode.split(',')), opt_raw_bam)
    elif suffix == "gz":
        cmd = 'zcat {}/*fastq.gz | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | samtools view -bS - > {}'.format(sample, db_dest, threads, ' '.join(m_mode.split(',')), opt_raw_bam)
    elif suffix == "fastq":
        cmd = 'cat {}/*fastq | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | samtools view -bS - > {}'.format(sample, db_dest, threads, ' '.join(m_mode.split(',')), opt_raw_bam)
    else:
        sys.exit("Reads have to be in the form of .bz2, .gz or .fastq!")

    run_cmd_in_shell(cmd)
    logger.info("Mapping single sample {}: completed!".format(sample))
    return opt_raw_bam

def bwt2_batch_mapping(sample_list, db_dest, threads, m_mode, processors, if_clean):

    # logger.info("Batch mapping samples:\n {}".format("\n".join(sample_list)))

    proc_num = len(sample_list)
    output_dir = '/'.join(db_dest.split('/')[:-1])
    params = [[output_dir, db_dest, sample_list[i], threads, m_mode, if_clean] for i in range(proc_num)]

    return utils.multi_map(processors, single_mapping, params)


def run_cmd_in_shell(cmd):
    logger.info('Running command: {}'.format(cmd))
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
    filter_path = "filter.py"
    cmd = "samtools view -h -F 0x4 {} | {} --minqual {} --minlen {} --maxsnps {} > {}".format(raw_bam, filter_path,\
     min_q, min_l, max_snp_edist, opt_filtered_bam)
    run_cmd_in_shell(cmd)
    return opt_filtered_bam

def batch_bam_filter(bams, min_q, min_l, max_snp_edist, processors):
    """
    Parallelize single_bam_filter()
    """
    # logger.info("Batch filter bams:\n {}".format("\n".join(bams)))

    proc_num = len(bams)
    output_dir = '/'.join(bams[0].split('/')[:-1])
    params = [[output_dir, bams[i], min_q, min_l, max_snp_edist] for i in range(proc_num)]


    return utils.multi_map(processors, single_bam_filter, params)

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
    consensus_path = "consensus.py"

    if t_dist != None:

        cmd = "{} {} --sortindex --mincov {} --trim {} --dominant_frq_thrsh {} > {}".format(consensus_path, filtered_bam, min_c, t_dist, domi_ale_frq, opt_filtered_consensus)
        run_cmd_in_shell(cmd)
    else:
        cmd = "{} {} --sortindex --mincov {} --dominant_frq_thrsh {} > {}".format(consensus_path, filtered_bam, min_c, domi_ale_frq, opt_filtered_consensus)
        run_cmd_in_shell(cmd)

    return  opt_filtered_consensus

def batch_consensus_builder(filtered_bams, min_c, t_dist, domi_ale_frq, processors):
    """
    Parallelize single_consensus_builder()
    """
    # logger.info("Batch building consensus:\n {}".format("\n".join(filtered_bams)))
    proc_num = len(filtered_bams)
    output_dir = '/'.join(filtered_bams[0].split('/')[:-1])
    params = [[output_dir, filtered_bams[i], min_c, t_dist, domi_ale_frq] for i in range(proc_num)]


    return utils.multi_map(processors, single_consensus_builder, params)


def output_trimmed_reads(trim_pos, bam_file):
    """
    This feature is optional.
    If chosen, trimmed reads extracted from filtered bams and re-direct to fastq files which
    are stored in outputs/trimmed_reads
    """

    in_samfile = pysam.AlignmentFile(bam_file, 'rb')

    reads_opt = bam_file.replace('.bam.sorted', '.fastq')

    out_fastq = open(reads_opt, 'w')

    left, right = trim_pos.split(':')

    for read in in_samfile.fetch():
        if read.is_reverse:
            read_name = read.query_name
            read_seq = str(Seq(read.query_alignment_sequence[int(left): -int(right)])\
                .reverse_complement())
            read_qual = read.qual[int(left): -int(right)][::-1]

            out_fastq.write('@'+read_name+'\n')
            out_fastq.write(read_seq+'\n')
            out_fastq.write("+\n")
            out_fastq.write(read_qual+'\n')
        else:
            read_name = read.query_name
            read_seq = str(Seq(read.query_alignment_sequence[int(left): -int(right)]))
            read_qual = read.qual[int(left): -int(right)]

            out_fastq.write('@'+read_name+'\n')
            out_fastq.write(read_seq+'\n')
            out_fastq.write("+\n")
            out_fastq.write(read_qual+'\n')
    out_fastq.close()


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

def blast_genomes(output_dir, db_dest, query_set, threads, if_clean):

    raw_blastn_opt_tab = get_raw_blastn_opt(output_dir, db_dest, query_set, threads)

    if os.path.exists(raw_blastn_opt_tab) and not if_clean:
        return raw_blastn_opt_tab

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

    gappy_idx = utils.find(homo_sub) # find gappy positions in subject sequence
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
    mydna = Seq(seq)
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
    query_genome = utils.unique([mp_dict[i[0]] for i in ipt_mtx])
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
            recon_seq = seq_proofreading(recon_genome_dict[c[0]])
            ctig_record = SeqRecord(Seq(recon_seq), id = c[0] + '_consensus', description = '')
            contigs_list.append(ctig_record)
        else:
            ctig_seq = len(ref_fna_dict[c[0]]) * '-'
            ctig_record = SeqRecord(Seq(ctig_seq), id = c[0] + '_consensus', description = '')
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
    seq_records.append(SeqRecord(Seq(ref_seq_concatenated), id = ref_fna.split('/')[-1], description= ''))
    for file in files_lst:
        file_to_dict = SeqIO.to_dict(SeqIO.parse(open(file), "fasta"))
        seq_name = file.split('/')[-1].split('____')[-1]
        concatenated_seq = ''
        for header in sorted_headers:
            header = header + '_consensus'
            concatenated_seq += str(file_to_dict[header].seq)
        concatenated_seq = seq_proofreading(concatenated_seq)
        single_seq_record = SeqRecord(Seq(concatenated_seq), id = seq_name, description = '')
        seq_records.append(single_seq_record)
    SeqIO.write(seq_records, output_file, 'fasta')

    return output_file

def abspath_finder(file):
    # Convert wahtever path into abs path
    return os.path.abspath(file)

def create_folder(name):
    folder_name = abspath_finder(name) # determine the abs path
    if os.path.exists(folder_name):
        logger.info('{} exists: pass!\n'.format(folder_name))
        pass
    else:
        os.makedirs(folder_name)
        logger.info('Creating folder: {}\n'.format(folder_name))
    return folder_name

def get_bowtie2_db_files(intermediate_path, reference_genome):

    opt_1 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.1.bt2'
    opt_2 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.2.bt2'
    opt_3 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.3.bt2'
    opt_4 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.4.bt2'
    opt_5 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.rev.1.bt2'
    opt_6 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.rev.2.bt2'
    opt_prefix = intermediate_path + '/' + reference_genome.split('/')[-1]

    return {'db_files': [opt_1, opt_2, opt_3, opt_4, opt_5, opt_6], 'db_prefix': opt_prefix}



def get_single_bam_file(sample_4_bowtie2, intermediate_path):
    opt_raw_bam = intermediate_path +'/'+ sample_4_bowtie2.split('/')[-1] + '.bam'
    return opt_raw_bam



def get_blastn_db_files(intermediate_path, reference_genome):

    opt_1 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.nhr'
    opt_2 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.nin'
    opt_3 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.nsq'
    opt_prefix = intermediate_path + '/' + reference_genome.split('/')[-1]

    return {'db_files': [opt_1, opt_2, opt_3], 'db_prefix': opt_prefix}


def get_raw_blastn_opt(intermediate_path, db_dest, query_set, threads):
    outfmt = '6 qaccver saccver pident length mismatch gapopen qstart qend\
                sstart send evalue bitscore qseq sseq'

    cmd = "blastn -db {} -query {} -outfmt '{}' -num_threads {} -word_size 9 -out {}/raw_blastn_opt.tab".format(db_dest, query_set, outfmt, threads, intermediate_path)

    return '{}/raw_blastn_opt.tab'.format(intermediate_path)

def seq_proofreading(Sequence):
    legit_bases = ['A', 'T', 'G', 'C', '-'] # five bases allowed to occur in the reconstructed sequence
    legalized_seq = [b.upper() if b.upper() in legit_bases else '-' for b in Sequence]
    return "".join(legalized_seq)


def authenticate(intermediate_path, refseq, filtered_bams):
    mapdamage2_results = create_folder(intermediate_path + '/mapdamage')


    G2A = []
    C2T = []

    for single_filetered_bam in filtered_bams:

        mp_opt_dir = mapdamage2_results + '/' + single_filetered_bam.split('/')[-1].replace('.bam', '_mp')
        sorted_bam =  single_filetered_bam + '.sorted'
        cmd = 'mapDamage -r {} -i {} -d {} --no-stats'.format(refseq, sorted_bam, mp_opt_dir)
        run_cmd_in_shell(cmd)
        G2A.append(mp_opt_dir + '/3pGtoA_freq.txt')
        C2T.append(mp_opt_dir + '/5pCtoT_freq.txt')

    return G2A, C2T


if __name__ == "__main__":

    main()
