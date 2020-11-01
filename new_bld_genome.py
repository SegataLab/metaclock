def main():
    # parse parameters
    pass

    inter_results = []
    # deal with working directory
    for params in workflows:
        dest = workflow(params)
        inter_results.append(dest)

    # merge fiiles in inter_results
    merge(dest)


def workflow(params):
    mode = params['mode']
    # build database
    db_dest = build_db_file(mode, ref_genome, db_dir)

    # build mapping with database files
    inter_file = build_mapping(mode, db_dest, param_set)

    # workflow 1
    if mode == 'reads':
        pass

    # workflow 2
    if mode == 'contigs':
        pass

    return file_dir


def build_db_file(mode, ref_genome, db_dir):
    """
    Args:
      ref_genome: input file path
      db_dir: destination directory path

    Return:
      destination files path
    """
    params = {
        'ref': ref_genome,
        'dest': db_dir + '/' + ref_genome,
    }
    if mode == 'reads':
        # build command
        cmd = "bowtie2-build {ref} {dest}".format(**params)
    elif mode == 'contigs':
        cmd = "makeblastdb -in {ref} -dbtype nucl -out {dest}".fromat(**params)

    # run the command
    # TODO: exception handling
    run_cmd_in_shell(cmd)

    return params['dest']


def build_mapping(mode, db_dest, param_set):
    if mode == reads:
        bwt2_batch_mapping()
    elif mode == 'configs':
        blast_genomes()

    return dest_dir


def single_mapping(multi_args):
    ref, s, t, k = multi_args
    k = k.replace("*", "")
    suffix = detect_reads_suffix(s)
    labeled_modern_opt_bam = '/'.join(s.split('/')[:-1])+'/m__'+ s.split('/')[-1] + '.bam'
    if suffix == "bz2":
        cmd = 'bzcat {}/*bz2 | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), labeled_modern_opt_bam)
        run_cmd_in_shell(cmd)
    elif suffix == "gz":
        cmd = 'zcat {}/*gz | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), labeled_modern_opt_bam)
        run_cmd_in_shell(cmd)
    elif suffix == "fastq":
        cmd = 'cat {}/*fastq | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), labeled_modern_opt_bam)
        run_cmd_in_shell(cmd)
    else:
        sys.exit("Reads have to be in the form of .bz2, .gz or .fastq!")

def bwt2_batch_mapping(age_type, ):
    """
    Args:
      args_type (int): 1 -> ancient, 2 -> modern
    """

    if ',' in samples:
        sample_list = [os.getcwd() + '/' +  i for i in samples]
    elif '.txt' in samples:
        sample_list = [os.getcwd() + '/' +  i.rstrip() for i in open(samples).readlines()]
    else:
        sample_list = [os.getcwd() + '/' + samples]

    threads = [threads] * len(sample_list)
    ref_lst = [ref_genome] * len(sample_list)
    mode = [search_mode] * len(sample_list)

    multi_map(processors, single_mapping, zip(ref_lst, sample_list, threads, mode))

def blast_genomes():

    outfmt = '6 qaccver saccver pident length mismatch gapopen qstart qend\
                    sstart send evalue bitscore qseq sseq'
    cmd = "blastn -db {} -query {} -outfmt '{}' -num_threads {}\
 			-word_size 9 -out {}/INTER_blast_opt_tmp.tab".format(ref_fna, query, outfmt, threads, opt_dir)

    run_cmd_in_shell(cmd)

def run_cmd_in_shell(cmd):
    subprocess.call(cmd, shell=True)
