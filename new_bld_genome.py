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
        cmd = 'bzcat {}/*fastq.bz2 | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), labeled_modern_opt_bam)
        run_cmd_in_shell(cmd)
    elif suffix == "gz":
        cmd = 'zcat {}/*fastq.gz | bowtie2 -x {} -p {} --end-to-end {} --no-unal -U - -S - | samtools view -bS - > {}'.format(s, ref, t, ' '.join(k.split(',')), labeled_modern_opt_bam)
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

def QC_on_blastn(blast_tab, length, pid, opt_dir = os.getcwd()):

    """
    It applies QC on blast output
    Input: 1)'INTER_blast_opt_tmp.tab', 2) length of hits, 3) identity percentage of hits
    Program: call script 'bo6_screen.py'
    Output: filtered blast output, 'INTER_blastn_opt_tmp_cleaned.tab' 
    """

    cmd = 'cut -f 1-14 {} | bo6_screen.py --length {} --pid {} > {}/INTER_blastn_opt_tmp_cleaned.tab'.format(blast_tab, length, pid, opt_dir)
    subprocess.call(cmd, shell = True)


def run_cmd_in_shell(cmd):
    subprocess.call(cmd, shell=True)

def detect_reads_suffix(reads_foler):
    """
    Args: reads folder
    Return: suffix of reads file inside the folder
    """
    reads_file = subprocess.getoutput('ls {}/*fastq*'.format(reads_foler)).split('\n')[0]
    return reads_file.split('.')[-1]

def single_bam_filter(multi_args):
    """
    Args: a given bam coupled with parameters for QC
    Return: filtered bam 
    """

    bam, mq, ml, m_snp = multi_args
    filter_path = os.path.dirname(__file__)+"/cmseq/filter.py"
    cmd = "samtools view -h -F 0x4 {} | python3 {} --minqual {} --minlen {} --maxsnps {} > {}".format(bam, filter_path,\
     mq, ml, m_snp, "filtered_INTER_"+bam)
    subprocess.call(cmd, shell =True)

def batch_bam_filter(args):
    """
    Parallelize single_bam_filter()
    """

    bams = subprocess.getoutput("ls *.bam").split('\n')
    mq_lst = [args.minimum_quality] * len(bams)
    ml_lst = [args.minimum_length] * len(bams)
    m_snp = [args.max_snps] * len(bams)

    multi_map(args.processor, single_bam_filter, zip(bams, mq_lst, ml_lst, m_snp))

def single_consensus_builder(multi_args):
    """
    Args: a filtered bam coupled with parameters for QC
    return: consensus sequences in a fasta file
    """

    f_bam, mc, trim_end, domall = multi_args
    consensus_path = os.path.dirname(__file__)+"/cmseq/consensus.py"

    if trim_end != None:

        cmd = "python3 {} {} --sortindex --mincov {} --trim {} --dominant_frq_thrsh {} > {}".format(consensus_path, f_bam, mc, trim_end, domall, f_bam.split(".")[0]+'_consensus.fna')
        subprocess.call(cmd, shell = True)
    else:
        cmd = "python3 {} {} --sortindex --mincov {} --dominant_frq_thrsh {} > {}".format(consensus_path, f_bam, mc, domall, f_bam.split(".")[0]+'_consensus.fna')
        subprocess.call(cmd, shell = True)

def batch_consensus_builder(args):

    """
    Parallelize single_consensus_builder()
    """

    f_bams = subprocess.getoutput("ls filtered_INTER_*bam").split('\n')
    mc_lst = [args.minimum_coverage] * len(f_bams)
    trim_lst = [args.trim_reads_end] * len(f_bams)
    domall_lst = [args.minimum_dominant_allele] * len(f_bams)
    multi_map(args.processor, single_consensus_builder, zip(f_bams, mc_lst, trim_lst, domall_lst))

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


###############################################################################################################
# So far almost all work was handled by bowtie2, blastn, and other dependencies
###############################################################################################################


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
    

##############################################################################################################
# Output reconstructed genomes in individual fasta files same as how consensus works
##############################################################################################################