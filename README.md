# Ittiada - Inferring Time-Tree integrating Ancient DNA and Assemblies #

**Main features**

* Take assembled genomes (including metagenomically-assembled genomes) or NGS reads as input
* Use single genome or multiple genomes as reference
* Reduce noises from damaged ancient DNA by trimming highly frequent degenerated sites
* Support multiprocessing
* Maximize traceble evolutionary signal searching for all possible homologous genomic regions (more than coding sequences)
* Automatedly merging ancient genomes with modern ones

**Dependencies:**

* cmseq
* bo6_screen.py
* numpy
* samtools
* pysam
* biopython
--------------------------------
## Description

Evolution of bateria plays a pivotal role in shaping human gut microbiome.
Phylogenetic analysis has been recognized as the most powerful molecular approach for recontructing bacterial evolutionary histories.
With the advent of computational metagenomics, well-assembled bacterial genomes have been uncovered at an exponential rate. Therefore, an automated and efficient approach is required to build strain-level phologenies of
massive number of bacterial species. Here, we introduce Ittiada - Inferring Time-Tree integrating Ancient DNA and Assemblies - an automated tool for building multiple sequence alignment at the genome level,
particularly for accurately integrating ancient lineages with modern relatives. The genome alignment is used for downstrean
phylogentics analysis, such as molecular clock and highly-resolved strain-level phylogeny. 

Graphic pipeline for contigs-based approach:
![Fig. 1](https://bitbucket.org/CibioCM/ittiada/raw/1a25119c4c0ba42fb8206ee267b8a2294690b57d/images/Contigs_based.png "Fig. 1")

Graphic pipeline for reads-based approach:
![Fig. 2](https://bitbucket.org/CibioCM/ittiada/raw/1a25119c4c0ba42fb8206ee267b8a2294690b57d/images/Reads_based.png "Fig. 2")

```
usage: contigs_based, reads_based [-h] {contigs_based,reads_based} ...

positional arguments:
  {contigs_based,reads_based}
                        program mode
    contigs_based       Use contigs_based approach.
    reads_based         Use reads_based approach.

optional arguments:
  -h, --help            show this help message and exit
  
usage: contigs_based, reads_based contigs_based [-h] [-trim TRIM_READS_END]
                                                [-minqual MINIMUM_QUALITY]
                                                [-minlen MINIMUM_LENGTH]
                                                [-maxsnps MAX_SNPS]
                                                [-mincov MINIMUM_COVERAGE]
                                                [-m SEARCH_MODE]
                                                [-aln_len ALIGNMENT_LENGTH]
                                                [-pid IDENTITY]
                                                [-blast_t BLAST_THREADS]
                                                [-t THREADS] [-p PROCESSOR]
                                                [AncientSamples]
                                                [contigs_isorefs_mags]
                                                [RefSeq_fasta]

positional arguments:
  AncientSamples        input folder names contains ancient microbiome reads
                        separated by comma or input a text file each line
                        indicates a folder name
  contigs_isorefs_mags  specify the folder name which contains all query
                        genomes (either RefSeqs or MAGs), in FASTA form, used
                        for blast.
  RefSeq_fasta          input fasta file of reference sequence

optional arguments:
  -h, --help            show this help message and exit
  -trim TRIM_READS_END, --trim_reads_end TRIM_READS_END
                        Trim the reads before computing the consensus. A value
                        of 10:10 means that the first and last 10 positions of
                        each read will be ignored. Default: None
  -minqual MINIMUM_QUALITY, --minimum_quality MINIMUM_QUALITY
                        specify the minimum quality of aligned reads to keep.
  -minlen MINIMUM_LENGTH, --minimum_length MINIMUM_LENGTH
                        specify the minimum length of aligned reads to keep.
  -maxsnps MAX_SNPS, --max_snps MAX_SNPS
                        maximum edit distance on the alignment for a read to
                        pass.
  -mincov MINIMUM_COVERAGE, --minimum_coverage MINIMUM_COVERAGE
                        minimum position coverage to keep for reconstruction.
  -m SEARCH_MODE, --search_mode SEARCH_MODE
                        specify search mode in bowtie2. e.g. -a or -k,3. Look
                        for bowtie2 manual for details
  -aln_len ALIGNMENT_LENGTH, --alignment_length ALIGNMENT_LENGTH
                        minimum alignment length.
  -pid IDENTITY, --identity IDENTITY
                        minimum identity of alignment.
  -blast_t BLAST_THREADS, --blast_threads BLAST_THREADS
                        threading number of blastn
  -t THREADS, --threads THREADS
                        threading number of bowtie2.
  -p PROCESSOR, --processor PROCESSOR
                        multiple processor number. Suggest input similar
                        number as samples for the maximum computation trade-
                        off.

usage: contigs_based, reads_based reads_based [-h] [-trim TRIM_READS_END]
                                              [-a_minqual A_MINIMUM_QUALITY]
                                              [-a_minlen A_MINIMUM_LENGTH]
                                              [-a_maxsnps A_MAX_SNPS]
                                              [-a_mincov A_MINIMUM_COVERAGE]
                                              [-a_bt_t A_THREADS]
                                              [-a_m A_SEARCH_MODE]
                                              [-m_minqual M_MINIMUM_QUALITY]
                                              [-m_minlen M_MINIMUM_LENGTH]
                                              [-m_maxsnps M_MAX_SNPS]
                                              [-m_mincov M_MINIMUM_COVERAGE]
                                              [-m_bt_t M_THREADS]
                                              [-m_m M_SEARCH_MODE]
                                              [-a_p A_PROCESSOR]
                                              [-m_p M_PROCESSOR] [-tr]
                                              [RefSeq_fasta] [AncientSamples]
                                              [ModernSamples]

positional arguments:
  RefSeq_fasta          input a fasta file of reference sequence
  AncientSamples        input folder names containing metagenomic reads
                        separated by comma or input a text file each line
                        indicates a folder name
  ModernSamples         input folder names containing metagenomic reads
                        separated by comma or input a text file each line
                        indicates a folder name

optional arguments:
  -h, --help            show this help message and exit
  -trim TRIM_READS_END, --trim_reads_end TRIM_READS_END
                        Trim the reads before computing the consensus. A value
                        of 10:10 means that the first and last 10 positions of
                        each read will be ignored. Default: None
  -a_minqual A_MINIMUM_QUALITY, --a_minimum_quality A_MINIMUM_QUALITY
                        specify the minimum quality of aligned reads to keep
                        [default:30].(for ancient samples).
  -a_minlen A_MINIMUM_LENGTH, --a_minimum_length A_MINIMUM_LENGTH
                        specify the minimum length of aligned reads to keep
                        [default:30].(for ancient samples)
  -a_maxsnps A_MAX_SNPS, --a_max_snps A_MAX_SNPS
                        maximum edit distance on the alignment for a read to
                        pass [default:0.03].(for ancient samples)
  -a_mincov A_MINIMUM_COVERAGE, --a_minimum_coverage A_MINIMUM_COVERAGE
                        minimum position coverage to keep for reconstruction
                        [default:5].(for ancient samples)
  -a_bt_t A_THREADS, --a_threads A_THREADS
                        threading number of bowtie2 set for ancient
                        samples.[default:1]
  -a_m A_SEARCH_MODE, --a_search_mode A_SEARCH_MODE
                        specify search mode in bowtie2. e.g. -a or -k,3. Look
                        for bowtie2 manual for details [default:-k 3].(for
                        ancient samples)
  -m_minqual M_MINIMUM_QUALITY, --m_minimum_quality M_MINIMUM_QUALITY
                        specify the minimum quality of aligned reads to keep
                        [default:30].(for modern samples).
  -m_minlen M_MINIMUM_LENGTH, --m_minimum_length M_MINIMUM_LENGTH
                        specify the minimum length of aligned reads to keep
                        [default:30].(for modern samples)
  -m_maxsnps M_MAX_SNPS, --m_max_snps M_MAX_SNPS
                        maximum edit distance on the alignment for a read to
                        pass [default:0.03].(for modern samples)
  -m_mincov M_MINIMUM_COVERAGE, --m_minimum_coverage M_MINIMUM_COVERAGE
                        minimum position coverage to keep for reconstruction
                        [default:5].(for modern samples)
  -m_bt_t M_THREADS, --m_threads M_THREADS
                        threading number of bowtie2 set for modern
                        samples.[default:1]
  -m_m M_SEARCH_MODE, --m_search_mode M_SEARCH_MODE
                        specify search mode in bowtie2. e.g. -a or -k,3. Look
                        for bowtie2 manual for details [default:-k 3].(for
                        modern samples)
  -a_p A_PROCESSOR, --a_processor A_PROCESSOR
                        multiple processor number. *the production of threads
                        and processor should be less than available
                        CPUs[default:1]
  -m_p M_PROCESSOR, --m_processor M_PROCESSOR
                        multiple processor number. *the production of threads
                        and processor should be less than available
                        CPUs[default:1]
  -tr, --trimmed_reads  generate a file contains reads with possible
                        degenerated sites trimmed
```

## Example 1 (contigs_based mode): ##
```
Input: folder - 2179a_merged (Merged reads from ERR3278761 and ERR3278819 on ENA database)
       folder - 2180a_merged (Merged reads from ERR3278820 and ERR3278821 on ENA database)
       folder - 224_Iceman (ERR3294837 on ENA database)
       file - GCA_002834165.fna (Ruminococcus bromii RefSeq GCA_002834165.1 on NCBI)
       folder - genome_contigs (A folder contains a collection of assembled genomes, each fasta file represents one genome)
```

### Build a genome alignment, using assembled contigs and reads of ancient stool metagenomes (2179a_merged, 2180a_merged, 224_Iceman) with default options.

* Without trimming two ends reads where DNA degradation occurs frequently
* Minimum alignment quality [30] (i.e. aligned reads are most likely unique)
* Maximum edit distance on the alignment [0.03] (The edit distance can be defined as the number of single letter (nucleotide) changes that have to be made to one string (read) for it to be equal to another string (reference genome))
* Minimum alignment length [50nt].
* Minimum coverage [5X]
* Search mode in bowtie2 [-k,3] (i.e. search for 3 alignments, report each)
* Minimum BLASTn hit length [500bp]
* Minimum BLASTn hit identity [95.0]
* CPU used for blast phase [1]
* Thread used for bowtie2 [1]
* Processor used for reconstructing ancient genome [1]


```
BuildGenomeAln.py contigs_based 2179a_merged,2180a_merged,224_Iceman genome_contigs GCA_002834165.fna
```

### Output files


```
Parameters_setting.txt - Report for parameter setting in analysis.
GCA_002834165-GenomeAln_ContigsBased.fna - Genome alignment output [RefSeq-GenomeAln_Approach.fna]
intermediates - A folder of intermediate files.
aln_stats.tsv - Basic statistics of reconstructed alignment.
```

## Example 2(reads_based mode): ##

```
Input: folder - 2179a_merged (Merged reads from ERR3278761 and ERR3278819 on ENA database)
       folder - 2180a_merged (Merged reads from ERR3278820 and ERR3278821 on ENA database)
       folder - 224_Iceman (ERR3294837 on ENA database)
       file - GCA_002834165.fna (Ruminococcus bromii RefSeq GCA_002834165.1 on NCBI)
       folder - AsnicarF_2017__MV_FEI5_t3Q1 (metagenome NGS reads compressed in bz2 inside)
       folder - BritoIL_2016__M2.46.ST (metagenome NGS reads compressed in bz2 inside)
       folder - ChengpingW_2017__AS106raw (metagenome NGS reads compressed in bz2 inside)
       folder - CosteaPI_2017__mickey1-11-30-0 (metagenome NGS reads compressed in bz2 inside)
       folder - IjazUZ_2017__S153_a_WGS (metagenome NGS reads compressed in bz2 inside)
       file - metagenomes.txt (a text file which contains metagenome folder names provided)
```

metagenome.txt should look like:

```
AsnicarF_2017__MV_FEI5_t3Q1
BritoIL_2016__M2.46.ST
ChengpingW_2017__AS106raw
CosteaPI_2017__mickey1-11-30-0
IjazUZ_2017__S153_a_WGS
```

### Build a genome alignment, using metagenome samples and ancient stool metagenomes (2179a_merged, 2180a_merged, 224_Iceman) with default options.

* Without trimming two ends of reads where DNA degradation occurs frequently in ancient samples
* Minimum alignment quality for ancient samples [30] (i.e. aligned reads are most likely unique)
* Minimum alignment quality for modern samples [30] (i.e. aligned reads are most likely unique)
* Maximum edit distance on the alignment for ancient samples [0.03] (The edit distance can be defined as the number of single letter (nucleotide) changes that have to be made to one string (read) for it to be equal to another string (reference genome))
* Maximum edit distance on the alignment for modern samples [0.03] (The edit distance can be defined as the number of single letter (nucleotide) changes that have to be made to one string (read) for it to be equal to another string (reference genome))
* Minimum alignment length for ancient samples [50nt].
* Minimum alignment length for modern samples [50nt].
* Minimum coverage for ancient samples [5X]
* Minimum coverage for modern samples [5X]
* Search mode in bowtie2 [-k,3] for ancient samples (i.e. search for 3 alignments, report each)
* Search mode in bowtie2 [-k,3] for modern samples (i.e. search for 3 alignments, report each)
* Thread used for bowtie2 for ancient samples [1]
* Thread used for bowtie2 for modern samples [1]
* Processor used for processing ancient samples [1]
* Processor used for processing modern samples [1]

```
BuildGenomeAln.py reads_based GCA_002834165.fna 2179a_merged,2180a_merged,224_Iceman metagenomes.txt 
```

### Output files

```
Parameters_setting.txt - Report for parameter setting in analysis.
GCA_002834165-GenomeAln_ReadsBased.fna - Genome alignment output [RefSeq-GenomeAln_Approach.fna]
intermediates - A folder of intermediate files.
aln_stats.tsv - Basic statistics of reconstructed alignment.
```


