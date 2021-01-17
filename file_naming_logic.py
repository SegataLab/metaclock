#!/usr/bin/env python


samples_4_bowtie2 = ['/shares/CIBIO-Storage/CM/news/users/kun.huang/tmp_Mac_test/wiki/MetaClock/ancient_metagenomes/ERR3003614',
                     '/shares/CIBIO-Storage/CM/news/users/kun.huang/tmp_Mac_test/wiki/MetaClock/ancient_metagenomes/ERR3003615',
                     '/shares/CIBIO-Storage/CM/news/users/kun.huang/tmp_Mac_test/wiki/MetaClock/ancient_metagenomes/ERR3003619']

samples_4_blastn = ['/shares/CIBIO-Storage/CM/news/users/kun.huang/tmp_Mac_test/wiki/MetaClock/modern_genome_assemblies/GCF_000529525.fna',
                    '/shares/CIBIO-Storage/CM/news/users/kun.huang/tmp_Mac_test/wiki/MetaClock/modern_genome_assemblies/GCF_900289035.fna',
                    '/shares/CIBIO-Storage/CM/news/users/kun.huang/tmp_Mac_test/wiki/MetaClock/modern_genome_assemblies/GCF_902384065.fna']

contigs_folder_path = '/shares/CIBIO-Storage/CM/news/users/kun.huang/tmp_Mac_test/wiki/MetaClock/modern_genome_assemblies'


intermediate_path = '/shares/CIBIO-Storage/CM/news/users/kun.huang/tmp_Mac_test/wiki/MetaClock/intermediates'
reference_genome = '/shares/CIBIO-Storage/CM/news/users/kun.huang/tmp_Mac_test/wiki/MetaClock/GCA_001639275.fna'


def building_bowtie2_db(intermediate_path, reference_genome):

	opt_1 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.1.bt2'
	opt_2 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.2.bt2'
	opt_3 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.3.bt2'
	opt_4 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.4.bt2'
	opt_5 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.rev.1.bt2'
	opt_6 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.rev.2.bt2'

	return opt_1, opt_2, opt_3, opt_4, opt_5, opt_6



def reads_2_bam(sample_4_bowtie2, intermediate_path):
	opt_raw_bam = intermediate_path +'/'+ sample_4_bowtie2.split('/')[-1] + '.bam'
	return opt_raw_bam


def filter_bam(raw_bam, intermediate_path):
	filter_bam = intermediate_path + '/' + 'filtered____'+ raw_bam.split('/')[-1]
	return filter_bam

def FilteredBam_2_consensus(filtered_bam, intermediate_path):
	opt_filtered_consensus = intermediate_path + '/' + 'consensus_'+ filtered_bam.split('/')[-1].replace('.bam', '.fna')
	sorted_filtered_bam = filtered_bam + '.sorted'
	indexed_sorted_bam = sorted_filtered_bam + '.bai'
	return opt_filtered_consensus, sorted_filtered_bam, indexed_sorted_bam


print("Generating intermediate files through bowtie2 process.......\n.\n.\n.\n")
bt2_db = building_bowtie2_db(intermediate_path, reference_genome)
print("bowtie2 database files:\n")
print(reference_genome + " ------> " + bt2_db[0])
print(reference_genome + " ------> " + bt2_db[1])
print(reference_genome + " ------> " + bt2_db[2])
print(reference_genome + " ------> " + bt2_db[3])
print(reference_genome + " ------> " + bt2_db[4])
print(reference_genome + " ------> " + bt2_db[5])
print('\n')
print('\n')
for s in samples_4_bowtie2:
	raw_bam = reads_2_bam(s, intermediate_path)
	print("The first intermediate file (raw bam file):\n{} -----> {}".format(s,raw_bam))
	filtered_bam = filter_bam(raw_bam, intermediate_path)
	print("The second intermediate file (filtered bam file):\n{} -----> {}".format(raw_bam, filtered_bam))
	filtered_consensus = FilteredBam_2_consensus(filtered_bam, intermediate_path)
	print("The third intermediate file1 (filtered consensus file):\n{} -----> {}".format(filtered_bam, filtered_consensus[0]))
	print("The third intermediate file2 (filtered sorted bam file):\n{} -----> {}".format(filtered_bam, filtered_consensus[1]))
	print("The third intermediate file3 (filtered sorted indexed bam file):\n{} -----> {}".format(filtered_bam, filtered_consensus[2]))

print('\n')
print('\n')
print('\n')

print('Generating intermediate files through blastn....\n.\n.\n.\n.\n ')


def building_blastn_db(intermediate_path, reference_genome):
	
	opt_1 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.nhr'
	opt_2 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.nin'
	opt_3 = intermediate_path + '/' + reference_genome.split('/')[-1] + '.nsq'
	opt_prefix = intermediate_path + '/' + reference_genome.split('/')[-1]

	return opt_1, opt_2, opt_3, opt_prefix

blastn_db = building_blastn_db(intermediate_path, reference_genome)
print("blastn database files:\n")
print(reference_genome + ' -----> '+ blastn_db[0])
print(reference_genome + ' -----> '+ blastn_db[1])
print(reference_genome + ' -----> '+ blastn_db[2])


def generate_query_set_and_mp_file(intermediate_path, contigs_folder):

	mapping_file = intermediate_path + '/genome_contigs_mp.csv'
	cmd = 'cat {}/* > {}/QuerySet.fna'.format(contigs_folder, intermediate_path)
	return mapping_file, '{}/QuerySet.fna'.format(intermediate_path)

mp_QuerySet = generate_query_set_and_mp_file(intermediate_path, contigs_folder_path)

print("Generate the mapping file to link genome name and corresponding contigs:\n{}".format(mp_QuerySet[0]))
print("Generate the collective set of all genomes for blastn:\n{}".format(mp_QuerySet[1]))


def generate_raw_blastn_results(intermediate_path, db_dest, query_set, threads):
	outfmt = '6 qaccver saccver pident length mismatch gapopen qstart qend\
				sstart send evalue bitscore qseq sseq'

	cmd = "blastn -db {} -query {} -outfmt '{}' -num_threads {} -word_size 9 -out {}/raw_blastn_opt.tab".format(db_dest, query_set, outfmt, threads, intermediate_path)

	return '{}/raw_blastn_opt.tab'.format(intermediate_path)

raw_blastn_opt = generate_raw_blastn_results(intermediate_path, blastn_db[3], mp_QuerySet[1], 10)
print("Generate the raw blastn results:\n{}".format(raw_blastn_opt))

def QC_on_blastn(raw_blastn_opt, length, pid, intermediate_path):
	cmd = 'cut -f 1-14 {} | bo6_screen.py --length {} --pid {} > {}/cleaned_blastn_opt.tab'.format(raw_blastn_opt, length, pid, intermediate_path)
	return '{}/cleaned_blastn_opt.tab'.format(intermediate_path)

cleaned_blastn_opt = QC_on_blastn(raw_blastn_opt, 500, 95, intermediate_path)	
print("Generate cleaned blastn output:\n{} -----> {}".format(raw_blastn_opt, cleaned_blastn_opt))
print("A long process of puzzling up cleaned blastn results......\n.\n.\n.\n")

def generate_blastn_genoems(sample, intermediate_path):
	opt = intermediate_path + '/' + 'homolog____' + sample.split('/')[-1]
	return opt

for s in samples_4_blastn:
	blastn_final_opt = generate_blastn_genoems(s, intermediate_path)
	print(s + ' ----> ' + blastn_final_opt)


