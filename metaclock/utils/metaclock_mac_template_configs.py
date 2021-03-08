#!/usr/bin/env python

import json
import argparse
import os, sys
def main():
    configs = """ 
    {
      "ancient_reads": {
      "input_type":"reads",
      "reference_genome":"",
      "age_type":1,
      "intermediate":"",
      "samples":"",
      "parameter_set":{
          "search_report_mode":"-k,5",
          "bowtie2_threads":3,
          "minimum_mapping_quality":30,
          "minimum_mapping_length":30,
          "maximum_snp_edit_distance":0.03,
          "nproc":3,
          "minimum_coverage":5,
          "trim_distance":"5:5",
          "dominant_allele_frequency":0.8,
          "output_trimmed_reads":0
        }
      },
      "modern_reads": {
      "input_type":"reads",
      "reference_genome":"",
      "age_type":2,
      "intermediate":"",
      "samples":"",
      "parameter_set":{
          "search_report_mode":"-k,1",
          "bowtie2_threads":3,
          "minimum_mapping_quality":30,
          "minimum_mapping_length":30,
          "maximum_snp_edit_distance":0.03,
          "nproc":3,
          "minimum_coverage":5,
          "dominant_allele_frequency":0.8
        }
      },
      "contigs": {
      "input_type":"contigs",
      "reference_genome":"",
      "intermediate":"",
      "samples":"",
          "parameter_set":{
          "homolog_length":500,
          "homolog_identity":95.0,
          "blastn_threads":6
        }
      }
    }
    """
    parser = argparse.ArgumentParser(description="This script will generate a template configuration file for metaclock_mac. And parameters can be tuned in users' need.")
    parser.add_argument('-d', '--directory', help='Specify the path of directory for storing configuration file. default: current working directory.', default=os.getcwd())

    args = parser.parse_args()
    configs_file_name = os.path.abspath(args.directory) + '/configs.json' 
    opt = open(configs_file_name, 'w')
    opt.write(configs)
    opt.close()

    sys.stdout.write("The template file has been generated in: {}\n".format(configs_file_name))
    sys.stdout.write("The template file can be opened and edited using text editors.\n")
    sys.stdout.write("Welcome back!\n")

if __name__ == "__main__":
  main()
