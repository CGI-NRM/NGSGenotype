# Usage:
# python ultimate_locus_parser.py [folder containing fastq files] > example_output_name.csv

import os
import gzip
import numpy
import sys

# Functions:
def calculate_allele_dict(file_name):
    next_line = ''
    allele_dict = {}
    with gzip.open(data_folder + '/' + file_name) as locus_fastq:
        for line in locus_fastq:
            if next_line == 'yes':
                str_line = line.decode().strip()
                if str_line not in allele_dict.keys():
                    allele_dict[str_line] = 1
                else:
                    allele_dict[str_line] = allele_dict[str_line] + 1
                next_line = ''
            if line.startswith(b'@') and len(line) > 2:
                next_line = 'yes'
        return allele_dict

def determine_genotype(locus_dict):
    if '' in locus_dict.keys():
        locus_dict.pop('') # remove value with key '' (nothing)
    all_alleles = locus_dict.keys()
    all_counts = [locus_dict[i] for i in all_alleles]
    sorted_alleles = [list(all_alleles)[i] for i in numpy.argsort(all_counts)][::-1]
    if len(sorted_alleles) > 1:
        allele_sum = locus_dict[sorted_alleles[0]] + locus_dict[sorted_alleles[1]]
    elif len(sorted_alleles) == 1:
        allele_sum = locus_dict[sorted_alleles[0]]
    else:
        allele_sum = 0
    if allele_sum >= min_reads:
        if len(sorted_alleles) > 1: # also check if no allele is one of the top ones
            if ((locus_dict[sorted_alleles[1]] / allele_sum) * 100) >= min_het_perc:
                top_alleles = sorted(sorted_alleles[0:2]) # sort top alleles alphabetically to make order consistent
                return top_alleles[0] + top_alleles[1]
            else:
                return sorted_alleles[0] + sorted_alleles[0]
        elif len(sorted_alleles) == 0:
            return ''
        elif sorted_alleles[0] != '':
            return sorted_alleles[0] + sorted_alleles[0]
        else:
            return ''
    else:
        return ''

def recurring_csv(gt_dict, print_order_list, out_string):
    if out_string != 'EMPTY':
        out_string = out_string + ',' + gt_dict[print_order_list.pop(0)]
    else:
        out_string = gt_dict[print_order_list.pop(0)]
    if len(print_order_list) == 0:
        return out_string
    else:
        return recurring_csv(gt_dict, print_order_list, out_string)

# Settings:
min_reads = 10
min_het_perc = 20
data_folder = sys.argv[1] # folder containing fastq files

# Get names of files and samples:
all_fastq_files = [i for i in os.listdir(data_folder) if i.endswith('fastq.gz')]
all_samples = list(set([i.split('_')[2].split('-')[0] for i in all_fastq_files])) # this will summarise samples with the same names (happens in Becky was in a prank mood while naming samples)

# Load loci print order:
full_loci_list = []
with open('../bear_snp_loci.csv') as bear_snp_csv:
    for locus in bear_snp_csv:
        full_loci_list.append(locus.strip())

# Calculate and print all genotypes:
print(','.join(['Sample'] + full_loci_list))
for i in all_samples:
    sample_dict = {}
    for j in all_fastq_files:
        if i in j:
            sample_dict[j.split('_')[0]] = calculate_allele_dict(j)
    sample_genotypes = {}
    for cur_locus_key in sample_dict.keys():
        sample_genotypes[cur_locus_key] = determine_genotype(sample_dict[cur_locus_key])
    print(i.split('BMK')[0] + ',' + recurring_csv(sample_genotypes, full_loci_list.copy(), 'EMPTY'))
