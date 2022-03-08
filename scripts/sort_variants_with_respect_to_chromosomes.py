#!/usr/bin/env python

from mimetypes import guess_type
from sys import stderr
from sys import stdout
from collections import defaultdict
import argparse
import gzip

def sort_variant(scf, line_to_sort):
    if scf2chr[scf] == 'A':
        variant_file_A.write(line_to_sort)
    elif scf2chr[scf] == 'X':
        variant_file_X.write(line_to_sort)
    else:
        variant_file_other.write(line_to_sort)

parser = argparse.ArgumentParser(description='Sorting variants by chromosomal assignment')
parser.add_argument('vcf_file', help='The vcf file to be sorted')
parser.add_argument('asignment_table', help='tsv file with header that contain "scf" and "chr" columns')
parser.add_argument('-o', '-output', help='pattern used to generate new vcf files (with _<Chr>.vcf suffix)', default = 'sorted_variants')

args = parser.parse_args()

encoding = guess_type(args.vcf_file)[1]
_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

stderr.write("Sorting " + args.vcf_file + " vcf file using " + args.asignment_table + " asignment table.")

scf2chr = defaultdict(lambda: 'other')

# args.asignment_table='data/Afus1/asm/asm_01_spades_filt_reads/scaffold_lengths_and_covs.tsv'
with open(args.asignment_table) as asn_tab:
    header = asn_tab.readline().rstrip('\n').split('\t')
    scf_col = [i for i,h in enumerate(header) if h == 'scf'][0]
    chr_col = [i for i,h in enumerate(header) if h == 'chr'][0]
    # if there is no scf/chr, it prints an obscure error
    for line in asn_tab:
        scf_info_vec = line.rstrip('\n').split('\t')
        # this could be more general
        if scf_info_vec[chr_col] in ['A', 'X']:
            scf2chr[scf_info_vec[scf_col]] = scf_info_vec[chr_col]

### sort header
# args.o = 'data/resequencing/SNP_calls/freebayes_filtered_sorted'
variant_file_A = open(args.o + "_A.vcf", "w")
variant_file_X = open(args.o + "_X.vcf", "w")
variant_file_other = open(args.o + "_other.vcf", "w")

# args.vcf_file='data/resequencing/SNP_calls/freebayes_all_samples_filt.vcf'
with _open(args.vcf_file) as vcf_file:
    header_line = vcf_file.readline()
    while header_line.startswith("#"):
        if header_line.startswith('##contig'):
            scf=header_line.split('ID=')[1].split(',')[0]
            # lines that should be sorted to individual headers
            sort_variant(scf, header_line)
        else:
            # lines that should be present in all vcf headers
            variant_file_A.write(header_line)
            variant_file_X.write(header_line)
            variant_file_other.write(header_line)
        header_line = vcf_file.readline()

    variant = header_line
    scf = variant.rstrip('\n').split('\t')[0]
    sort_variant(scf, variant)

    for variant in vcf_file:
        scf = variant.rstrip('\n').split('\t')[0]
        sort_variant(scf, variant)

variant_file_A.close()
variant_file_X.close()
variant_file_other.close()
