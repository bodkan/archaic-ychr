import argparse
import sys
from collections import defaultdict
import vcf


def valid_sample_names(samples_given, samples_in_header):
    '''Check if all sample names provided at the command line are
    present in the VCF header.
    '''
    return all(s in samples_in_header for s in samples_given)


parser = argparse.ArgumentParser(description='Convert genotypes from a specified set of samples to FASTA')
parser.add_argument('--vcf-file', help='VCF file to parse', required=True)
parser.add_argument('--fasta-file', help='FASTA output file', required=True)
parser.add_argument('--sample-names', help='List of samples to take from VCF', nargs='*', default=[])
parser.add_argument('--chrom', help='VCF file to parse', default='Y')
args = parser.parse_args()

vcf_reader = vcf.Reader(open(args.vcf_file, 'rb'))

# check if the provided samples are indeed present in the VCF
if not valid_sample_names(args.sample_names, vcf_reader.samples):
    sys.exit('Not all provided sample names are present in the VCF file!')

# if the sample names were not provided, use all samples in the VCF
if not args.sample_names:
    args.sample_names = vcf_reader.samples

# iterate through the VCF and accumulate called bases for all samples
samples_dict = defaultdict(list)
for record in vcf_reader.fetch(args.chrom):
    for sample in args.sample_names:
        call = record.genotype(sample)
        called_base = call.gt_bases if call.gt_bases else 'N'
        samples_dict[call.sample].append(called_base)

# write out the called bases for each sample in a FASTA format
with open(args.fasta_file, 'w') as output_fasta:
    for sample_name in args.sample_names:
        print('>' + sample_name, file=output_fasta)
        print(''.join(samples_dict[sample_name]), file=output_fasta)
        print('', file=output_fasta)
