#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict
import vcf
import pandas as pd


def valid_sample_names(samples_given, samples_in_header):
    """Check if all sample names provided at the command line are
    present in the VCF header.
    """
    return all(s in samples_in_header for s in samples_given)


parser = argparse.ArgumentParser(description="Convert genotypes from VCF to FASTA")
parser.add_argument("--vcf", help="VCF file to parse", required=True)
parser.add_argument("--fasta", help="FASTA output file", required=True)
parser.add_argument("--include", help="List of samples to include from the VCF", nargs="*", default=[])
parser.add_argument("--exclude", help="List of samples to exclude from the VCF", nargs="*", default=[])
parser.add_argument("--chimp-unique", help="Include sites with unique chimp alleles?", action="store_true", default=False)
args = parser.parse_args()
# args = parser.parse_args("--vcf data/vcf/merged_full.vcf.gz --fasta asd.fa".split())

vcf_reader = vcf.Reader(open(args.vcf, "rb"))

args.include, args.exclude = set(args.include), set(args.exclude)
# check if the provided samples are indeed present in the VCF
if not valid_sample_names(args.include, vcf_reader.samples):
    sys.exit("Not all samples to include are present in the VCF file!")
if not valid_sample_names(args.exclude, vcf_reader.samples):
    sys.exit("Not all samples to exclude are present in the VCF file!")
if len(args.include & args.exclude):
    sys.exit("Overlap between samples to include and exclude not allowed")

# if the sample names were not provided, use all samples in the VCF
samples = args.include if args.include else set(vcf_reader.samples) - args.exclude

# iterate through the VCF and accumulate called bases for all samples
samples_dict = defaultdict(list)
for record in vcf_reader.fetch("Y"):
    for name in samples:
        call = record.genotype(name)
        called_base = call.gt_bases if call.gt_bases else "N"
        samples_dict[call.sample].append(called_base)
breakpoint()

gt_df = pd.DataFrame(samples_dict)
if args.chimp_unique: drop = gt_df.columns.str.extract("(chimp.*)", flags=re.IGNORECASE, expand=False).dropna()
allele_counts = gt_df.drop(columns=drop).apply(lambda row: len(set(i for i in row if i != "N")), axis=1)
gt_df = gt_df.loc[allele_counts > 1]

# write out the called bases for each sample in a FASTA format
with open(args.fasta, "w") as output:
    for name in samples:
        print(">" + name, file=output)
        print("".join(gt_df[name]), file=output)

