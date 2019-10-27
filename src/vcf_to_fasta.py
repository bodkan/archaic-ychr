#!/usr/bin/env python3

import re
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
parser.add_argument("--vcf", help="Input VCF file with biallelic variants", required=True)
parser.add_argument("--fasta", help="FASTA output file", required=True)
parser.add_argument("--include", help="List of samples to include from the VCF", nargs="*", default=[])
parser.add_argument("--exclude", help="List of samples to exclude from the VCF", nargs="*", default=[])
parser.add_argument("--variable", help="Output only variable sites?", action="store_true", default=False)
parser.add_argument("--no-damage", help="Remove potentional aDNA damage SNPs?", action="store_true", default=False)

args = parser.parse_args()
#args = parser.parse_args("--vcf data/vcf/exome_modern.vcf.gz --fasta out.fa --variable".split())

vcf_reader = vcf.Reader(open(args.vcf, "rb"))

args.include, args.exclude = set(args.include), set(args.exclude)
# check if the specified samples are present in the VCF
if not valid_sample_names(args.include, vcf_reader.samples):
    sys.exit("Not all samples to include are present in the VCF file!")
if len(args.include & args.exclude):
    sys.exit("Overlap between samples to include and exclude not allowed")

# if the sample names were not provided, use all samples in the VCF
samples = args.include if args.include else set(vcf_reader.samples) - args.exclude

# iterate through the VCF and accumulate called bases for all samples
samples_dict = defaultdict(list)
ref_bases = []
for i, record in enumerate(vcf_reader.fetch("Y")):
    if i % 1000000 == 0: print(f"\r{i} positions processed", end="")
    if args.no_damage and (record.REF == "C" and record.ALT[0] == "T" or record.REF == "G" and record.ALT[0] == "A"):
        continue
    ref_bases.append(record.REF)
    for name in samples:
        call = record.genotype(name)
        called_base = call.gt_bases if call.gt_bases else "N"
        samples_dict[call.sample].append(called_base)

gt_df = pd.DataFrame(samples_dict)
ref_bases = pd.Series(ref_bases)

if args.variable:
    # count the number of alleles observed at each site and subset the GT
    # dataframe only to biallelic sites
    allele_counts = gt_df.apply(lambda row: len(set(i for i in row if i != "N")), axis=1)
    gt_df = gt_df.loc[allele_counts > 1]

    const_sites = dict(ref_bases[allele_counts == 1].value_counts())
    with open(re.sub(".fa$", ".counts", args.fasta), "w") as output:
        print(" ".join(str(const_sites[i]) for i in "ACGT"), file=output)

# write out the called bases for each sample in a FASTA format
with open(args.fasta, "w") as output:
    for name in samples:
        print(">" + re.sub("-", "_", name), file=output)
        print("".join(gt_df[name]), file=output)
