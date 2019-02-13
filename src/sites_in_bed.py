#!/usr/bin/env python3

import sys
import argparse
from pybedtools import BedTool

parser = argparse.ArgumentParser(
    description="Generate a list of all positions in given BED file")
parser.add_argument("--bed", help="Input file in a BED format", required=True)
parser.add_argument("--output", help="Output file (stdout if not specified)")
parser.add_argument("--format", help="Output format (BED - chrom, start, end; POS - chrom, pos)",
                    choices=["bed", "pos"], required=True)
args = parser.parse_args()

bed_regions = BedTool(args.bed)

# get a list of pairs of (chrom, position) for each position
# in regions in a given BED file
positions = ((r.chrom, pos) for r in bed_regions
                            for pos in range(r.start + 1, r.end + 1))

output = open(args.output, "w") if args.output else sys.stdout

for chrom, pos in positions:
    print(chrom, end="\t", file=output)

    if args.format == "BED":
        print(pos - 1, end="\t", file=output)

    print(pos, file=output)

if output != sys.stdout:
    output.close()
