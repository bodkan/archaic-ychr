#!/usr/bin/env python3

import argparse
import re
import os.path

parser = argparse.ArgumentParser()

parser.add_argument("--xml", help="Path to the original XML file", required=True)
parser.add_argument("--counts", help="Path to the file with invariant base counts", required=True)

args = parser.parse_args()

counts = open(args.counts, "r").readline().strip()
xml = open(args.xml, "r").read()

alignment_id = os.path.splitext(os.path.basename(args.counts))[0]

xml = re.sub(rf'<data\nid="({alignment_id})"', r'<data\nid="\g<1>_original"', xml)

xml = re.sub(r'</data>', rf'</data>\n<data id="{alignment_id}" spec="FilteredAlignment" filter="-" data="@{alignment_id}_original" constantSiteWeights="{counts}"/>', xml)

edited_xml = open(re.sub(r'.xml', r'_edited.xml', args.xml), "w")
print(xml, file=edited_xml)
edited_xml.close()
