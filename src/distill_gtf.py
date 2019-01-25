#!/usr/bin/python

# PSEUDOCODE:
# already grep for the ENSG of interest before feeding to stdin
# where the 3rd col is "gene" (probably 1st row), grab the gene_name field and make that the output filename
# only care about rows where 3rd col is "CDS"
# extract the transcript_id value
# write out a table of gene_id, transcript_id, cds exon number, start, stop

# then, tasks to do (in Python or R?)
# collapse all transcripts into one, with different spans labeled constitutive or non-constitutive. did I already do this back in the day?
# map actual chromosomal position onto an x axis that trims out all but +/- 10 bp of each intron
# a start on this is in src/dbLOF/src/dataviz.R ~line 186 but it doesn't leave any intron length at all for plotting
# and the collapsing concept has not yet been implemented

import sys
import re

for line in sys.stdin:
	cells = line.strip().split('\t')
	if (cells[1] == 'ensembl' or cells[1] == 'ensembl_havana') and cells[2] == 'CDS':
		uglylist = cells[8].strip().split(';')
		transcript_id_field = uglylist[2]
		transcript_id = re.findall('"([^"]*)"',transcript_id_field)[0]
		chrom = cells[0]
		start = cells[3]
		stop = cells[4]
		print transcript_id, chrom, start, stop


# gunzip -c raw_data/Homo_sapiens.GRCh37.87.gtf.gz | grep ENSG00000186868 | python src/distill_gtf.py > raw_data/MAPT.distilled.gtf

