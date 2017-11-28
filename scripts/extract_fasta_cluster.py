#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, csv, os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import argparse
import pandas as pd
buffer_size = 10
#this script creates miRNA clusters
#separated for maximum this variable
#then creates a fasta containing only the section of this clusters

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cluster", help="Input file (csv)",  required=True)
parser.add_argument("-f", "--fasta", help="Input file (.fasta format)",  required=True)
parser.add_argument("-wz", "--windows_size", help="Number of bp to split to left and right",  required=True, type=int)
args = parser.parse_args()
#in clusters leave windows_size for each side
windows_size = args.windows_size
allhits = []
output_file = open(args.cluster + ".fasta", 'wb')
allhit = pd.read_csv(args.cluster , sep=",")

count = 0
dumped = False
curr_file = 0
buffer = []
for record in SeqIO.parse(args.fasta,'fasta'):
	clean_seq = ''.join(str(record.seq).splitlines())
	clean_seq_len = len(clean_seq)
	index_file = 0
	for index, row in allhit[allhit['scaffold'] == record.id].iterrows():
		start = int(row['bpfrom'])
		end = int(row['bpto'])
		if start > end:
			aux_end = end
			end = start
			start = aux_end
		if start - windows_size <= 0:
			nstart = 0
		else:
			nstart = start - windows_size
		if end + windows_size >= clean_seq_len:
			nend = clean_seq_len
		else:
			nend = end + windows_size
		new_seq = clean_seq[int(nstart):int(nend)]
		new_id = record.id + "_" + str(start) + "_" + str(end)
		seq = SeqRecord(Seq(new_seq), id=new_id ,description = "new_from:" + str(nstart) + " new_to:" + str(nend) + " cluster_index: " + str(index))
		buffer.append(seq)
SeqIO.write(buffer, output_file, "fasta")
