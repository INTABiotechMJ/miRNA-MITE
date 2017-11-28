#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, csv, os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import argparse
import pandas as pd
#this script creates miRNA clusters
#separated for maximum this variable
#then creates a fasta containing only the section of this clusters

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--allhits", help="Input file (blastn file)",  required=True)
parser.add_argument("-csz", "--cluster_sep_size", help="Number of max bp to separate each cluster",  required=True, type=int)
parser.add_argument("-o", "--output_cluster", help="Output file (.csv format)",  required=True)
args = parser.parse_args()
#if 2 miRNA are separated for less than cluster_sep_max_size, count on the same cluster
cluster_sep_size = args.cluster_sep_size

allhits = pd.read_csv(args.allhits, sep="\t",header=None, usecols=[0,1,8,9], names=['miRNA', 'scaffold', 'bpfrom', 'bpto'])
allhits['checked'] = 0
cluster = pd.DataFrame([], columns=['scaffold','bpfrom','bpto','miRNAs'])

#if some result is inverted, change from for to
#for index, row in allhits.iterrows():
#	if row['bpfrom'] > row['bpto']:
#		rfrom = row['bpfrom']
#		allhits.loc[index,'bpfrom'] = row['bpto']
#		allhits.loc[index,'bpto'] = rfrom
#		allhits.drop(index)

changed = True
for index, row in allhits.iterrows():
	newbpfrom = None
	newbpto = None
	if allhits.at[index,'checked']:
		continue
	scaffold = row['scaffold']
	bpfrom = row['bpfrom']
	miRNA = row['miRNA']
	bpto = row['bpto']
	if bpto < bpfrom:
		bpaux = bpfrom
		bpfrom = bpto
		bpto = bpaux
	query = allhits.query('checked == 0 and scaffold == @scaffold and ((bpto + @cluster_sep_size >= @bpfrom and bpfrom  - @cluster_sep_size <= @bpfrom) or (bpto + @cluster_sep_size >= @bpto and bpfrom  - @cluster_sep_size <= @bpto))')
	allhits.loc[index,'checked'] = 1
	#if query.empty:
	#	print allhits[allhits.checked==0].shape
	#	cluster = cluster.append({'scaffold': scaffold, 'bpto': bpto, 'bpfrom':bpfrom,'miRNAs':row['miRNA']}, ignore_index=True)
	miRNAs = []
	for indexq, rowq in query.iterrows():
		miRNAs.append(rowq['miRNA'])
		if newbpto == None:
			newbpto = rowq['bpto']
		if newbpfrom == None:
			newbpfrom = rowq['bpfrom']
		if rowq['bpto'] < rowq['bpfrom']:
			aux_to = rowq['bpfrom']
			aux_from = rowq['bpto']
		else:
			aux_from = rowq['bpfrom']
			aux_to = rowq['bpto']
			#print rowq['miRNA'],newbpfrom, newbpto
		newbpfrom = min(newbpfrom,aux_from)
		newbpto = max(newbpto,aux_to)
		allhits.loc[indexq,'checked'] = 1

		#allhits.loc[rowq['scaffold'],'bpfrom'] = min(row['bpfrom'],bpfrom_2)
		#allhits.loc[rowq['scaffold'],'bpto'] = max(row['bpto'],bpto_2)
		#allhits.drop(index, axis=0, inplace=True)
		#allhits.index = range(len(allhits))
	if not query.empty:
		cluster = cluster.append({'scaffold': scaffold, 'bpto': newbpto, 'bpfrom':newbpfrom,'miRNAs':";".join(miRNAs)}, ignore_index=True)
cluster.to_csv(args.output_cluster,index=False)
