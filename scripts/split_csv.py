#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
splits input file (commonly csv) into n files equal size
leaving header and creates n files with the same
input name but adding .$i to the end
Author: Juan Manuel Crescente <crescente.juan@inta.gob.ar>
'''
import argparse
from math import floor
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input file (\n separated)",  required=True)
parser.add_argument("-n", "--nfiles",default=1, help="How much files to be splitted into (default 1)", type=int)
args = parser.parse_args()

#create output files
output_files = []
for i in range(0,args.nfiles):
	output_files.append(open(args.input + "." + str(i), 'wb'))

#open input csv
input_file_handler = open(args.input, 'r')
input_file = input_file_handler.read().split('\n')
#get header
header = input_file[0]
#get size of each file
#total lenght
count = len(input_file)
#each n-1 output file lenght
division = int(floor(count / args.nfiles))
#last file lenght, minus header
rest = count - (int(division) * args.nfiles)
curr_file = 0
curr_cursor = 0
for files in range(0, args.nfiles - 1):
	output_files[curr_file].write(header + "\n")
	for i in range(0, division):
		curr_cursor += 1
		output_files[curr_file].write(input_file[curr_cursor] + "\n")
	curr_file += 1

output_files[curr_file].write(header + "\n")
for i in range(0, division + rest - 2):
	curr_cursor += 1
	output_files[curr_file].write(input_file[curr_cursor] + "\n")

#close all files
for i in range(0,args.nfiles):
	output_files[i].close()
input_file_handler.close()
