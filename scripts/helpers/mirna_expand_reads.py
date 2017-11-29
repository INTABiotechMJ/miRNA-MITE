#!/usr/bin/env python
# -*- coding: utf-8 -*-
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def mirna_expand(libraries, output_fasta):
    """
    """
    buffer_seqs = []
    count = 0
    for library in libraries:
        lib_file = open(library)
        for line in lib_file:
            if len(line.split('\t')) < 2:
                continue
            sequence = line.split('\t')[0].strip()
            number = line.split('\t')[1]
            if sequence:
                for i in range(int(number)):
                    count += 1
                    record = SeqRecord(Seq(sequence), id=str(count), description=" ")
                    buffer_seqs.append(record)
    SeqIO.write(buffer_seqs, output_fasta, "fasta")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()#pylint: disable=invalid-name
    parser.add_argument("-l", "--libraries", help="(.fasta)", required=True, nargs='+')
    parser.add_argument("-o", "--output", help="Output file name (.fasta format)", required=True)
    args = parser.parse_args()#pylint: disable=invalid-name
    mirna_expand(args.libraries, args.output)