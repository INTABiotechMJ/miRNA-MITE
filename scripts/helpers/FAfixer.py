#!/usr/bin/env python
# -*- coding: utf-8 -*-
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def FAfixer(input_fasta, output_fasta):
    """
    """
    buffer_seqs = []
    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        record = SeqRecord(seq_record.seq, id=seq_record.id, description=seq_record.description)
        buffer_seqs.append(record)
    SeqIO.write(buffer_seqs, output_fasta, "fasta")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()#pylint: disable=invalid-name
    parser.add_argument("-i", "--input_fasta", help="(.fasta)", required=True)
    parser.add_argument("-o", "--output_fasta", help="(.fasta)", required=True)
    args = parser.parse_args()#pylint: disable=invalid-name
    FAfixer(args.input_fasta, args.output_fasta)