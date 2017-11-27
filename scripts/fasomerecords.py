#!/usr/bin/env python
# -*- coding: utf-8 -*-
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def exclude(sequence, records, output_fasta):
    fasta_seq = SeqIO.parse(sequence,'fasta')
    buffer = []
    for record in fasta_seq:
        if record.id in records:
            continue
        buffer.append(record)
    SeqIO.write(buffer, output_fasta, "fasta")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sequence", help="Sequence file (.fasta)",  required=True)
    parser.add_argument("-e", "--exclude", help="Exclude some ids", action='append')
    parser.add_argument("-o", "--output", help="Output file name (.fasta format)",  required=True)
    args = parser.parse_args()
    exclude(args.sequence, args.exclude, args.output)