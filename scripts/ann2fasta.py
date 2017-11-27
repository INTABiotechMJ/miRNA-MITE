#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def fasta_from_ann(annotation, sequence, windows, output_fasta):
    """Extract fasta files from annotations
    """
    df_gff = pd.read_csv(annotation, index_col=False, sep='\t', header=None)
    df_gff.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    fasta_seq = SeqIO.parse(sequence, 'fasta')
    buffer_seqs = []
    for record in fasta_seq:
        df_exctract = df_gff[(df_gff.seqname == record.id)]
        for key,val in df_exctract.iterrows():
            clean_seq = ''.join(str(record.seq).splitlines())
            if int(val.start) - windows < 0:
                start = 0
            else:
                start = int(val.start) - windows
            if int(val.end) + windows > len(clean_seq):
                end = len(clean_seq)
            else:
                end = int(val.end) + windows
            new_seq = clean_seq[start:end]
            new_id = record.id + "_from_" + str(val.start) + "_to_" + str(val.end) + "_feature_" + val.feature
            desc = "attribute: " + val.attribute + " strand: " + val.strand
            seq = SeqRecord(Seq(new_seq), id=new_id, description=desc)
            buffer_seqs.append(seq)
    SeqIO.write(buffer_seqs, output_fasta, "fasta")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()#pylint: disable=invalid-name
    parser.add_argument("-a", "--annotation", help="Annotation file (.gff3 format)", required=True)
    parser.add_argument("-s", "--sequence", help="Sequence file (.fasta)", required=True)
    parser.add_argument("-w", "--windows", help="+- nt to cut from and to", type=int, required=True)
    parser.add_argument("-o", "--output", help="Output file name (.fasta format)", required=True)
    args = parser.parse_args()#pylint: disable=invalid-name
    fasta_from_ann(args.annotation, args.sequence, args.windows, args.output)
