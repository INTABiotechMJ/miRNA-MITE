#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def fasta_from_ann(annotation, sequence, feature, windows, output_fasta):
    df_gff = pd.read_csv(annotation, index_col=False, sep='\t',header=None)
    df_gff.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    fasta_seq = SeqIO.parse(sequence,'fasta')
    buffer = []
    for record in fasta_seq:
        df_exctract = df_gff[(df_gff.seqname == record.id) & (df_gff.feature == feature)]
        for k,v in df_exctract.iterrows():
            clean_seq = ''.join(str(record.seq).splitlines())
            if int(v.start) - windows < 0:
                start = 0
            else:
                start = int(v.start) - windows
            if int(v.end) + windows > len(clean_seq):
                end = len(clean_seq)
            else:
                end = int(v.end) + windows
            new_seq = clean_seq[start:end]
            new_id = record.id + "_from_" + str(v.start) + "_to_" + str(v.end) + "_feature_" + v.feature
            desc = "attribute: " + v.attribute + " strand: " + v.strand
            seq = SeqRecord(Seq(new_seq), id=new_id,description = desc)
            buffer.append(seq)
    SeqIO.write(buffer, output_fasta, "fasta")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--annotation", help="Annotation file (.gff3 format)",  required=True)
    parser.add_argument("-s", "--sequence", help="Sequence file (.fasta)",  required=True)
    parser.add_argument("-f", "--feature", help="Feature to extract (ie. gene)",  required=True)
    parser.add_argument("-w", "--windows", help="+- nt to cut from and to",  required=True)
    parser.add_argument("-o", "--output", help="Output file name (.fasta format)",  required=True)
    args = parser.parse_args()
    fasta_from_ann(args.annotation, args.sequence, args.feature, args.windows, args.output)