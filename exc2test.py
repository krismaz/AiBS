#! /usr/bin/python3
import sys
import numpy
import AiBS
from Bio import SeqIO


fasta_sequences = SeqIO.to_dict(SeqIO.parse(open('exc1.fasta'),'fasta'))




seq1 = fasta_sequences['Seq1'].seq.tostring()
seq2 = fasta_sequences['Seq2'].seq.tostring()

print(AiBS.global_linear(seq1, seq2, AiBS.complex_y, -5))
