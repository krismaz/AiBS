! /usr/bin/python3
import sys
import numpy
import AiBS
from Bio import SeqIO
import argparse
from itertools import permutations

parser = argparse.ArgumentParser(description='Genome Sequence Analyzor')
parser.add_argument('cost', type=str, help='Cost file to read')
parser.add_argument('fasta', type=str, help='Fasta file to read')
parser.add_argument('seq1', type=str, help='Sequence 1')
parser.add_argument('seq2', type=str, help='Sequence 2')
parser.add_argument('seq3', type=str, help='Sequence 3')

args = parser.parse_args()

fasta_sequences = SeqIO.to_dict(SeqIO.parse(open(args.fasta),'fasta'))

seq1 = str(fasta_sequences[args.seq1].seq)
seq2 = str(fasta_sequences[args.seq2].seq)
seq3 = str(fasta_sequences[args.seq2].seq)

with open(args.cost, 'r') as cost_file:
	letters = cost_file.readline().lower()
	complex_ycost = numpy.matrix(cost_file.read())




def input_y(a, b):
	try:
		return complex_ycost[letters.index(a), letters.index(b)]
	except:
		print('WTF ' + a + ' - ' + b + ' is not a real cost')
		sys.exit(1)
		return 0 
