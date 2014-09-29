#!/usr/bin/python3
import sys
import numpy
import AiBS
from Bio import SeqIO
import argparse
from itertools import permutations

parser = argparse.ArgumentParser(description='Genome Sequence Analyzor')
parser.add_argument('cost', type=str, help='Cost file to read')
parser.add_argument('--backtrack', dest='backtrack', action='store_true')
parser.add_argument('fastas', type=str, nargs='+', help='Fasta')

args = parser.parse_args()

def input_y(a, b):
	try:
		return complex_ycost[letters.index(a), letters.index(b)]
	except:
		print('WTF ' + a + ' - ' + b + ' is not a real cost')
		sys.exit(1)
		return 0

print("#len\texact\tapprox")

for fasta in args.fastas:

	fasta_sequences = SeqIO.to_dict(SeqIO.parse(open(fasta),'fasta'))

	seq1 = str(list(fasta_sequences.values())[0].seq).lower()
	seq2 = str(list(fasta_sequences.values())[1].seq).lower()
	seq3 = str(list(fasta_sequences.values())[2].seq).lower()

	with open(args.cost, 'r') as cost_file:
		letters = cost_file.readline().lower()
		complex_ycost = numpy.matrix(cost_file.read())

	print(str(len(seq1)) + "\t" + str(AiBS.exact_3(seq1, seq2, seq3, input_y, 5, args.backtrack)) + "\t" + str(AiBS.approx_score(input_y, 5, seq1, seq2, seq3)))
