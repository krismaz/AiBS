#! /usr/bin/python3
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
parser.add_argument('--backtrack', dest='backtrack', action='store_true')
parser.set_defaults(backtrack=False)
parser.add_argument('--affine', dest='affine', action='store_true')
parser.set_defaults(affine=False)
parser.add_argument('--matrix', dest='matrix', action='store_true')
parser.set_defaults(matrix=False)

args = parser.parse_args()

fasta_sequences = SeqIO.to_dict(SeqIO.parse(open(args.fasta),'fasta'))

seq1 = str(fasta_sequences[args.seq1].seq)
seq2 = str(fasta_sequences[args.seq2].seq)

with open(args.cost, 'r') as cost_file:
	letters = cost_file.readline().lower()
	complex_ycost = numpy.matrix(cost_file.read())




def input_y(a, b):
	try:
		return complex_ycost[letters.index(a), letters.index(b)]
	except:
		print('WTF ' + a + ' - ' + b)
		return 0 


if args.matrix:
	keys = list(fasta_sequences.keys())
	keys.sort()
	res = numpy.zeros(shape=(len(keys),len(keys)))
	print(keys)
	for s1, s2 in list(permutations(keys,2)):
		if args.affine:
			res[keys.index(s1), keys.index(s2)] = AiBS.global_affine(fasta_sequences[s1], fasta_sequences[s2],  input_y, 5, 5, False)
		else: 
			res[keys.index(s1), keys.index(s2)] = AiBS.global_linear(fasta_sequences[s1], fasta_sequences[s2],  input_y, 5, False)
	print(res)
	sys.exit()

if args.affine:
	print(AiBS.global_affine(seq1, seq2,  input_y, 5, 5, args.backtrack))
else: 
	print(AiBS.global_linear(seq1, seq2,  input_y, 5, args.backtrack))
