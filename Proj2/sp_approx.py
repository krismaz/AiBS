#!/usr/bin/python3
import sys
import numpy
import AiBS
from Bio import SeqIO
import argparse
from itertools import permutations

parser = argparse.ArgumentParser(description='Genome Sequence Analyzor')
parser.add_argument('cost', type=str, help='Cost file to read')
parser.add_argument('fasta', type=str, help='Fasta file to read')
parser.add_argument('--permutations', dest='permutations', action='store_true')
parser.set_defaults(permutations=False)
parser.add_argument('seqs', type=str, nargs='+', help='Sequences')

args = parser.parse_args()

fasta_sequences = SeqIO.to_dict(SeqIO.parse(open(args.fasta), 'fasta'))

seqs = []
for seq in args.seqs: 
    seqs.append(str(fasta_sequences[seq].seq).lower())

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

if args.permutations:
    AiBS.multi_approx(input_y, 5, *seqs)
else: 
    alignment = AiBS.approx(input_y, 5, *seqs)
    for i in range(len(alignment)):
        print(">seq" + str(i+1))
        print(alignment[i])