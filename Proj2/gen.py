import argparse
import random

parser = argparse.ArgumentParser(description='Genome Sequence Generator')
parser.add_argument('chars', type=str, help='Fasta file to read')
parser.add_argument('N', type=int, help='Cost file to read')

args = parser.parse_args()

print('>seq1')
print(''.join([random.choice(args.chars) for i in range(args.N)]))
print('>seq2')
print(''.join([random.choice(args.chars) for i in range(args.N)]))
print('>seq3')
print(''.join([random.choice(args.chars) for i in range(args.N)]))
