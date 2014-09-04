#! /usr/bin/python3
import sys
import numpy
from Bio import SeqIO

complex_ycost = numpy.matrix('10 2 5 2; 2 10 2 5; 5 2 10 2; 2 5 2 10')

letters = ['A', 'C', 'G', 'T']

def complex_y(a, b):
	return complex_ycost[letters.index(a), letters.index(b)]


def global_linear(seq1, seq2, y, gapcost, bactrack = False):
	D = numpy.zeros(shape=(len(seq1)+1,len(seq2)+1))
	Dopt = numpy.zeros(shape=(len(seq1)+1,len(seq2)+1))

	for i in range(1,len(seq1)+1):
		D[i,0] = D[i-1,0] + gapcost

	for j in range(1,len(seq2)+1):
		D[0,j] = D[0,j-1] + gapcost

	for i in range(1,len(seq1)+1):
		for j in range(1,len(seq2)+1):
			a = seq1[i-1]
			b = seq2[j-1]
			v1 = D[i-1, j-1] + y(a, b)
			v2 = D[i, j-1] + gapcost
			v3 = D[i-1, j] + gapcost
			D[i, j] = max(v1, v2, v3)

	if not bactrack:
		return D[len(seq1), len(seq2)]

	upperAlign = ''
	lowerAlign = ''

	i, j = len(seq1), len(seq2)
	while i != 0 or j != 0:
		a = seq1[i-1]
		b = seq2[j-1]
		if j>0 and i>0 and D[i, j] == D[i-1, j-1] + y(a, b):
			upperAlign = upperAlign + a
			lowerAlign = lowerAlign + b
			i, j = i-1, j-1
		elif j>0 and D[i, j] == D[i, j-1] + gapcost:
			upperAlign = upperAlign + '-'
			lowerAlign = lowerAlign + b
			j = j-1
		elif i>0 and D[i, j] == D[i-1, j] + gapcost:
			upperAlign = upperAlign + a
			lowerAlign = lowerAlign + '-'
			i = i-1

	print(upperAlign)
	print(lowerAlign)

	return D[len(seq1), len(seq2)]
