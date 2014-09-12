#! /usr/bin/python3 
import sys
import numpy
from Bio import SeqIO

complex_ycost = numpy.matrix('0 5 2 5; 5 0 5 2; 2 5 0 5; 5 2 5 0')

letters = 'ACGT'.lower()

def complex_y(a, b):
	return complex_ycost[letters.index(a), letters.index(b)]

undef = float("inf")


def global_linear(seq1, seq2, y, gapcost, bactrack = False):
	seq1 = seq1.lower()
	seq2 = seq2.lower()
	D = numpy.zeros(shape=(len(seq1)+1,len(seq2)+1))

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
			D[i, j] = min(v1, v2, v3)

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

	print('>seq1')
	print(upperAlign[::-1])
	print('>seq2')
	print(lowerAlign[::-1])

	return D[len(seq1), len(seq2)]


def global_affine(seq1, seq2, y, a, b, bactrack = False):
	seq1 = seq1.lower()
	seq2 = seq2.lower()
	D = numpy.zeros(shape=(len(seq1)+1,len(seq2)+1))
	S = numpy.zeros(shape=(len(seq1)+1,len(seq2)+1))
	I = numpy.zeros(shape=(len(seq1)+1,len(seq2)+1))

 
	for i in range(0,len(seq1)+1):
		for j in range(0,len(seq2)+1):

			dv1, dv2 = undef, undef
			if i>0 and j>=0:
				dv1 = S[i-1, j] + (a+b)
			if i>1 and j>=0:
				dv2 = D[i-1, j] + a
			D[i,j] = min(dv1, dv2)

			iv1, iv2 = undef, undef
			if i>=0 and j>0:
				iv1 = S[i, j-1] + (a+b)
			if i>=0 and j>1:
				iv2 = I[i, j-1] + a
			I[i,j] = min(iv1, iv2)

			sv1, sv2, sv3, sv4 = undef, undef, undef, undef

			if i==0 and j==0:
				sv1 = 0
			if i>0 and j>0:
				sv2 = S[i-1, j-1] + y(seq1[i-1], seq2[j-1])
			if i>0 and j>=0:
				sv3 = D[i,j]
			if i>=0 and j>0:
				sv4 = I[i,j]

			S[i,j] = min(sv1, sv2, sv3, sv4)

	if not bactrack:
		return S[len(seq1), len(seq2)]

	i, j = len(seq1), len(seq2)
	upperAlign = ''
	lowerAlign = ''
	while i > 0 or j > 0:
		if (i>0 and j>0) and S[i,j] == S[i-1, j-1] + y(seq1[i-1], seq2[j-1]):
			upperAlign += seq1[i-1]
			lowerAlign += seq2[j-1]
			i -= 1
			j -= 1
		else:
			k = 1
			while True:
				if i>=k and S[i,j] == S[i-k, j] + (a*k + b):
					upperAlign += seq1[i-1:(i-k)-1:-1]
					lowerAlign += '-'*k
					i -= k
					break
				elif j>=k and S[i,j] == S[i, j-k] + (a*k + b):
					upperAlign += '-'*k
					lowerAlign += seq2[j-1:(j-k)-1:-1]
					j -= k
					break
				else:
					k += 1

	print('>seq1')
	print(upperAlign[::-1])
	print('>seq2')
	print(lowerAlign[::-1])

	return S[len(seq1), len(seq2)]



