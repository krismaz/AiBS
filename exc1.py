#! /usr/bin/python3
import sys
import numpy
from Bio import SeqIO


fasta_sequences = SeqIO.to_dict(SeqIO.parse(open('exc1.fasta'),'fasta'))



gapcost = -5

cost = numpy.matrix('10 2 5 2; 2 10 2 5; 5 2 10 2; 2 5 2 10')

letters = ['A', 'C', 'G', 'T']

def y(a, b):
	return cost[letters.index(a), letters.index(b)]

#seq1 = str(sequences['Seq1']).replace(' ','')
#seq2 = str(sequences['Seq2']).replace(' ','')

#seq1 = 'AATAAT'
#seq2 = 'AAGG'

seq1 = fasta_sequences['Seq1'].seq.tostring()
seq2 = fasta_sequences['Seq2'].seq.tostring()


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

Dopt[0,0] = 1

for i in range(1,len(seq1)+1):
	Dopt[i,0] = 1

for j in range(1,len(seq2)+1):
	Dopt[0,j] = 1

for i in range(1,len(seq1)+1):
	for j in range(1,len(seq2)+1):
		count = 0
		a = seq1[i-1]
		b = seq2[j-1]
		if D[i,j] == D[i-1, j-1] + y(a, b):
			count += Dopt[i-1, j-1]
		if D[i, j] == D[i-1, j] + gapcost:
			count += Dopt[i-1, j]
		if D[i, j] == D[i, j-1] + gapcost:
			count += Dopt[i, j-1]
		Dopt[i, j] = count


upperAlign = upperAlign[::-1]
lowerAlign = lowerAlign[::-1]

print(seq1)

print(D)
print("Cost " + str(D[len(seq1), len(seq2)]))
print(upperAlign)
print(lowerAlign)
print("nalign " + str(Dopt[len(seq1), len(seq2)]))