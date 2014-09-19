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
        return D

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

    print(upperAlign[::-1])
    print(lowerAlign[::-1])

    return D[len(seq1), len(seq2)]

def exact_3(seq1, seq2, seq3, y, gapcost, backtrack = False): 
    D = numpy.zeros(shape=(len(seq1)+1, len(seq2)+1, len(seq3)+1))
    D12 = global_linear(seq1, seq2, y, gapcost)
    D13 = global_linear(seq1, seq3, y, gapcost)
    D23 = global_linear(seq2, seq3, y, gapcost)
    for i in range(0, len(seq1)):
        for j in range(0, len(seq2)): 
            D[i, j, 0] = D12[i, j] + (i + j) * gapcost
            D[i, 0, j] = D13[i, j] + (i + j) * gapcost
            D[0, i, j] = D23[i, j] + (i + j) * gapcost

    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1): 
            for k in range(1, len(seq3)+1):
                cij = y(seq1[i-1], seq2[j-1])
                cik = y(seq1[i-1], seq3[k-1])
                cjk = y(seq2[j-1], seq3[k-1])

                d1 = D[i-1, j-1, k-1] + cij + cik + cjk
                d2 = D[i-1, j-1, k] + cij + 2 * gapcost
                d3 = D[i-1, j, k-1] + cik + 2 * gapcost
                d4 = D[i, j-1, k-1] + cjk + 2 * gapcost
                d5 = D[i-1, j, k] + 2 * gapcost
                d6 = D[i, j-1, k] + 2 * gapcost
                d7 = D[i, j, k-1] + 2 * gapcost

                D[i, j, k] = min(d1, d2, d3, d4, d5, d6, d7)
    if not backtrack: 
        return D[len(seq1), len(seq2), len(seq3)]

    upperAlign = ''
    middleAlign = ''
    lowerAlign = ''

    i, j, k = len(seq1), len(seq2), len(seq3)
    while i != 0 or j != 0 or k != 0:
        a = seq1[i-1] if i > 0 else 0
        b = seq2[j-1] if j > 0 else 0
        c = seq3[k-1] if k > 0 else 0

        cij = y(seq1[i-1], seq2[j-1]) if i > 0 and j > 0 else 0
        cik = y(seq1[i-1], seq3[k-1]) if i > 0 and k > 0 else 0 
        cjk = y(seq2[j-1], seq3[k-1]) if j > 0 and k > 0 else 0

        d5 = D[i-1, j, k] + 2 * gapcost
        d6 = D[i, j-1, k] + 2 * gapcost
        d7 = D[i, j, k-1] + 2 * gapcost


        if j>0 and i>0 and k>0 and D[i, j, k] == D[i-1, j-1, k-1] + cij + cik + cjk:
            upperAlign = upperAlign + a
            middleAlign = middleAlign + b
            lowerAlign = lowerAlign + c
            i, j, k = i-1, j-1, k-1
        elif i > 0 and j > 0 and D[i, j, k] == D[i-1, j-1, k] + cij + 2 * gapcost:
            upperAlign = upperAlign + a
            middleAlign = middleAlign + b
            lowerAlign = lowerAlign + '-'
            i, j = i-1, j-1
        elif i > 0 and k > 0 and D[i, j, k] == D[i-1, j, k-1] + cik + 2 * gapcost:
            upperAlign = upperAlign + a
            middleAlign = middleAlign + '-'
            lowerAlign = lowerAlign + c
            i, k = i-1, k-1
        elif j > 0 and k > 0 and D[i, j, k] == D[i, j-1, k-1] + cjk + 2 * gapcost:
            upperAlign = upperAlign + '-'
            middleAlign = middleAlign + b
            lowerAlign = lowerAlign + c
            j, k = j-1, k-1
        elif i > 0 and D[i, j, k] == D[i-1, j, k] + 2 * gapcost:
            upperAlign = upperAlign + a
            middleAlign = middleAlign + '-' 
            lowerAlign = lowerAlign + '-'
            i = i-1
        elif j > 0 and D[i, j, k] == D[i, j-1, k] + 2 * gapcost:
            upperAlign = upperAlign + '-'
            middleAlign = middleAlign + b 
            lowerAlign = lowerAlign + '-'
            j = j-1
        elif k > 0 and D[i, j, k] == D[i, j, k-1] + 2 * gapcost:
            upperAlign = upperAlign + '-'
            middleAlign = middleAlign + '-' 
            lowerAlign = lowerAlign + c
            k = k-1
        else: 
            print("LOL SYSTEM EXIT ONE!")
            sys.exit(1)
    print(upperAlign[::-1])
    print(middleAlign[::-1])
    print(lowerAlign[::-1])

    return D[len(seq1), len(seq2), len(seq3)]        
