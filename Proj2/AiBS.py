#! /usr/bin/python3 
import sys
import numpy
from itertools import permutations
from Bio import SeqIO

complex_ycost = numpy.matrix('0 5 2 5; 5 0 5 2; 2 5 0 5; 5 2 5 0')

letters = 'ACGT'.lower()


def complex_y(a, b):
    return complex_ycost[letters.index(a), letters.index(b)]


undef = float("inf")


def global_linear(seq1, seq2, y, gapcost, backtrack = False, matrix = False):
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

    if not backtrack:
        return D if matrix else D[len(seq1), len(seq2)]

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

    return (upperAlign[::-1], lowerAlign[::-1])

def exact_3(seq1, seq2, seq3, y, gapcost, backtrack=False): 
    D = numpy.zeros(shape=(len(seq1)+1, len(seq2)+1, len(seq3)+1))
    D12 = global_linear(seq1, seq2, y, gapcost, matrix=True)
    D13 = global_linear(seq1, seq3, y, gapcost, matrix=True)
    D23 = global_linear(seq2, seq3, y, gapcost, matrix=True)

    for i in range(0, len(seq1)+1):
        for j in range(0, len(seq2)+1):
            D[i, j, 0] = D12[i, j] + (i + j) * gapcost
    for i in range(0, len(seq1)+1):
        for k in range(0, len(seq3)+1):
            D[i, 0, k] = D13[i, k] + (i + k) * gapcost
    for j in range(0, len(seq2)+1):
        for k in range(0, len(seq3)+1):
            D[0, j, k] = D23[j, k] + (j + k) * gapcost

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
            print("LOL SYSTEM EXIT ONE!",upperAlign, middleAlign, lowerAlign)
            sys.exit(1)
    print(upperAlign[::-1])
    print(middleAlign[::-1])
    print(lowerAlign[::-1])

    return D[len(seq1), len(seq2), len(seq3)]        


def approx(y, gapcost, *seqs):
    best_seq = None
    best_score = -1

    for seq1 in seqs: 
        score = sum([global_linear(seq1, seq2, y, gapcost) for seq2 in seqs if seq2 != seq1])
        if score > best_score:
            best_score = score
            best_seq = seq1

    columns = [ [c] for c in best_seq ]
    seqs = list(seqs)
    seqs.remove(best_seq)
    for seq in seqs: 
        A = global_linear(best_seq, seq, y, gapcost, backtrack=True)
        columns = extend(columns, *A)
    return transpose_columns(columns)

def approx_score(y, gapcost, *seqs):
    best_seq = None
    best_score = -1

    for seq1 in seqs: 
        score = sum([global_linear(seq1, seq2, y, gapcost) for seq2 in seqs if seq2 != seq1])
        if score > best_score:
            best_score = score
            best_seq = seq1

    columns = [ [c] for c in best_seq ]
    seqs = list(seqs)
    seqs.remove(best_seq)
    for seq in seqs: 
        A = global_linear(best_seq, seq, y, gapcost, backtrack=True)
        columns = extend(columns, *A)
    return compute_score(transpose_columns(columns))

def multi_approx(y, gapcost, *seqs):
    best_seq = None
    best_score = -1

    for seq1 in seqs: 
        score = sum([global_linear(seq1, seq2, y, gapcost) for seq2 in seqs if seq2 != seq1])
        if score > best_score:
            best_score = score
            best_seq = seq1

    seqs = list(seqs)
    print("index of best seq: " + str(seqs.index(best_seq)))
    seqs.remove(best_seq)
    for perm in permutations(list(range(len(seqs)))):
        columns = [ [c] for c in best_seq ]
        print("permutation: " + str(perm))
        for i in perm:
            A = global_linear(best_seq, seqs[i], y, gapcost, backtrack=True)
            columns = extend(columns, *A)
        alignment = list(transpose_columns(columns))
        print("Alignment: " + str(alignment))
        score = compute_score(alignment)
        print("Score: " + str(score))


def transpose_columns(columns):
    #return ["".join(x for y in columns for x in y)]
    return (''.join(col[j] for col in columns) for j,_ in enumerate(columns[0])) #Now Dan can be happy :D

def extend(columns, upperAlign, lowerAlign):
    c, i = 0, 0
    results = []
    col_len = len(columns[0])
    while(c < len(columns) or i < len(upperAlign)):
        if c < len(columns) and columns[c][0] == '-':
            columns[c].append('-')
            results.append(columns[c])
            c = c + 1
            continue
        if upperAlign[i] != '-' and lowerAlign[i] != '-':
            columns[c].append(lowerAlign[i])
            results.append(columns[c])
            c, i = c + 1, i + 1
            continue
        if upperAlign[i] == '-':
            temp = list('-' * col_len)
            temp.append(lowerAlign[i])
            results.append(temp)
            i = i + 1
            continue
        if lowerAlign[i] == '-':
            columns[c].append('-')
            results.append(columns[c])
            c, i = c + 1, i + 1
            continue
        print("SYSTEM FAILURE!!! PANIC")
        sys.exit(1)
    return results


def compute_score(data):
    cost = [[0, 5, 2, 5, 5],  # A
            [5, 0, 5, 2, 5],  # C
            [2, 5, 0, 5, 5],  # G
            [5, 2, 5, 0, 5],  # T
            [5, 5, 5, 5, 0]]  #-'
    # Compute the score of each induced pairwise alignment
    score = 0
    row = [str2seq(s) for s in list(data)]
    for i in range(len(row)):
        for j in range(i+1, len(row)):
            if len(row[i]) != len(row[j]):
                print("ERROR: Rows " + str(i) + " and " + str(j) + " have different lengths.")
                sys.exit(1)
            for c in range(len(row[i])):
                score = score + cost[row[i][c]][row[j][c]]
    return score

def str2seq(s):
    dict_str2seq = {'a':0, 'c':1, 'g':2, 't':3, 'A':0, 'C':1, 'G':2, 'T':3, '-':4, 'N':0, 'R':0, 'S':0}
    try:
        seq = [dict_str2seq[c] for c in s]
        return seq
    except KeyError as h:
        print("ERROR: Illegal character " + str(h) + " in input string.")
        sys.exit(1)
