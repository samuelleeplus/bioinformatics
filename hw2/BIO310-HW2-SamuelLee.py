from pylab import *
from Bio import SeqIO
from Bio.Seq import Seq
import linecache


class ScoreParam:
    def __init__(self, match, mismatch, gap, gap_start=0):
        self.gap_start = gap_start
        self.gap = gap
        self.match = match
        self.mismatch = mismatch

    def matchchar(self, a, b):
        assert len(a) == len(b) == 1
        if a == b:
            return self.match
        else:
            return self.mismatch

    def __str__(self):
        return "match = %d; mismatch = %d; gap_extend = %d" % (
            self.match, self.mismatch, self.gap)


def make_matrix(sizex, sizey):
    return [[0] * sizey for i in range(sizex)]


def print_matrix(x, y, A):
    # decide whether there is a 0th row/column
    if len(x) == len(A):
        print("%5s" % (" "), end=' ')
    else:
        print("%5s %5s" % (" ", "*"), end=' ')
        y = "-" + y

    # print the top row
    for c in x:
        print("%5s" % (c), end=' ')
    print()

    for j in range(len(A[0])):
        print("%5s" % (y[j]), end=' ')
        for i in range(len(A)):
            print("%5.0f" % (A[i][j]), end=' ')
        print()


def local_align(x, y, score=ScoreParam(gapScore, matchScore, mismatchScore)):
    A = make_matrix(len(x) + 1, len(y) + 1)
    A[0][0] = 0
    rows = len(rowSequence)
    cols = len(columnSequence)
    for i in range(1, rows):
        A[i][0] = A[i - 1][0] + gapScore

    for j in range(1, cols):
        A[0][j] = A[0][j - 1] + gapScore

    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            A[i][j] = max(A[i][j - 1] + gapScore, A[i - 1][j] + gapScore,
                          A[i - 1][j - 1] + score.matchchar(x[i - 1], y[j - 1]), 0)
            # track the cell with the largest score

    print("Scoring:", str(score))
    print("A matrix =")
    print_matrix(x, y, A)
    # return the opt score and the best location


# Import pairwise2 module
from Bio import pairwise2
# Import format_alignment method
from Bio.pairwise2 import format_alignment
import sys

fileName = input("Please enter the name of the file: ")
mismatchScore = int(input("Please enter the mismatch score: "))
gapScore = int(input("Please enter the gap score: "))
matchScore = int(input("Please enter the match score: "))
# open file in read mode
file = open(fileName, "r")
# store the sequences into a file, given it will be separated by a line
rowSequence = linecache.getline(fileName, 1)
columnSequence = linecache.getline(fileName, 2)

sys.stdout = open('output.txt', 'wt')

alignments = pairwise2.align.globalms(rowSequence[:-1], columnSequence[:-1], matchScore, mismatchScore, gapScore,
                                      gapScore)

# Use format_alignment method to format the alignments in the list
for a in alignments:
    print(format_alignment(*a))

local_align(rowSequence, columnSequence, ScoreParam(gap=gapScore, match=matchScore, mismatch=mismatchScore))