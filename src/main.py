# import statements
import numpy as np
import os
import sys

# parwise alignement
def pairwise_align(seq1, seq2):
    #compute pairwise alignment using needleman-wunsch algorithm

def run_alignments(input):
    #compute pairwise_align on all pairs of seqs
    #each thread computes an alignment

def create_similarity_matrix(alignments):
    #takes all the alignments and creates a similarity matrix

def neighbour_joining(matrix):
    #creates neighbour joining tree using the similarity matrix from create_similarity_matrix
    #NJ more optimal than UPGMA (Onm instead of On^3)

def multiple_sequence_alignment(sequences,tree):
    #runs the actual alignment, returns aligned sequences