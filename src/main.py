import numpy as np
import itertools

# define main functions

GAP_PENALTY = -4


def read_fasta(fasta_file):
    seq_id_counter = 0
    sequence_dict = {}
    with open(fasta_file, "r") as filin:
        seq_id = ""
        for line in filin:
            if line.startswith(">"):
                seq_id_counter += 1
                seq_id = seq_id_counter
                sequence_dict[seq_id] = ""
            else:
                sequence_dict[seq_id] += line.strip()
    return sequence_dict


def read_blosum(matrix_file):
    # initialize the dictionary
    blosum_dict = {}
    with open(matrix_file, "r") as filein:
        lines = filein.readlines()

    # Creates a list of AA from the first line of the file
    amino_acids = lines[0].split()[0:]

    # Splits the matrix starting second line
    for line in lines[1:]:
        content = line.split()
        amino_acid1 = content[0]  # row names are AA

        for i, amino_acid2 in enumerate(amino_acids):  # rows, columns
            key = (amino_acid1, amino_acid2)
            value = int(content[i + 1])
            blosum_dict[key] = value

    return blosum_dict


def pairwise_alignment(seq_m, seq_n, matrix, algt_score=False):
    # m rows, n columns as per convention
    m = len(seq_m)
    n = len(seq_n)

    # generates a matrix of size n+1 m+1 to store gap penalty
    scores = np.zeros((m+1, n+1))

    # fill the first row and column with the gap penalty
    scores[:, 0] = [GAP_PENALTY * i for i in range(scores.shape[0])]
    scores[0, :] = [GAP_PENALTY * i for i in range(scores.shape[1])]

    # fill the score matrix
    for i in range(1, scores.shape[0]):
        for j in range(1, scores.shape[1]):
            match = scores[i-1, j-1] + matrix[(seq_m[i-1], seq_n[j-1])]
            gap1 = scores[i-1, j] + GAP_PENALTY
            gap2 = scores[i, j-1] + GAP_PENALTY
            scores[i, j] = max(match, gap1, gap2)

    # TODO: kmers ?

    # Traceback to find the optimal alignement
    # TODO : the traceback lol


    if algt_score == True:
        return scores
    else:
        return scores[m, n]


def calculate_score(sequences, matrix):
    seq_ids = list(sequences.keys())
    size = len(seq_ids)

    # initialize the matrix with zeros
    scores_matrix = np.zeros((size, size))

    # iterate over the possible indexes
    for i in range (size):
        seq1 = sequences[i+1]
        # j iterates between i+1 and size of matrix so that we only calculate the diagonal (limits compute time)
        for j in range(i+1, size):
            seq2 = sequences[j+1]
            alt_score = pairwise_alignment(seq1, seq2, matrix)
            scores_matrix[i, j] = alt_score
    print(scores_matrix)
    return scores_matrix


# def create_guide_tree(sequences, matrix):
# creates guide tree

# def run_multiple_alignment(sequences, tree)
# uses tree to align the sequences

blosum62 = read_blosum("blossum_62.txt")

my_seqs = read_fasta("test.fasta")
matrice = calculate_score(my_seqs, blosum62)

