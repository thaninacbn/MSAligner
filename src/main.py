import numpy as np

# define main functions

GAP_PENALTY = -4


def read_fasta(fasta_file):
    # code from previous python exos
    with open(fasta_file, "r") as filin:
        sequence_list = []
        for line in filin:
            if line.startswith(">"):
                pass
            else:
                read = filin.readline()
                read = read.strip()
                sequence_list.append(read)
                sequence = "".join(sequence_list)
    # TODO: deal with the fact that there might be more than 1 sequence per file
    return sequence


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


def pairwise_alignment(seq1, seq2, matrix):
    # store the length of the sequences
    n = len(seq1)
    m = len(seq2)

    # generates a matrix of size n m
    scores = np.zeros(m + 1, n + 1)

    # fill the first row and column with the gap penalty
    for j in range(0, m + 1):
        scores[0][j] = GAP_PENALTY * j

    for i in range(0, n + 1):
        scores[0][i] = GAP_PENALTY * i

    # fill the score matrix
    for i in range (1, m+1):
        for j in range (1, n+1):
            alignscore = scores[i-1][j-1]+matrix[()]



# def calculate_score (sequences):
# creates scores matrix


# def create_guide_tree(sequences, matrix):
# creates guide tree

# def run_multiple_alignment(sequences, tree)
# uses tree to align the sequences

blosum62 = read_blosum("blossum_62.txt")
print(blosum62)
