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
    n = len(seq1)+1
    m = len(seq2)+1

    # generates a matrix of size n m
    scores = np.zeros((m, n))

    # fill the first row and column with the gap penalty
    for i in range(0, n):
        scores[0][i] = GAP_PENALTY * i

    for j in range(0, m):
        scores[0][j] = GAP_PENALTY * j

    # fill the score matrix
    for i in range(1, n):
        for j in range(1, m):
            match = scores[i-1][j-1]+matrix[(seq1[i-1], seq2[j-1])]
            gap1 = scores[i-1][j] + GAP_PENALTY
            gap2 = scores[i][j-1] + GAP_PENALTY # index ot of bounds here??
            scores[i][j] = max(match, gap1, gap2)

    # find the maximum alignment score in the matrix
    list_scores = scores[m, 0:n]+scores[0:m, n]
    alignment_score = max(list_scores)
    # TODO: kmers ?

    return alignment_score


# def calculate_score (sequences):
# creates scores matrix

# def create_guide_tree(sequences, matrix):
# creates guide tree

# def run_multiple_alignment(sequences, tree)
# uses tree to align the sequences

blosum62 = read_blosum("blossum_62.txt")

sequence1 = "HEAGAWGHEE"
sequence2 = "PAWHEAE"

test = pairwise_alignment(sequence1, sequence2, blosum62)
print(test)
