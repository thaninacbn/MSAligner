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
    # TODO: deal with the fact that there might be more than 1 sequence per file: dictionary
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


def pairwise_alignment(seq_m, seq_n, matrix):
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
    print(scores)

    # find the maximum alignment score in the matrix
    list_scores = scores[m, 0:n]+scores[0:m, n]
    alignment_score = max(list_scores)
    print(alignment_score)
    # TODO: kmers ?

    return scores


# def calculate_score (sequences):
# creates scores matrix


# def create_guide_tree(sequences, matrix):
# creates guide tree

# def run_multiple_alignment(sequences, tree)
# uses tree to align the sequences

blosum62 = read_blosum("blossum_62.txt")
#print(blosum62)

sequence1 = "HEAGAWGHEE"
sequence2 = "PAWHEAE"

test = pairwise_alignment(sequence1, sequence2, blosum62)
print(test)
