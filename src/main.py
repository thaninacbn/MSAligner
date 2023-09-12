import argparse
import numpy as np

# define main functions


GAP_PENALTY = -5


def read_fasta(fasta_file):
    """
    Reads a fasta file, gives each sequence an ID and stores them in a dictionary

    Parameters
    ----------
    fasta_file: str
        The path to the file containing the sequences (in fasta format)

    Returns
    -------
    sequence_dict : dict
        Dictionary with key = sequence ID (int, attributed from 1 to n= number of sequences) and value = sequence

    """

    seq_id_counter = 0
    # will store sequences in a dictionary
    sequence_dict = {}
    with open(fasta_file, "r") as filin:
        seq_id = ""
        for line in filin:
            if line.startswith(">"):
                # assigns a numerical ID to every sequence (for easy accession later)
                seq_id_counter += 1
                seq_id = seq_id_counter
                sequence_dict[seq_id] = ""
            else:
                sequence_dict[seq_id] += line.strip()
    return sequence_dict


def read_blosum(matrix_file):
    """
    Reads a txt file containing a blosum matrix and returns a dictionary of all the amino-acid match scores

    Parameters
    ----------
    matrix_file: str
        The path to the file containing the blosum matrix

    Returns
    -------
    blosum_dict : dict
        Dictionary with key = tuple containing a pair of amino acids and value = match score

    """
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


def pairwise_alignment(seq_m, seq_n, blosum_matrix, traceback=False):
    """
    Performs a global alignment on two sequences using the Needleman and Wunsch algorithm.

    Parameters
    ----------
    seq_m: str
        First sequence to align

    seq_n: str
        Second sequence to align

    blosum_matrix: dict
        Dictionary containing scoring values from a blosum matrix read using the read_blosum() function

    traceback: bool, optional
        If set to true, returns the alignment (as a list of strings) of the 2 sequences
        Default = false


    Returns
    -------
    alignment_score: int
        Alignement score of the two sequences

    alignmed_seq : list
        List containing the aligned sequences (seq_m as index 0 and seq_n as index 1)
        returned only if argument traceback is set to True

    """

    # m rows, n columns as per convention
    m = len(seq_m)
    n = len(seq_n)

    # generates a matrix of size n+1 m+1 to store gap penalty
    scores = np.zeros((m + 1, n + 1))

    # fill the first row and column with the gap penalty
    scores[:, 0] = [GAP_PENALTY * i for i in range(scores.shape[0])]
    scores[0, :] = [GAP_PENALTY * i for i in range(scores.shape[1])]

    # fill the score matrix
    for i in range(1, scores.shape[0]):
        for j in range(1, scores.shape[1]):
            match = scores[i - 1, j - 1] + blosum_matrix[(seq_m[i - 1], seq_n[j - 1])]
            gap1 = scores[i - 1, j] + GAP_PENALTY
            gap2 = scores[i, j - 1] + GAP_PENALTY
            scores[i, j] = max(match, gap1, gap2)

    alignment_score = scores[m, n]

    if traceback == False:
        return alignment_score
    else:

        # Traceback to find the optimal alignment

        align_m = ""
        align_n = ""

        # starting point in the bottom right cell in the matrix
        i = m
        j = n

        while i > 0 and j > 0:
            current_score = scores[i, j]
            diagonal_score = scores[i - 1, j - 1]
            upwards_score = scores[i, j - 1]
            left_score = scores[i - 1, j]

            # check which cell the current score comes from, update i and j accordingly
            if current_score == diagonal_score + blosum_matrix[(seq_m[i - 1], seq_n[j - 1])]:
                align_m += seq_m[i - 1]
                align_n += seq_n[j - 1]
                i -= 1
                j -= 1
            elif current_score == upwards_score + GAP_PENALTY:
                align_n += seq_n[j - 1]
                align_m += '-'
                j -= 1
            elif current_score == left_score + GAP_PENALTY:
                align_n += '-'
                align_m += seq_m[i - 1]
                i -= 1

        while j > 0:
            align_n += seq_n[j - 1]
            align_m += '-'
            j -= 1
        while i > 0:
            align_n += '-'
            align_m += seq_m[i - 1]
            i -= 1

        align_m = align_m[::-1]
        align_n = align_n[::-1]

        aligned_seq = [align_m, align_n]

        return aligned_seq


def calculate_score(sequences, blosum_matrix):
    """
    Calculates the pairwise alignment scores of all non-redundant sequence pairs and stores them into a 2D numpy array.
    Calls for the pairwise_alignment() function.

    Parameters
    ----------
    sequences: dict
        A dictionnary that contains the sequences to align (value) and their ID (key)

    blosum_matrix: dict
        A dictionary that contains the match-scores of all amino-acids from a blosum matrix (usually the output of the
        read_blosum function)

    Returns
    -------
    scores_matrix: matrix
        A 2D numpy array that contains the alignment scores of all possible non-redundant pairs of sequences on the
        lower triangular half, and np.nan values on the upper triangular half.

    """

    seq_ids = list(sequences.keys())
    size = len(seq_ids)

    # initialize the matrix with zeros
    scores_matrix = np.empty((size, size))
    scores_matrix.fill(np.nan)

    # iterate over the possible indexes
    for i in range(size):
        seq1 = sequences[i + 1]
        # j iterates between i+1 and size of matrix so that we only calculate the diagonal (limits compute time)
        for j in range(i + 1, size):
            seq2 = sequences[j + 1]
            alt_score = pairwise_alignment(seq1, seq2, blosum_matrix)
            scores_matrix[j, i] = alt_score
    print(scores_matrix)  # TODO delete when everything else works
    return scores_matrix


def turn_scores_into_distance(scores_matrix):
    """
    Turns a scores matrix into a distances matrix

    Parameters
    ----------
    scores_matrix: matrix
        A numpy 2D matrix that contains the pairwise alignment scores of the sequences on the lower triangular half.
        Usually the result of calling the function calculate_score().

    Returns
    ----------
    dist_mat : matrix
        A 2D numpy matrix containing the distance between the sequences in the lower triangular half and np.nan in
        the upper triangular half.

    """

    # initialize the matrix to the same size as the scores matrix
    dist_mat = np.empty((scores_matrix.shape[0], scores_matrix.shape[1]))
    dist_mat.fill(np.nan)
    # set the min and max values from the scores matrix
    maximum = np.nanmax(scores_matrix)
    minimum = np.nanmin(scores_matrix)  # Nanamin

    # fill the distance matrix
    for i in range(scores_matrix.shape[0]):
        for j in range(i + 1, scores_matrix.shape[1]):
            dist_mat[j, i] = maximum - (scores_matrix[j, i] + minimum)

    print(dist_mat)  # TODO delete when everything else works
    return dist_mat


def create_guide_tree(sequences, dist_matrix):
    """
    Takes a distance matrix and a dictionary of sequences and creates a guide tree using the UPGMA method.

    Parameters
    ----------
    sequences: dict
        A dictionnary that contains the sequences to align (value) and their ID (key)

    dist_matrix: matrix
        A numpy 2D array containing the distances between each sequence. This matrix should be lower-diagonal. Usually
        the output of the turn_scores_into_distances() function.

    Returns
    -------
    tree_structure : tuple
        The structure of the UPGMA guide tree in the form of a tuple of tuples.

    """

    seq_ids = list(sequences.keys())
    distances = dist_matrix.tolist()  # literally lost my mind and decided to drop the np arrays for the time being

    def get_lowest_value(matrix):
        # est ce que je dois faire des docstring pour ces inner fonctions ? vu qu'elles sont pas accessible
        # set the default lowest possible value to some ridiculously high value
        min_cell = float("inf")
        x, y = -1, -1

        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                # updates the coordinates if the value in the cell is lower than the one currently stored
                if matrix[i][j] < min_cell:
                    min_cell = matrix[i][j]
                    x, y = i, j

        return x, y

    def merge_seq_ids(seq_list, a, b):
        if b < a:
            a, b = b, a

        # orders the sequences in a tuple
        seq_list[a] = (seq_list[a], seq_list[b])
        del seq_list[b]

    def reduce_matrix(matrix, a, b):
        # index out of range error if this isn't there
        if b < a:
            a, b = b, a

        # updates all the coordinates with the arithmetic mean of the a and b indices
        new_row = []
        for i in range(0, a):
            new_row.append((matrix[a][i] + matrix[b][i]) / 2)
        matrix[a] = new_row

        for i in range(a + 1, b):
            matrix[i][a] = (matrix[i][a] + matrix[b][i]) / 2

        for i in range(b + 1, len(matrix)):
            matrix[i][a] = (matrix[i][a] + matrix[i][b]) / 2
            del matrix[i][b]

        del matrix[b]

    def run_upgma(matrix, seq_list):

        while len(seq_list) > 1:
            x, y = get_lowest_value(matrix)
            reduce_matrix(matrix, x, y)
            merge_seq_ids(seq_list, x, y)
        return seq_list[0]

    tree_structure = run_upgma(distances, seq_ids)
    return tree_structure

    # im literally so sick rn that the very thought of coding this thing is making my headache 10 times worse
    # update from sometime later: i think this tree is LITERALLY the root of my sicknesses
    # (yes i am sick in plural)


def run_multiple_alignment(sequence_dict, tree, blosum_matrix):
    # flattens the tree to reduce complexity bc i don't have time for this anymore
    def flatten_tree(tree_tuple):
        tree_list = []
        for item in tree_tuple:
            if isinstance(item, tuple):
                flatten_tree(item)
            else:
                tree_list.append(item)
        return tree_list

    flat_tree = flatten_tree(tree)

    index_m = flat_tree[0]
    index_n = flat_tree[1]
    seq_m = sequence_dict[index_m]
    seq_n = sequence_dict[index_n]

    start_align = pairwise_alignment(seq_m, seq_n, blosum_matrix, traceback=True)

    for i in range(2, len(flat_tree)):
        seq_n = sequence_dict[i]


    return 1


# uses tree to align the sequences

blosum62 = read_blosum("blosum_62.txt")

my_seqs = read_fasta("melanie.fasta")
mat_scores = calculate_score(my_seqs, blosum62)

sequence1 = "HEAGAWGHEEPAH"
sequence2 = "PAWHEAE"

# test = pairwise_alignment(sequence1, sequence2, blosum62)
