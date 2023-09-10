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


def pairwise_alignment(seq_m, seq_n, matrix, algt_score=False):
    """
    Reads a fasta file, gives each sequence an ID and stores them in a dictionary

    Parameters
    ----------
    seq_m: str
        First sequence to align

    seq_n: str
        Second sequence to align

    matrix: dict
        Dictionary containing scoring values from a blosum matrix read using the read_blosum() function

    algt_score: bool, optional



    Returns
    -------
    alignment_score: int
        Alignement score of the two sequences

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
            match = scores[i - 1, j - 1] + matrix[(seq_m[i - 1], seq_n[j - 1])]
            gap1 = scores[i - 1, j] + GAP_PENALTY
            gap2 = scores[i, j - 1] + GAP_PENALTY
            scores[i, j] = max(match, gap1, gap2)

    alignment_score = scores[m, n]

    # TODO: kmers ?

    # Traceback to find the optimal alignement
    # TODO : the traceback lol

    # me trying to be smart and have 1 function do multiple things
    if algt_score == True:
        return scores
    else:
        return alignment_score


def calculate_score(sequences, matrix):
    seq_ids = list(sequences.keys())
    size = len(seq_ids)

    # initialize the matrix with zeros
    scores_matrix = np.zeros((size, size))

    # iterate over the possible indexes
    for i in range(size):
        seq1 = sequences[i + 1]
        # j iterates between i+1 and size of matrix so that we only calculate the diagonal (limits compute time)
        for j in range(i + 1, size):
            seq2 = sequences[j + 1]
            alt_score = pairwise_alignment(seq1, seq2, matrix)
            scores_matrix[j, i] = alt_score
    print(scores_matrix)  # TODO delete when everything else works
    return scores_matrix


def turn_scores_into_distance(scores_matrix):
    dist_mat = np.zeros((scores_matrix.shape[0], scores_matrix.shape[1]))  # TODO: not zeros here actually
    maximum = np.nanmax(scores_matrix)
    minimum = np.nanmin(scores_matrix)  # Nanamin

    for i in range(scores_matrix.shape[0]):
        for j in range(i + 1, scores_matrix.shape[1]):
            dist_mat[j, i] = maximum - (scores_matrix[j, i] + minimum)

    print(dist_mat)  # TODO delete when everything else works
    return dist_mat


def create_guide_tree(sequences, dist_matrix):
    # UPGMA method hopefully

    seq_ids = list(sequences.keys())

    def join_rows(dist_matrix, a, b):
        #if b < a:
            #a, b = b, a

        row = []
        for i in range(0, a):
            row.append((dist_matrix[i][a] + dist_matrix[i][b]) / 2)
        dist_matrix = np.vstack([dist_matrix, row])

        for i in range(a + 1, b):
            dist_matrix[i, a] = (dist_matrix[i, a] + dist_matrix[b, i]) / 2

        for i in range(b + 1, dist_matrix.shape[1]):
            dist_matrix[i, a] = (dist_matrix[i, a] + dist_matrix[i, b]) / 2
            del dist_matrix[i, b]

        del dist_matrix[b]

    def get_structure(labels, a, b):
        #if b < a:
            #a, b = b, a

        labels[a]=(labels[a], labels[b])

        del labels[b]


    while dist_matrix.shape[0] > 1:
        coords = np.where(dist_matrix == np.nanmin(dist_matrix)) # there's smth wrong with the where here, idk what yet (maybe bc there's 0s in the matrix still)
        x = coords[0]
        y = coords[1]
        join_rows(dist_matrix,x,y)
        get_structure(seq_ids,x,y)

    return seq_ids[0]

    # im literally so sick rn that the very thought of coding this thing is making my headache 10 times worse
    # update from sometime later: i think this tree is LITERALLY the root of my sicknesses
    # (yes i am sick in plural)



# def run_multiple_alignment(sequences, tree)
# uses tree to align the sequences

blosum62 = read_blosum("blosum_62.txt")

my_seqs = read_fasta("Fichier_fasta.txt")
matrice = calculate_score(my_seqs, blosum62)
dist = turn_scores_into_distance(matrice)

test = create_guide_tree(my_seqs,dist)
print(test)
