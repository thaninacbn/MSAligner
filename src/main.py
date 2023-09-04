# define main functions

def read_fasta(fasta_file):
    #code from previous python exos
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
    return sequence




def pairwise_alignment(seq1, seq2):
    # runs pairwise alignment using needleman-wunsch algorithm

def calculate_score (sequences):
    # creates a scores matrix

def create_guide_tree(sequences, matrix):
    # creates guide tree

def run_multiple_alignment(sequences, tree)
    # uses tree to align the sequences

