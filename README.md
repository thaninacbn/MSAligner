# MSAligner

A multiple sequence aligner written in python 🐍

Student project for the Python programming class of the Msc in Bioinformatics at the University of Paris Cité

MSAligner is a simple re-implementation of Clustal and uses dynamic programming and a guide UPGMA tree to 
quickly align protein sequences.

## Setting things up

Clone this repository in your local machine

```
git clone https://github.com/thaninacbn/MSAligner
```

Create a conda environment using the given .yml file

```
conda env create -f mseqalign.yml
```

And activate the environment

```
conda activate mseqalign
```


## Getting started 

Start by moving to the source code directory
```
cd src
```


### Sample alignment

You can try MSAligner on sample fasta files provided in the /data directory.

| ribosomes.fasta                                  | calmodulins.fasta                                   | wee1.fasta                                       | 
|--------------------------------------------------|-----------------------------------------------------|--------------------------------------------------|
| 50 AA                                            | 150 AA                                              | 600 AA                                           |   
| Homologs to human RL31 (accession number P62891) | Homologs to human CALM1 (accession number P0DP23)   | Homologs to human WEE1 (accession number P30291) |  
                                                                                                                                    
To run MSAligner on ribosomes.fasta, run the following command:

```
python msaligner.py ../data/ribosomes.fasta output_ribosomes.fasta
```

To try it on the other files, simply run:

```
python msaligner.py ../data/calmodulins.fasta output_calmo.fasta
```

or

```
python msaligner.py ../data/wee1.fasta output_wee1.fasta
```


### Run on your own sequences
To run MSAligner in default interactive mode with your own sequences, simply call from the command line:

````
python msaligner.py path/to/your/sequences.fasta output.txt
````

No matter the sequences you chose to use, a new /out directory will be automatically created and will contain the output.txt file.


## Using MSAligner as a module

You can also use the MSAligner functions by importing it as a module in your own python scripts 🐍

To do so start by importing MSAligner:

``` 
import msaligner as msa
```

There are seven functions you can call from MSAligner:

``read_fasta()`` and ``read_blosum()``: read a fasta file and a blosum matrix respectively. Both return dictionaries

``pairwise_alignment()``: takes 2 sequences as an input, a blossum matrix and asks you if you wish
to proceed to traceback (traceback="yes" or "no"). Performs the pairwise alignment of these two sequences
according to the Needleman and Wunsch algorithm.

``calculate_score()``: creates a scores matrix for all the pairs of sequences in a fasta file
by calling the pairwise_alignment() function

``turn_scores_into_distance()``: turns a scores matrix (contains negative and positive values)
to a distance matrix (contains only positive values)

``create_guide_tree()``: creates a guide tree for multiple alignment using the UPGMA method

``run_multiple_alignment()`` performs multiple sequence alignment on all proteic sequences from
a dictionary of sequences using an UPGMA guide tree