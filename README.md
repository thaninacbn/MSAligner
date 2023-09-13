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

Activate the environment

```
conda activate mseqalign
```


## Getting started 

To run the aligner in default mode, cd into the src directory and simply call from the command line:
````
python msaligner.py path/to/sequences.fasta path/to/output.txt
````

The output file will be a plain .txt file containing your aligned sequences: each line corresponds to a sequence!
