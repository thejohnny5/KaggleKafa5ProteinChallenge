# Kaggle Kafa 5 Protein Challenge
By: Jonathan Gray

## Challenge guidelines and objective
[Kaggle Cafa](https://www.kaggle.com/competitions/cafa-5-protein-function-prediction/data)

## Requirements:
NCBI Blast tools: [Latest Release](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
Python Version 3.9+

Preliminary Data Analysis performed in "Data Exploration.ipynb"

Modules:
src/orthologs -> Make blast protein database from list of fasta of proteins. Run the list against the database to find orthologs within the file.

src/models -> contains all machine learning models and params for testing.
