# CAKL

CAKL: Commutative Algebra K-mer Learning

CAKL is an alignment-free method that computes k-mer algebraic representations of sequences based on k-mer analysis leveraging tools from commutative algebra, as detailed in the main manuscript. Though the code is designed for DNA sequences, the method can be generalized to arbitary finite sequences.

Copyright © 2025 Faisal Suwayyid.

## Table of Contents

- [CAKL](#CAKL)
  - [Table of Contents](#table-of-contents)
  - [Description](#description)
  - [Basic Usage](basic-usage)
  - [License](#license)
  - [Citation](#citation)
  - [Required Packages](#required-packages)
  - [Respository File description](#respository-file-description)

## Description

CAKL is an alignment-free method that computes k-mer algebraic representations of sequences based on k-mer analysis leveraging tools from commutative algebra, as detailed in the main manuscript. Though the code is designed for DNA sequences, the method can be generalized to arbitary finite sequences.

## Basic Usage

- **`example.py`**  
  Provides a minimal example for generating features.  
  Simply replace the placeholder sequence with your own sequence and run the script.  
  You can also modify it as needed to suit your workflow.  

- **`example.ipynb`**  
  A complete example pipeline:  
  - Reads sequences from a FASTA file and their labels from a CSV file.  
  - Generates features for each sequence (saved as individual NumPy arrays).  
  - Stacks the features into a single feature matrix.  
  - Computes the distance matrix from the generated features.
  - Construct the tree, view it, and save a Newick version of the tree.

## License

CAKL is licensed under the MIT license (COPYING.txt), with an extra clause (CONTRIBUTING.txt) clarifying the license for modifications released without an explicit written license agreement.

## Citation

When the paper for this software is available on the Arxiv and/or published, we will provide an appropriate bibtex entry for those who would like to cite this software.

## Required Packages
The inidivudal packages that we utilized is listed below:
```
numpy
scikit-learn
scipy
gudhi
biopython
matplotlib
pandas
```

## Repository File Description

- **`CAKL/`**  
  Contains the core source code and examples:  
  - `psrt.py`: main implementation of the CAKL framework.  
  - `example.py`: usage example on arbitrary sequences.  
  - `example.ipynb`: Jupyter notebook demonstrating how to reproduce feature extraction and distance computations for two sample datasets.  

- **`data2/`**  
  Includes two example datasets used in the experiments and demonstrations.  

