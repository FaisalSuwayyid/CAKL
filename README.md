# CAKL: Commutative Algebra K-mer Learning

![Alt text](concept.png)

CAKL is an alignment-free method that computes k-mer algebraic representations of sequences based on k-mer analysis leveraging tools from commutative algebra, as detailed in the main manuscript. Though the code is designed for DNA sequences, the method can be generalized to arbitrary finite sequences.

Copyright © 2025 Faisal Suwayyid.

## Table of Contents

- [CAKL](#CAKL)
  - [Table of Contents](#table-of-contents)
  - [Description](#description)
  - [Basic Usage](basic-usage)
  - [Instructions for Usage](instructions-for-usage)
  - [License](#license)
  - [Citation](#citation)
  - [Required Packages](#required-packages)
  - [Respository File description](#respository-file-description)

## Description

CAKL is an alignment-free method that computes k-mer algebraic representations of sequences based on k-mer analysis leveraging tools from commutative algebra, as detailed in the main manuscript. Though the code is designed for DNA sequences, the method can be generalized to arbitrary finite sequences.

## Required Packages
The codes have been tested on macOS: Sonoma (14.6.1). The implementation is written in Python. The individual packages can be installed by following the instructions in their webpages, and are listed below:
```
# numpy==1.26.4
# scikit-learn ==1.4.2
# scipy==1.13.1
# gudhi==3.10.1
# biopython==1.84
# matplotlib==3.10.0
# pandas==2.2.2
# seaborn==0.13.2
# ete3==3.1.3
```

## Instructions for Usage
To use the code on your data, please ensure your input files are properly formatted:

1. **FASTA file**
   - Contains DNA sequences.
   - Each sequence should include its accession or accession with version.
   - If you are working with non-DNA data, you may need to modify the nucleotide bases in the code (psrt.py).

2. **CSV file**
   - Must contain a column named `Accession (version)` that exactly matches the accessions (or accession+version identifiers) in the FASTA file.  
   - The labels column can vary; adjust the column name in the code to match the one used in your CSV file.
  
3. **Featurization**
   - `Featurization.ipynb` reads sequences from the FASTA file and their corresponding accessions from the CSV file.  
   - Features are computed only for sequences whose accessions are listed in the CSV file, and in the same order as they appear in the CSV. 
   - For each sequence, the resulting feature array is stored and saved using the accession as the filename.
  
4. **Stacking**
   - `Stacking.ipynb` stacks the extracted features into a single feature matrix.  
   - Features are aligned in a consistent **sorted order** to ensure reproducibility.  

5. **Distances**
   - `Distances.ipynb` computes the pairwise distance matrix from the stacked features.  
   - The distance matrix follows the same **sorted order** of accessions used during stacking.  

6. **Applications**
   - **Classification**:  
     - `5nn_classification.ipynb` and `1nn_classification.ipynb` perform k-NN classification tasks (with k=5 and k=1, respectively).  
     - The resulting scores can be used to evaluate model performance.  
   - **Phylogenetic trees**:  
     - Phylogenetic trees can be constructed directly from the computed distance matrices.  



## Demo and reproduction

- **`example.py`**  
  Provides a minimal example for generating features.  
  Simply replace the placeholder sequence with your own sequence and run the script.  
  The code will generate the features for the provided sequence and save them in a file.
  The generated features include the f-, h-, and facet curves.
  You can also modify it as needed to suit your workflow.  

- **`example.ipynb`**  
  A complete example pipeline:  
  - Reads sequences from a FASTA file and their labels from a CSV file.  
  - Generates features for each sequence (saved as individual NumPy arrays).  
  - Stacks the features into a single feature matrix.  
  - Computes the distance matrix from the generated features.  
  - Constructs the tree, visualizes it as a dendrogram, and saves a Newick version.  
  - This example reproduces the tree construction used in the paper, and uses dendrogram visualization.  

## Repository File Description

- **`CAKL/`**  
  Contains the core source code and example workflows:  
  - `psrt.py`: Main implementation of the CAKL framework.  
  - `example.py`: Minimal usage example on arbitrary sequences.  
  - `example.ipynb`: Jupyter notebook demonstrating feature extraction, distance computation, and tree construction for sample datasets.  
  - `Featurization.ipynb`: Extracts features from sequences listed in a FASTA file and CSV metadata.  
  - `Stacking.ipynb`: Stacks individual feature arrays into a single matrix in sorted order for reproducibility.  
  - `Distances.ipynb`: Computes pairwise distance matrices from the stacked features.  
  - `1nn_classification.ipynb`: Performs leave-one-out 1-NN classification and reports evaluation metrics.  
  - `5nn_classification.ipynb`: Performs 5-fold 5-NN classification and reports evaluation metrics.  
  - `purity.ipynb`: Computes purity scores of phylogenetic trees.
  - `complete_workflow.ipynb`: Demonstrates the full workflow — **Featurization → Stacking → Distances** — on the NCBI datasets.  This notebook can then be extended with classification tasks for testing purposes.  Running it on the complete datasets, available in [link](https://weilab.math.msu.edu/DataLibrary/2D/Downloads/genetics.zip) allows full reproduction of the results presented in the paper.  


- **`CAKL/data/`**  
  Includes accession lists of NCBI datasets used in the main manuscript. It also includes samples of them to run the code on. The complete datasets are available in [link](https://weilab.math.msu.edu/DataLibrary/2D/Downloads/genetics.zip).

- **`CAKL/data2/`**  
  Includes datasets used for phylogenetic tree construction and genetic identification.  

- **`CAKL/trees_purity/`**  
  Contains phylogenetic trees generated from the datasets, used for purity computations.  

## License

CAKL is licensed under the MIT license (COPYING.txt), with an extra clause (CONTRIBUTING.txt) clarifying the license for modifications released without an explicit written license agreement.

## Citation

If you wish to cite this work, please use the following citation:
```
@misc{suwayyid2025cakl,
      title={CAKL: Commutative Algebra $k$-mer Learning of Genomics}, 
      author={Faisal Suwayyid and Yuta Hozumi and Hongsong Feng and Mushal Zia and JunJie Wee and Guo-Wei Wei},
      year={2025},
      eprint={2508.09406},
      archivePrefix={arXiv},
      url={https://arxiv.org/abs/2508.09406}, 
}
```
