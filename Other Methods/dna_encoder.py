

## Alfpy

import numpy as np

def base_base_correlation(seq, k, alphabet=None):
    """
    Compute the base-base correlation (BBC) vector for a sequence.

    Args:
        seq (str): DNA sequence.
        k (int): Maximum distance to observe correlation between bases.
        alphabet (str or list, optional): List of possible characters. 
            Used to filter the sequence if provided.

    Returns:
        numpy.ndarray: BBC vector of shape (1, L*L), where L is the size of the alphabet.

    Raises:
        ValueError: If the sequence is too short for the given k.

    Examples:
        >>> print(base_base_correlation('ATGCATGC', 1, 'ATGC'))
        [[-0.12547302 -0.12547302  0.2281059   0.17169665  0.01815213
          -0.12547302 -0.12547302  0.04258163  0.04258163  0.17169665
          -0.12547302 -0.12547302 -0.12547302  0.2281059   0.17169665
          -0.12547302]]
    """
    s = seq.upper()

    if k > len(s) - 2:
        raise ValueError(f"Sequence too short to compute BBC with k={k}")

    if alphabet is None:
        alphabet = sorted(set(s))
    else:
        # Ensure alphabet is a sorted list for consistent ordering
        alphabet = sorted(alphabet)
        # Filter the sequence to include only characters in the alphabet
        s = ''.join([c for c in s if c in alphabet])

    L = len(alphabet)
    if L == 0:
        raise ValueError("Alphabet is empty after filtering the sequence.")

    # Create a mapping from nucleotide to index
    nuc_to_idx = {nuc: idx for idx, nuc in enumerate(alphabet)}

    # Compute the base probabilities for each character
    p = np.zeros(L)
    for c in s:
        p[nuc_to_idx[c]] += 1
    p_sum = np.sum(p)
    if p_sum == 0:
        raise ValueError("No valid nucleotides found in the sequence after filtering.")
    p /= p_sum
    p = p.reshape(1, L)  # Shape (1, L)

    bbc = np.zeros((L, L))

    for l in range(1, k + 2):
        l_dist_correlations = np.zeros((L, L))
        for i in range(len(s) - l):
            nuc1 = nuc_to_idx[s[i]]
            nuc2 = nuc_to_idx[s[i + l]]
            l_dist_correlations[nuc1][nuc2] += 1
        total_pairs = np.sum(l_dist_correlations)
        if total_pairs == 0:
            continue  # Avoid division by zero if no pairs are found for this distance
        l_dist_correlations /= total_pairs

        # Compute the deviation from statistical independence
        D = l_dist_correlations - np.dot(p.T, p)

        # Accumulate the correlation measures
        bbc += D + (D ** 2 / 2 * (p.T ** 2) * (p ** 2)) + (D ** 3)

    # Flatten the BBC matrix into a vector
    bbc = bbc.flatten().reshape(1, L * L)

    return bbc


import numpy as np

def fcgr_vector(dnaseq, word_size):
    """
    Create a FCGR (Frequency Chaos Game Representation) vector representing a DNA sequence.

    Args:
        dnaseq (str or list): DNA sequence composed of characters 'A', 'T', 'G', 'C'.
        word_size (int): Word size (>= 1). Determines the resolution of the FCGR.
                         The resulting vector will have a length of 4^word_size.

    Returns:
        list: FCGR vector with length equal to 4^word_size.

    Examples:
        >>> s = 'ATGCTGATGGATG'
        >>> print(fcgr_vector(s, 1))
        [5, 3, 5, 0]

        >>> print(fcgr_vector(s, 2))
        [1, 0, 1, 0, 0, 0, 4, 0, 2, 2, 0, 0, 1, 3, 0, 0]
    """
    # Validate word_size
    if word_size < 1:
        raise ValueError("word_size must be at least 1.")

    # Initialize parameters
    ndata = 4 ** word_size  # Total number of k-mers
    genlen = len(dnaseq)
    if genlen == 0:
        raise ValueError("DNA sequence is empty.")

    # Define CGR points for each nucleotide
    Apoint = np.array([0.0, 1.0])
    Tpoint = np.array([1.0, 1.0])
    Gpoint = np.array([1.0, 0.0])
    Cpoint = np.array([0.0, 0.0])

    # Initialize CGR matrix
    CGRs = np.zeros((genlen + 1, 2))
    CGRs[0] = np.array([0.5, 0.5])  # Start at the center

    # Populate CGR points
    for i in range(genlen):
        nucleotide = dnaseq[i].upper()
        if nucleotide == 'A':
            CGRs[i + 1] = 0.5 * (CGRs[i] + Apoint)
        elif nucleotide == 'T':
            CGRs[i + 1] = 0.5 * (CGRs[i] + Tpoint)
        elif nucleotide == 'G':
            CGRs[i + 1] = 0.5 * (CGRs[i] + Gpoint)
        elif nucleotide == 'C':
            CGRs[i + 1] = 0.5 * (CGRs[i] + Cpoint)
        else:
            # Handle unknown nucleotides by keeping the previous point
            CGRs[i + 1] = CGRs[i]

    # Calculate grid resolution
    grid_size = 2 ** word_size
    cell_size = 1.0 / grid_size

    # Initialize FCGR vector
    vectors = [0.0] * ndata

    # Map CGR points to grid cells
    for point in CGRs[1:]:  # Exclude the initial point
        xx = int(point[0] / cell_size)
        yy = int(point[1] / cell_size)

        # Handle boundary conditions
        if xx >= grid_size:
            xx = grid_size - 1
        if yy >= grid_size:
            yy = grid_size - 1

        index = yy * grid_size + xx
        vectors[index] += 1

    return vectors


import numpy as np

def a_2DSGraphVector(seq):
    """Create 10-dimensional statistical vector to characterize a DNA sequence.

    Args:
        seq (str/list): DNA sequence

    Returns:
        numpy.ndarray (10,)

    Examples:
        >>> s = 'ATGCCTGACTGNATGAGAGAC'
        >>> print(_2DSGraphVector(s))
        [  2.31 -11.83  12.    -4.89   7.75   0.74  11.5   -3.85  9.5  -2.35]
    """
    nucleotides = [
        ("A", -(3**0.3333333)),
        ("T", 2**0.3333333),
        ("G", -(5**0.5)),
        ("C", 3**0.5)
    ]

    points = np.zeros((len(seq), 2))
    temppoint = np.array([0.0, 0.0])

    d = {}
    for nt, val in nucleotides:
        d[nt] = [np.array([1.0, val]), 0.0, 0.0, 0]

    i = 0
    for nt in seq:
        # Include only recognized nucleotides
        if nt in d:
            points[i] = d[nt][0] + temppoint
            d[nt][1] += points[i, 0]
            d[nt][2] += points[i, 1]
            d[nt][3] += 1
            temppoint = points[i]
            i += 1

    if i == 0:
        # Handle case with no valid nucleotides
        return np.zeros(10)

    peak = points[:i].max(axis=0)[1]
    lowest = points[:i].min(axis=0)[1]
    l = [peak, lowest]

    for nt, _ in nucleotides:
        count = d[nt][3]
        if count > 0:
            avg_x = d[nt][1] / count
            avg_y = d[nt][2] / count
        else:
            avg_x = 0.0
            avg_y = 0.0
        l.extend([avg_x, avg_y])

    return np.array(l)



def a_2DMGraphVector(seq, n):
    """Create n-dimensional moment vector to characterize a DNA sequence.

    Args:
        seq (str/list): DNA sequence
        n (int): Number of moments

    Returns:
        numpy.ndarray (n,)

    Examples:
        >>> s = 'ATGCCTGACTGNATGAGAGAC'
        >>> print(_2DMGraphVector(s, 10))
        [21. 13.44  13.44 16.15  21.16   29.01  40.87  58.58  84.99  124.39]

        >>> print(_2DMGraphVector(s, 5))
        [ 21.    13.44  13.44  16.15  21.16]
    """
    nucleotides = [
        ("A", -(3**0.3333333)),
        ("T", 2**0.3333333),
        ("G", -(5**0.5)),
        ("C", 3**0.5)
    ]

    points = np.zeros((len(seq), 2))
    temppoint = np.array([0.0, 0.0])

    d = {}
    for nt, val in nucleotides:
        d[nt] = np.array([1.0, val])

    i = 0
    for nt in seq:
        # Include only recognized nucleotides
        if nt in d:
            points[i] = d[nt] + temppoint
            temppoint = points[i]
            i += 1
    seqlen = i

    if seqlen == 0:
        return np.zeros(n)

    l = []
    for k in range(n):
        if k == 0:
            # Define the 0th moment as the sequence length
            v = seqlen
        else:
            # Calculate the k-th moment
            diff = points[:seqlen, 0] - points[:seqlen, 1]
            v = np.sum(diff**k)
            if seqlen > 0:
                v /= pow(seqlen, k)
            else:
                v = 0.0
        l.append(v)
    return np.array(l)

def a_2DNGraphVector(genomeseq):
    """Create 48-dimensional natural vector to characterize a DNA sequence.

    Args:
        genomeseq (str/list): DNA sequence

    Returns:
        numpy.ndarray (48,)

    Examples:
        >>> s = 'ATGCCTGACTGNATGAGAGAC'
        >>> print(_2DNGraphVector(s))
        [  6.00e+00   4.00e+00   6.00e+00   4.00e+00   1.17e+01   7.00e+00
           1.10e+01   8.75e+00   1.99e+00   9.52e-01   1.51e+00   2.18e+00
          -5.52e-01  -1.76e-01  -4.70e-01  -3.18e-01   5.73e-02   1.66e-02
           4.52e-02   4.10e-02  -5.12e-03  -1.35e-03  -3.82e-03  -4.10e-03
           4.77e-04   1.13e-04   3.35e-04   4.30e-04  -4.41e-05  -9.38e-06
          -2.92e-05  -4.47e-05   4.08e-06   7.81e-07   2.55e-06   4.66e-06
          -3.78e-07  -6.51e-08  -2.23e-07  -4.85e-07   3.50e-08   5.43e-09
           1.94e-08   5.05e-08  -3.24e-09  -4.52e-10  -1.70e-09  -5.26e-09]
    """
    genlen = len(genomeseq)

    nts = 'ATGC'

    idx = {nt: [] for nt in nts}

    for i, nt in enumerate(genomeseq):
        if nt in idx:
            idx[nt].append(i)

    counts = {nt: len(idx[nt]) for nt in nts}
    vector = [counts[nt] for nt in nts]

    means = {}
    for nt in nts:
        if counts[nt] > 0:
            means[nt] = np.mean(idx[nt])
        else:
            means[nt] = 0.0  # Assign default mean

    for nt in nts:
        vector.append(means[nt])

    templist = {nt: idx[nt].copy() for nt in nts}

    for k in range(2, 12):
        for nt in nts:
            if counts[nt] > 0:
                for i in range(counts[nt]):
                    denominator = pow(counts[nt] * genlen, k - 1)
                    if denominator != 0:
                        tempn = pow((idx[nt][i] - means[nt]), k) / denominator
                    else:
                        tempn = 0.0
                    templist[nt][i] = tempn
                moment_sum = np.sum(templist[nt])
            else:
                moment_sum = 0.0
            vector.append(moment_sum)
    return np.array(vector)











def complexity(s):
    """Calculate the Lempel-Ziv complexity c of a given sequence.

    As described in:
    Kaspar F, Schuster HG. Phys Rev A. 1987 36(2):842-848.
    doi: http://dx.doi.org/10.1103/PhysRevA.36.842

    Based on:
    http://stackoverflow.com/a/30694008

    Args:
        s (str/list): Sequence of any characters.

    Returns:
        int: Lempel-Ziv complexity.
    """
    if not s:
        return 0  # No complexity for empty sequence

    i, k, l = 0, 1, 1
    k_max = 1
    n = len(s)
    c = 1
    while l + k <= n:
        if s[i + k - 1] == s[l + k - 1]:
            k += 1
            if l + k > n:
                c += 1
                break
        else:
            if k > k_max:
                k_max = k
            i += 1
            if i == l:
                c += 1
                l += k_max
                if l >= n:
                    break
                i = 0
                k = 1
                k_max = 1
            else:
                k = 1
    return c

def is_reproducible(seq, index, hist_len, last_match=0):
    """Check whether a substring is reproducible within the sequence.

    A sequence R is reproducible from sequence S (denoted S > R) when
    R can be obtained from S by copying elements from p-th location in
    S to the end of S. For example, AACGT > AACGTCGTCG with p = 3 and
    AACGT > AACGTAC with p = 2.

    Based on:
    Hohl M, Ragan M. Systematic Biology. 2007. 56(2):206-21

    Args:
        seq (str/list): The entire sequence.
        index (int): Current position in the sequence to check for reproducibility.
        hist_len (int): Length of the history to consider for matching.
        last_match (int, optional): The starting position for searching matches. Defaults to 0.

    Returns:
        tuple:
            bool: True if reproducible, False otherwise.
            int: Position of the matching substring if reproducible, else 0.
    """
    hist_start = index - hist_len
    if hist_start < 0:
        # Not enough history to match
        return False, 0

    # Limit the search to the history segment
    for i in range(last_match, hist_start + 1):
        # Ensure that the comparison does not exceed sequence bounds
        if hist_start + hist_len > len(seq):
            return False, 0

        # Compare the substring in history with the current substring
        if seq[i:i + hist_len] == seq[hist_start:hist_start + hist_len]:
            return True, i  # Reproducible, return the match position

    return False, 0

def complexity1(seq):
    """Calculate the Lempel-Ziv complexity c of a given sequence.

    As described in:
    Otu HH, Sayood K. Bioinformatics. 2003 19(16):2122-30.
    PubMed PMID: 14594718.

    Based on:
    Hohl M, Ragan M. Systematic Biology. 2007. 56(2):206-21

    Args:
        seq (str/list): Sequence of any characters.

    Returns:
        int: Lempel-Ziv complexity.
    """
    start = 2
    length = len(seq)
    if length < start:
        # Complexity equals length for very short sequences.
        return length

    complexity = 2  # Starting complexity
    history_length = 1  # Length of the current history segment
    last_match = 0  # Position of the last match

    for index in range(start, length):
        is_reproducible_flag, last_match = is_reproducible(
            seq, index, history_length, last_match)
        if is_reproducible_flag:
            history_length += 1
        else:
            history_length = 1
            complexity += 1

    return complexity



############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################





import itertools
from typing import Dict, List
import pandas as pd


class kmer_counts_encoder:
    """
    A class to perform k-mer analysis on DNA sequences.

    Attributes:
        k (int): The length of k-mers.
        all_kmers (List[str]): List of all possible k-mers.
        counts (Dict[str, int]): Dictionary mapping k-mers to their counts in the DNA sequence.
        normalized_counts (Dict[str, float]): Dictionary mapping k-mers to their normalized frequencies.
    """

    def __init__(self, k: int):
        """
        Initializes the KmerAnalyzer with a specific k-mer length.

        Parameters:
            k (int): The length of k-mers.

        Raises:
            ValueError: If k is not a positive integer or if k is too large.
        """
        if not isinstance(k, int) or k <= 0:
            raise ValueError("k must be a positive integer.")
        if 4**k > 1e7:
            raise ValueError("k is too large. The number of possible k-mers exceeds 10 million.")

        self.k = k
        self.all_kmers = self.generate_kmers()
        self.counts: Dict[str, int] = {}
        # self.normalized_counts: Dict[str, float] = {}

    def generate_kmers(self) -> List[str]:
        """
        Generate all possible k-mers of length k using DNA letters.

        Returns:
            List[str]: A list of all possible k-mers.
        """
        return [''.join(p) for p in itertools.product('ACGT', repeat=self.k)]

    def count_kmers(self, dna: str) -> Dict[str, int]:
        """
        Count occurrences of each k-mer in a given DNA sequence, assigning zero to those not appearing.

        Parameters:
            dna (str): The DNA sequence to analyze.

        Returns:
            Dict[str, int]: A dictionary mapping each k-mer to its count.
        """
        dna = dna.upper()
        self.counts = {kmer: 0 for kmer in self.all_kmers}

        # Slide a window of size k across the DNA sequence
        for i in range(len(dna) - self.k + 1):
            kmer = dna[i:i + self.k]
            if kmer in self.counts:
                self.counts[kmer] += 1
            # If kmer contains invalid characters, it's ignored

        return self.counts

    def ffp_normalize_counts(self, dna_length: int) -> Dict[str, float]:
        """
        Normalize k-mer counts by DNA length - (k - 1).

        Parameters:
            dna_length (int): The length of the DNA sequence.

        Returns:
            Dict[str, float]: A dictionary mapping each k-mer to its normalized frequency.
        """
        if not self.counts:
            raise ValueError("Counts have not been calculated. Please run count_kmers() first.")

        total_positions = dna_length - self.k + 1
        if total_positions <= 0:
            raise ValueError("DNA length is shorter than k-mer length.")

        # self.normalized_counts = {kmer: count / total_positions for kmer, count in self.counts.items()}
        return self.counts

    def analyze(self, dna: str) -> Dict[str, float]:
        """
        Perform the full k-mer analysis: count k-mers and normalize their counts.

        Parameters:
            dna (str): The DNA sequence to analyze.

        Returns:
            Dict[str, float]: A dictionary mapping each k-mer to its normalized frequency.
        """
        self.count_kmers(dna)
        self.ffp_normalize_counts(len(dna))
        return self.counts

    def get_kmer_counts(self) -> Dict[str, int]:
        """
        Retrieve the k-mer counts.

        Returns:
            Dict[str, int]: A dictionary mapping each k-mer to its count.
        """
        if not self.counts:
            raise ValueError("Counts have not been calculated. Please run count_kmers() first.")
        return self.counts

    def get_normalized_counts(self) -> Dict[str, float]:
        """
        Retrieve the normalized k-mer counts.

        Returns:
            Dict[str, float]: A dictionary mapping each k-mer to its normalized frequency.
        """
        if not self.normalized_counts:
            raise ValueError("Normalized counts have not been calculated. Please run ffp_normalize_counts() first.")
        return self.normalized_counts

    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert the k-mer counts and normalized counts to a Pandas DataFrame.

        Returns:
            pd.DataFrame: A DataFrame containing k-mers, their counts, and normalized frequencies.
        """
        if not self.counts or not self.normalized_counts:
            raise ValueError("Counts and normalized counts have not been calculated. Please run analyze() first.")

        data = {
            'k-mer': self.all_kmers,
            'Count': [self.counts[kmer] for kmer in self.all_kmers],
            'Normalized Frequency': [self.normalized_counts[kmer] for kmer in self.all_kmers]
        }
        return pd.DataFrame(data)

    def display_top_kmers(self, top_n: int = 10) -> pd.DataFrame:
        """
        Display the top N k-mers based on their counts.

        Parameters:
            top_n (int): The number of top k-mers to display.

        Returns:
            pd.DataFrame: A DataFrame containing the top N k-mers.
        """
        df = self.to_dataframe()
        top_kmers = df.sort_values(by='Count', ascending=False).head(top_n)
        return top_kmers.reset_index(drop=True)

    def export_to_csv(self, filepath: str):
        """
        Export the k-mer analysis results to a CSV file.

        Parameters:
            filepath (str): The path to the output CSV file.
        """
        df = self.to_dataframe()
        df.to_csv(filepath, index=False)

    def plot_kmer_frequencies(self, top_n: int = 20):
        """
        Plot the frequencies of the top N k-mers using a bar chart.

        Parameters:
            top_n (int): The number of top k-mers to plot.
        """
        import matplotlib.pyplot as plt

        top_kmers = self.display_top_kmers(top_n)
        plt.figure(figsize=(10, 6))
        plt.bar(top_kmers['k-mer'], top_kmers['Normalized Frequency'], color='skyblue')
        plt.xlabel('k-mer')
        plt.ylabel('Normalized Frequency')
        plt.title(f'Top {top_n} k-mer Frequencies')
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.show()








############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################




import itertools
from typing import Dict, List
import pandas as pd


class FFP_encoder:
    """
    A class to perform k-mer analysis on DNA sequences.

    Attributes:
        k (int): The length of k-mers.
        all_kmers (List[str]): List of all possible k-mers.
        counts (Dict[str, int]): Dictionary mapping k-mers to their counts in the DNA sequence.
        normalized_counts (Dict[str, float]): Dictionary mapping k-mers to their normalized frequencies.
    """

    def __init__(self, k: int):
        """
        Initializes the KmerAnalyzer with a specific k-mer length.

        Parameters:
            k (int): The length of k-mers.

        Raises:
            ValueError: If k is not a positive integer or if k is too large.
        """
        if not isinstance(k, int) or k <= 0:
            raise ValueError("k must be a positive integer.")
        if 4**k > 1e7:
            raise ValueError("k is too large. The number of possible k-mers exceeds 10 million.")

        self.k = k
        self.all_kmers = self.generate_kmers()
        self.counts: Dict[str, int] = {}
        self.normalized_counts: Dict[str, float] = {}

    def generate_kmers(self) -> List[str]:
        """
        Generate all possible k-mers of length k using DNA letters.

        Returns:
            List[str]: A list of all possible k-mers.
        """
        return [''.join(p) for p in itertools.product('ACGT', repeat=self.k)]

    def count_kmers(self, dna: str) -> Dict[str, int]:
        """
        Count occurrences of each k-mer in a given DNA sequence, assigning zero to those not appearing.

        Parameters:
            dna (str): The DNA sequence to analyze.

        Returns:
            Dict[str, int]: A dictionary mapping each k-mer to its count.
        """
        dna = dna.upper()
        self.counts = {kmer: 0 for kmer in self.all_kmers}

        # Slide a window of size k across the DNA sequence
        for i in range(len(dna) - self.k + 1):
            kmer = dna[i:i + self.k]
            if kmer in self.counts:
                self.counts[kmer] += 1
            # If kmer contains invalid characters, it's ignored

        return self.counts

    def ffp_normalize_counts(self, dna_length: int) -> Dict[str, float]:
        """
        Normalize k-mer counts by DNA length - (k - 1).

        Parameters:
            dna_length (int): The length of the DNA sequence.

        Returns:
            Dict[str, float]: A dictionary mapping each k-mer to its normalized frequency.
        """
        if not self.counts:
            raise ValueError("Counts have not been calculated. Please run count_kmers() first.")

        total_positions = dna_length - self.k + 1
        if total_positions <= 0:
            raise ValueError("DNA length is shorter than k-mer length.")

        self.normalized_counts = {kmer: count / total_positions for kmer, count in self.counts.items()}
        return self.normalized_counts

    def analyze(self, dna: str) -> Dict[str, float]:
        """
        Perform the full k-mer analysis: count k-mers and normalize their counts.

        Parameters:
            dna (str): The DNA sequence to analyze.

        Returns:
            Dict[str, float]: A dictionary mapping each k-mer to its normalized frequency.
        """
        self.count_kmers(dna)
        self.ffp_normalize_counts(len(dna))
        return self.normalized_counts

    def get_kmer_counts(self) -> Dict[str, int]:
        """
        Retrieve the k-mer counts.

        Returns:
            Dict[str, int]: A dictionary mapping each k-mer to its count.
        """
        if not self.counts:
            raise ValueError("Counts have not been calculated. Please run count_kmers() first.")
        return self.counts

    def get_normalized_counts(self) -> Dict[str, float]:
        """
        Retrieve the normalized k-mer counts.

        Returns:
            Dict[str, float]: A dictionary mapping each k-mer to its normalized frequency.
        """
        if not self.normalized_counts:
            raise ValueError("Normalized counts have not been calculated. Please run ffp_normalize_counts() first.")
        return self.normalized_counts

    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert the k-mer counts and normalized counts to a Pandas DataFrame.

        Returns:
            pd.DataFrame: A DataFrame containing k-mers, their counts, and normalized frequencies.
        """
        if not self.counts or not self.normalized_counts:
            raise ValueError("Counts and normalized counts have not been calculated. Please run analyze() first.")

        data = {
            'k-mer': self.all_kmers,
            'Count': [self.counts[kmer] for kmer in self.all_kmers],
            'Normalized Frequency': [self.normalized_counts[kmer] for kmer in self.all_kmers]
        }
        return pd.DataFrame(data)

    def display_top_kmers(self, top_n: int = 10) -> pd.DataFrame:
        """
        Display the top N k-mers based on their counts.

        Parameters:
            top_n (int): The number of top k-mers to display.

        Returns:
            pd.DataFrame: A DataFrame containing the top N k-mers.
        """
        df = self.to_dataframe()
        top_kmers = df.sort_values(by='Count', ascending=False).head(top_n)
        return top_kmers.reset_index(drop=True)

    def export_to_csv(self, filepath: str):
        """
        Export the k-mer analysis results to a CSV file.

        Parameters:
            filepath (str): The path to the output CSV file.
        """
        df = self.to_dataframe()
        df.to_csv(filepath, index=False)

    def plot_kmer_frequencies(self, top_n: int = 20):
        """
        Plot the frequencies of the top N k-mers using a bar chart.

        Parameters:
            top_n (int): The number of top k-mers to plot.
        """
        import matplotlib.pyplot as plt

        top_kmers = self.display_top_kmers(top_n)
        plt.figure(figsize=(10, 6))
        plt.bar(top_kmers['k-mer'], top_kmers['Normalized Frequency'], color='skyblue')
        plt.xlabel('k-mer')
        plt.ylabel('Normalized Frequency')
        plt.title(f'Top {top_n} k-mer Frequencies')
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.show()








############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################







from itertools import product
from typing import Dict

class MarkovKString_encoder:
    """
    A class to handle k-mer generation, counting, normalization, and feature computation
    based on a Markov model for DNA sequences.
    """
    
    def __init__(self, k: int):
        """
        Initialize the MarkovKString with a specific k-mer length.
        
        Parameters:
        - k (int): The length of k-mers to be analyzed.
        """
        if k <= 0:
            raise ValueError("k must be a positive integer.")
        self.k = k
        self.kmers = self.generate_kmers(k)
        self.k1mers = self.generate_kmers(k - 1) if k > 1 else []
        self.k2mers = self.generate_kmers(k - 2) if k > 2 else []
    
    @staticmethod
    def generate_kmers(k: int) -> list:
        """
        Generate all possible k-mers of length k using DNA letters.
        
        Parameters:
        - k (int): Length of the k-mers to generate.
        
        Returns:
        - List[str]: A list of all possible k-mers.
        """
        if k <= 0:
            return []
        return [''.join(p) for p in product('ACGT', repeat=k)]
    
    def count_kmers(self, dna: str, k: int = None) -> Dict[str, int]:
        """
        Count occurrences of each k-mer in a given DNA sequence, assigning zero to those not appearing.
        
        Parameters:
        - dna (str): The DNA sequence to analyze.
        - k (int, optional): Length of k-mers to count. If None, uses the instance's k.
        
        Returns:
        - Dict[str, int]: A dictionary mapping each possible k-mer to its count in the DNA sequence.
        """
        if k is None:
            k = self.k
        if k <= 0:
            raise ValueError("k must be a positive integer.")
        dna_length = len(dna)
        if dna_length < k:
            raise ValueError("DNA sequence is shorter than k.")
        
        # Generate all possible k-mers with initial count 0
        kmers = self.generate_kmers(k)
        counts = {kmer: 0 for kmer in kmers}
        
        # Count k-mers in the DNA sequence
        for i in range(dna_length - k + 1):
            kmer = dna[i:i + k]
            if kmer in counts:
                counts[kmer] += 1
            else:
                # If the k-mer contains invalid characters, you might want to handle it
                pass  # Ignoring invalid k-mers
        
        return counts
    
    def markov_normalize_counts(self, counts: Dict[str, int], dna_length: int, k: int) -> Dict[str, float]:
        """
        Normalize k-mer counts by (DNA length - (k - 1)).
        
        Parameters:
        - counts (Dict[str, int]): Dictionary of k-mer counts.
        - dna_length (int): Length of the DNA sequence.
        - k (int): Length of the k-mers.
        
        Returns:
        - Dict[str, float]: Dictionary of normalized k-mer probabilities.
        """
        if k <= 0:
            raise ValueError("k must be a positive integer.")
        total_positions = dna_length - k + 1
        if total_positions <= 0:
            raise ValueError("Invalid total positions for normalization.")
        
        return {kmer: count / total_positions for kmer, count in counts.items()}
    
    def compute_feature_array(self, dna: str) -> Dict[str, float]:
        """
        Compute feature array for the instance's k using p(x)/p0(x) - 1 logic.
        
        Parameters:
        - dna (str): The DNA sequence to analyze.
        
        Returns:
        - Dict[str, float]: Feature array mapping each k-mer to its feature value.
        """
        k = self.k
        dna_length = len(dna)
        
        if dna_length < k:
            raise ValueError("DNA sequence is shorter than k.")
        if k < 1:
            raise ValueError("k must be at least 1.")
        
        # Step 1: Compute k-mers, (k-1)-mers, and (k-2)-mers probabilities
        kmer_counts = self.count_kmers(dna, k)
        kmer_probs = self.markov_normalize_counts(kmer_counts, dna_length, k)
    
        if k - 1 > 0:
            k1mer_counts = self.count_kmers(dna, k - 1)
            k1mer_probs = self.markov_normalize_counts(k1mer_counts, dna_length, k - 1)
        else:
            k1mer_probs = {}
    
        if k - 2 > 0:
            k2mer_counts = self.count_kmers(dna, k - 2)
            k2mer_probs = self.markov_normalize_counts(k2mer_counts, dna_length, k - 2)
        else:
            k2mer_probs = {}
    
        # Step 2: Compute p0(x) and feature array
        feature_array = {}
        for kmer, prob in kmer_probs.items():
            if k < 2:
                # When k < 2, prefix and suffix do not exist
                p0_x = 0
            else:
                prefix = kmer[:-1]  # (k-1)-mer prefix
                suffix = kmer[1:]   # (k-1)-mer suffix
                middle = kmer[1:-1] if k > 2 else ''
    
                # Compute p0(x), assigning zero if middle probability is zero
                if k > 2:
                    middle_prob = k2mer_probs.get(middle, 0)
                else:
                    middle_prob = 1  # When k=2, there is no (k-2)-mer, so set to 1 to avoid division by zero
    
                if k > 2 and middle_prob > 0:
                    p0_x = (k1mer_probs.get(prefix, 0) * k1mer_probs.get(suffix, 0)) / middle_prob
                elif k == 2:
                    p0_x = k1mer_probs.get(prefix, 0) * k1mer_probs.get(suffix, 0)
                else:
                    p0_x = 0  # When k=1, p0(x) is not defined
    
            # Avoid division by zero
            if p0_x > 0:
                feature_value = prob / p0_x - 1
            else:
                feature_value = 0
            
            feature_array[kmer] = feature_value
    
        return feature_array




############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################






# js and kl


import math
from itertools import product
class js_kl:
    def __init__(self, k):
        self.k = k

    def generate_kmers(self, k):
        """Generate all possible k-mers of length k using DNA letters."""
        return [''.join(p) for p in product('ACGT', repeat=k)]

    def count_kmers(self, dna, k):
        """Count occurrences of each k-mer in a given DNA sequence, assigning zero to those not appearing."""
        all_kmers = {''.join(p): 0 for p in product('ACGT', repeat=k)}
        for i in range(len(dna) - k + 1):
            kmer = dna[i:i + k]
            if kmer in all_kmers:
                all_kmers[kmer] += 1
            else:
                # Handle unexpected characters by ignoring or choose to raise an error
                pass
        return all_kmers


    def js_kl_normalize_counts(self, counts, total_positions):
        """Normalize k-mer counts by the total number of positions."""
        return {kmer: count / total_positions for kmer, count in counts.items()}


    def compute_kl_divergence(self, P, Q, base=2):
        """
        Compute the Kullback-Leibler divergence KL(P || Q).

        Parameters:
        - P (dict): The first probability distribution.
        - Q (dict): The second probability distribution.
        - base (int): The logarithm base to use. Default is 2.

        Returns:
        - float: The KL divergence value. Returns math.inf if divergence is infinite.
        """
        kl_div = 0.0
        for kmer in P:
            p = P[kmer]
            q = Q.get(kmer, 0)
            if p == 0:
                continue  # 0 * log(0 / q) is defined as 0
            if q == 0:
                return math.inf  # KL divergence is infinite
            kl_div += p * math.log(p / q, base)
        return kl_div

    def compute_js_divergence(self, P, Q, base=2):
        """
        Compute the Jensen-Shannon divergence JS(P || Q).

        Parameters:
        - P (dict): The first probability distribution.
        - Q (dict): The second probability distribution.
        - base (int): The logarithm base to use. Default is 2.

        Returns:
        - float: The JS divergence value.
        """
        # Compute the mixture distribution M
        M = {}
        for kmer in P:
            M[kmer] = 0.5 * (P[kmer] + Q.get(kmer, 0))
        
        # Compute KL(P || M) and KL(Q || M)
        kl_p_m = compute_kl_divergence(P, M, base)
        kl_q_m = compute_kl_divergence(Q, M, base)
        
        if kl_p_m == math.inf or kl_q_m == math.inf:
            return math.inf  # Divergence is infinite
        
        # Compute JS divergence
        js_div = 0.5 * kl_p_m + 0.5 * kl_q_m
        return js_div

    def featurize_ffp(self, dna):
        counts = self.count_kmers(dna, self.k)
        total_positions = len(dna) - self.k + 1
        return self.js_kl_normalize_counts(counts, total_positions)









############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################










import numpy as np

class FPS_encoder:
    """
    A class to extract a 12-dimensional feature vector from a DNA sequence
    using one-hot encoding, Discrete Fourier Transform (DFT), power spectrum,
    and moment calculations.
    """
    
    def __init__(self):
        """
        Initialize the DNAFeatureExtractor without a DNA sequence.
        """
        self.dna_sequence = None
        self.mapping = {
    'A': 0,  # Adenine
    'C': 1,  # Cytosine
    'G': 2,  # Guanine
    'T': 3,  # Thymine
    'U': 3,  # Uracil, treat same as T

    # Ambiguity codes (deterministic assumptions)
    'R': 0,  # A or G → assume A
    'Y': 1,  # C or T → assume C
    'S': 2,  # G or C → assume G
    'W': 0,  # A or T → assume A
    'K': 2,  # G or T → assume G
    'M': 0,  # A or C → assume A
    'B': 1,  # C or G or T → assume C
    'D': 0,  # A or G or T → assume A
    'H': 0,  # A or C or T → assume A
    'V': 0,  # A or C or G → assume A
    'N': 0,  # any base → assume A

    # Lowercase variants (optional, for robustness)
    'a': 0, 'c': 1, 'g': 2, 't': 3, 'u': 3,
    'r': 0, 'y': 1, 's': 2, 'w': 0,
    'k': 2, 'm': 0, 'b': 1, 'd': 0,
    'h': 0, 'v': 0, 'n': 0
}

        self.one_hot_matrix = None
        self.dft_matrix = None
        self.power_spectrum = None
        self.feature_vector = None

    def featurize_dna(self, dna_sequence: str):
        """
        Sets the DNA sequence and computes all necessary features.
        
        Parameters:
            dna_sequence (str): The DNA sequence to process.
        
        Raises:
            ValueError: If the DNA sequence contains invalid nucleotides.
        """
        self.dna_sequence = dna_sequence.upper()
        self._validate_dna_sequence()
        self.one_hot_matrix = self.one_hot_encode()
        self.dft_matrix = self.compute_dft()
        self.power_spectrum = self.compute_power_spectrum()
        # self.feature_vector = self.compute_moments()
        return self.compute_moments()
    
    def _validate_dna_sequence(self):
        """
        Validates that the DNA sequence contains only valid nucleotides (A, C, G, T).
        
        Raises:
            ValueError: If an invalid nucleotide is found.
        """
        if self.dna_sequence is None:
            raise ValueError("DNA sequence has not been set.")
        
        for idx, nucleotide in enumerate(self.dna_sequence):
            if nucleotide not in self.mapping:
                raise ValueError(
                    f"Invalid nucleotide '{nucleotide}' found at position {idx}. "
                    "Only A, C, G, T are allowed."
                )
    
    def one_hot_encode(self) -> np.ndarray:
        """
        Converts the DNA sequence into a one-hot encoded 4xN matrix.
        
        Rows correspond to A, C, G, T respectively.
        
        Returns:
            np.ndarray: A 4xN binary matrix.
        """
        N = len(self.dna_sequence)
        one_hot = np.zeros((4, N), dtype=int)
        
        for idx, nucleotide in enumerate(self.dna_sequence):
            row = self.mapping[nucleotide]
            one_hot[row, idx] = 1
        
        return one_hot
    
    def compute_dft(self) -> np.ndarray:
        """
        Computes the Discrete Fourier Transform (DFT) for each row of the one-hot matrix.
        
        Returns:
            np.ndarray: A 4xN matrix of complex numbers representing the DFT.
        """
        return np.fft.fft(self.one_hot_matrix, axis=1)
    
    def compute_power_spectrum(self) -> np.ndarray:
        """
        Computes the power spectrum (magnitude) of the DFT for each row,
        considering up to N/2 frequencies (excluding the zero frequency).
        
        Returns:
            np.ndarray: A 4x(M) matrix of power spectra, where M = N//2.
        """
        N = self.dft_matrix.shape[1]
        # Take frequencies from s=1 to s=N//2 (excluding the zero frequency)
        power_spectrum = np.abs(self.dft_matrix[:, 1:(N//2)+1])
        return power_spectrum
    
    def compute_moments(self) -> list:
        """
        Computes the first three moments for each nucleotide based on the power spectrum.
        
        The moments are calculated as follows:
            M1 = Sum(p^1) * a1
            M2 = Sum(p^2) * a2
            M3 = Sum(p^3) * a3
        where:
            a1 = 1 / (k * (N - k))^0 = 1
            a2 = 1 / (k * (N - k))^1
            a3 = 1 / (k * (N - k))^2
        
        Returns:
            list: A list containing 12 moments [M1_A, M2_A, M3_A, ..., M3_T].
        """
        moments = []
        M = self.power_spectrum.shape[1]
        N = self.one_hot_matrix.shape[1]
        
        for row_idx in range(self.one_hot_matrix.shape[0]):
            k = np.sum(self.one_hot_matrix[row_idx])
            denominator = k * (N - k)
            if denominator == 0:
                # Avoid division by zero; set moments to zero
                a1, a2, a3 = 0, 0, 0
            else:
                a1 = 1 / (denominator)**0  # j=1: a1 = 1 / denominator^0 = 1
                a2 = 1 / (denominator)**1  # j=2: a2 = 1 / denominator^1
                a3 = 1 / (denominator)**2  # j=3: a3 = 1 / denominator^2
            
            # Sum over s=1 to N/2
            sum_p1 = np.sum(self.power_spectrum[row_idx] ** 1)
            sum_p2 = np.sum(self.power_spectrum[row_idx] ** 2)
            sum_p3 = np.sum(self.power_spectrum[row_idx] ** 3)
            
            M1 = a1 * sum_p1
            M2 = a2 * sum_p2
            M3 = a3 * sum_p3
            
            moments.extend([M1, M2, M3])
        
        return moments
    
    def get_feature_vector(self) -> list:
        """
        Retrieves the 12-dimensional feature vector.
        
        Returns:
            list: A list of 12 features [M1_A, M2_A, M3_A, ..., M3_T].
        
        Raises:
            ValueError: If the DNA sequence has not been set.
        """
        if self.feature_vector is None:
            raise ValueError("Feature vector has not been computed. Please set a DNA sequence first.")
        return self.feature_vector
    
    def __str__(self):
        """
        Returns a string representation of the feature vector with nucleotide labels.
        
        Returns:
            str: Formatted string of features per nucleotide.
        
        Raises:
            ValueError: If the DNA sequence has not been set.
        """
        if self.feature_vector is None:
            return "Feature vector has not been computed. Please set a DNA sequence first."
        
        nucleotide_labels = ['A', 'C', 'G', 'T']
        feature_strings = []
        for i, nucleotide in enumerate(nucleotide_labels):
            M1, M2, M3 = self.feature_vector[3*i:3*(i+1)]
            feature_strings.append(
                f"Features for {nucleotide}: M1 = {M1}, M2 = {M2}, M3 = {M3}"
            )
        return "\n".join(feature_strings)






############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################
############################################################################################################################################################################################################





class NVM_encoder:
    def __init__(self, K):
        """
        Initializes the Kmer class with a specific k-mer size.

        Parameters:
        K (int): The size of the k-mer.
        """
        self.K = K
        self.index_map = {
            'a': 0, 'A': 0,
            'c': 1, 'C': 1,
            'g': 2, 'G': 2,
            't': 3, 'T': 3
        }

    def compute(self, sequence):
        """
        Computes the k-mer statistics for a given DNA sequence.

        Parameters:
        sequence (str): The DNA sequence to analyze.

        Returns:
        list: A list containing counts, average positions, and variance-like metrics for each k-mer.
        """
        K = self.K
        m = 4 ** K
        na_vect = [0] * (3 * m)
        pos_sum = [0] * m
        squa_sum = [0] * m
        n = len(sequence) - (K - 1)

        # First Pass: Count k-mers and sum their positions
        for i in range(n):
            flag = True
            for l in range(K):
                if sequence[i + l] not in self.index_map:
                    flag = False
                    break
            if not flag:
                continue

            # Calculate the k-mer index
            tem = self.index_map[sequence[i]]
            for l in range(1, K):
                tem = 4 * tem + self.index_map[sequence[i + l]]

            na_vect[tem] += 1
            pos_sum[tem] += i + 1

        # Compute Average Positions
        for k in range(m):
            if na_vect[k] != 0:
                na_vect[k + m] = pos_sum[k] / na_vect[k]
            else:
                na_vect[k + m] = 0

        # Second Pass: Compute Squared Deviations
        for i in range(n):
            flag = True
            for l in range(K):
                if sequence[i + l] not in self.index_map:
                    flag = False
                    break
            if not flag:
                continue

            # Calculate the k-mer index
            tem = self.index_map[sequence[i]]
            for l in range(1, K):
                tem = 4 * tem + self.index_map[sequence[i + l]]

            squa_sum[tem] += (i + 1 - na_vect[tem + m]) ** 2

        # Compute Variance-like Metric
        for k in range(m):
            if na_vect[k] != 0:
                na_vect[k + 2 * m] = squa_sum[k] / (n * na_vect[k])
            else:
                na_vect[k + 2 * m] = 0

        return na_vect

    def __repr__(self):
        return f"Kmer(K={self.K})"




import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import entropy

def kullback_leibler_divergence(p, q):
    """Compute KL divergence between two probability distributions."""
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

def jensen_shannon_divergence(p, q):
    """Compute the Jensen-Shannon divergence between two distributions."""
    m = 0.5 * (p + q)
    return 0.5 * kullback_leibler_divergence(p, m) + 0.5 * kullback_leibler_divergence(q, m)

# # Example data: rows are probability distributions
# data = np.random.dirichlet(alpha=np.ones(5), size=10)  # 10 rows of probability distributions over 5 categories

# # Compute pairwise KL divergences
# kl_distances = squareform(pdist(data, metric=lambda u, v: kullback_leibler_divergence(u, v)))

# # Compute pairwise Jensen-Shannon divergences
# js_distances = squareform(pdist(data, metric=lambda u, v: jensen_shannon_divergence(u, v)))

# print("KL Distance Matrix:")
# print(kl_distances)

# print("\nJensen-Shannon Distance Matrix:")
# print(js_distances)
