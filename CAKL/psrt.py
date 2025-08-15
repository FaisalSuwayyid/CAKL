import numpy as np
import gudhi as gd
import math
import bisect

""" The source code for the methods used for comparison in this study can be found in https://github.com/hozumiyu/KmerTopology, thanks to Dr. Yuta Hozumi """


class PH:
    """
    Persistent homology on a point cloud via Vietoris–Rips,
    with caching of SimplexTree builds and persistence computations.
    """



    def __init__(self, points, max_edge_length=1.0, max_dimension=1,
                min_edge_length=0.0, num_instances=200, specific_filtration=None):
        self.points = points
        self.max_edge_length = max_edge_length
        self.max_dimension = max_dimension
        self.min_edge_length = min_edge_length
        self.num_instances = num_instances
        self.specific_filtration = np.array(specific_filtration) if specific_filtration is not None else specific_filtration

        self._st_current = {"dim": -1,
                            "tree": None # dimension and tree
                            }
        self._persisted_st_current = {"dim": -1,
                            "tree": None # dimension and tree
                            }



    # persistent homology


    def rips(self, max_edge_length=None):
        if max_edge_length is None:
            max_edge_length = self.max_edge_length
        self.max_edge_length = max_edge_length
        return gd.RipsComplex(points=self.points, max_edge_length=max_edge_length)

    def simplex_tree(self, max_dimension=None):
        """
        Build (or reuse) the SimplexTree up to `max_dimension`.
        """
        if max_dimension is None:
            max_dimension = self.max_dimension
        self.max_dimension = max_dimension

        if max_dimension != self._st_current["dim"]:
            st = self.rips().create_simplex_tree(max_dimension=max_dimension)
            self._st_current["dim"] = max_dimension
            self._st_current["tree"] = st
        return self._st_current["tree"]

    def persisted_simplex_tree(self, max_dimension=None):
        """
        Compute persistence on the simplex tree (or reuse if already done).
        """
        if max_dimension is None:
            max_dimension = self.max_dimension
        self.max_dimension = max_dimension

        if max_dimension != self._persisted_st_current["dim"]:
            st = self.simplex_tree(max_dimension)
            st.compute_persistence()
            self._persisted_st_current["tree"] = st
            self._persisted_st_current["dim"] = max_dimension
            self._st_current["tree"] = st
            self._st_current["dim"] = max_dimension
        return self._persisted_st_current["tree"]

    def simplices_filtration(self, max_dimension=None):
        st = self.simplex_tree(max_dimension)
        return list(st.get_simplices()) # list(st.get_filtration())

    def simplices(self, max_dimension=None):
        return [s for s, _ in self.simplices_filtration(max_dimension)]

    def simplices_per_dimension(self, dim=0):
        return [
            s
            for s, _ in self.simplices_filtration(max_dimension=dim)
            if len(s) == dim + 1
        ]

    def simplices_all_per_dimension(self, max_dimension=None):
        if max_dimension is None:
            max_dimension = self.max_dimension
        self.max_dimension = max_dimension

        groups = [[] for _ in range(max_dimension + 1)]
        for s, _ in self.simplices_filtration(max_dimension):
            groups[len(s) - 1].append(s)
        return groups

    def persistence_intervals_by_dimension(self, max_dimension=None):
        st = self.persisted_simplex_tree(max_dimension)
        return {
            dim: st.persistence_intervals_in_dimension(dim)
            for dim in range((max_dimension or self.max_dimension) + 1)
        }

    def betti_at_filtration(self, dim=0, filtration=0.0, max_dimension=None):
        intervals = self.persistence_intervals_by_dimension(max_dimension)[dim]
        return sum(1 for b, d in intervals if b <= filtration < d)

    def betti_curves(self, max_dimension=None, r1=None, r2=None, steps=None, specific_filtration=None):
        if max_dimension is None:
            max_dimension = self.max_dimension
        self.max_dimension = max_dimension

        if self.specific_filtration is None:
            if r1 is None:
                r1 = self.min_edge_length
            self.min_edge_length = r1
            # determine r2 from constructed filtration if not given
            if r2 is None:
                r2 = self.max_edge_length
            self.max_edge_length = r2
            
            if steps is None:
                steps = self.num_instances
            self.num_instances = steps
            self.specific_filtration = np.linspace(r1, r2, steps)
        specific_filtration = self.specific_filtration

        intervals_by_dim = self.persistence_intervals_by_dimension(max_dimension)
        alphas = np.array(specific_filtration)

        curves = {}
        for d in range(max_dimension + 1):
            ivs = intervals_by_dim.get(d, [])
            curves[d] = np.array(
                [sum(1 for b, death in ivs if b <= α < death) for α in alphas],
                dtype=int
            )

        return alphas, curves
    


    # f- and h- vectors




    def compute_f_vector_at_t(self, simplices_with_filtration, t):
        """
        Given a list of (simplex, filtration) pairs, compute the f-vector for the
        subcomplex at time t (i.e. include all simplices with filtration value <= t).
        """
        f_vector = [0] * (self.max_dimension + 1)
        for simplex, filtration in simplices_with_filtration:
            if filtration <= t:
                dim = len(simplex) - 1
                if dim <= self.max_dimension:
                    f_vector[dim] += 1
        return f_vector

    def compute_f_vector_curves(self):
        """
        Compute the f-vector curves for a range of filtration values.

        Returns:
            t_values (np.ndarray): Filtration thresholds.
            f_curves (dict): Dictionary where f_curves[d] is a list of the number of 
                            d-dimensional simplices for each t in t_values.
        """

        if self._st_current["tree"] is None:
            simplex_tree = self.simplex_tree()
        # rips_complex = gd.RipsComplex(points=self.points, max_edge_length=self.max_edge_length)
        # simplex_tree = rips_complex.create_simplex_tree(max_dimension=self.max_dimension + 1)
        # simplices_with_filtration = list(simplex_tree.get_simplices())

        simplices_with_filtration = self.simplices_filtration()

        if self.specific_filtration is None:
            self.specific_filtration = np.linspace(self.min_edge_length, self.max_edge_length, self.num_instances)

        t_values = self.specific_filtration # np.linspace(self.min_edge_length, self.max_edge_length, self.num_instances)
        f_curves = {d: [] for d in range(self.max_dimension + 1)}
        for t in t_values:
            fv = self.compute_f_vector_at_t(simplices_with_filtration, t)
            for d in range(self.max_dimension + 1):
                f_curves[d].append(fv[d])
        return f_curves


    def compute_rate_curves(self):
        """
        Compute rate curves for each dimension using:
        Cumulative rate: f(t) / t.
        """
        f_curves = self.compute_f_vector_curves()
    

        if self.specific_filtration is None:
            self.specific_filtration = np.linspace(self.min_edge_length, self.max_edge_length, self.num_instances)

        t_values = self.specific_filtration

        cumulative_rates = {}
        for d, f_vals in f_curves.items():
            f_vals = np.array(f_vals)
            cum_rate = np.zeros_like(f_vals, dtype=float)
            nonzero = t_values > 0
            cum_rate[nonzero] = f_vals[nonzero] / t_values[nonzero]
            cumulative_rates[d] = cum_rate
        return cumulative_rates

    def compute_h_vector(self, f_counts):
        """
        Computes the h-vector from the f-vector.
        """
        last_nonzero = max(i for i, val in enumerate(f_counts) if val != 0)
        d = last_nonzero + 1

        h = (2 + self.max_dimension) * [0]
        for k in range(d + 1):
            h_k = 0
            for i in range(k + 1):
                f_val = 1 if i == 0 else f_counts[i - 1]
                h_k += (-1)**(k - i) * math.comb(d - i, k - i) * f_val
            h[k] = h_k
        return h

    def compute_h_vector_curves(self):
        """
        Apply compute_h_vector to each time slice of the f-vector curves.
        """
        f_curves = self.compute_f_vector_curves()
        d = self.max_dimension + 1
        num = self.specific_filtration.size

        h_curves = {k: [] for k in range(d + 1)}
        for idx in range(num):
            f_at_idx = [f_curves[j][idx] for j in range(self.max_dimension + 1)]
            h_vec = self.compute_h_vector(f_at_idx)
            for k in range(d + 1):
                h_curves[k].append(h_vec[k])
        return h_curves



    # def compute_h_rate_curves(self):
    #     """
    #     Apply compute_h_vector to each time slice of the cumulative rate curves.
    #     """
    #     cum_rates = self.compute_rate_curves()
    #     d = self.max_dimension + 1
    #     num = self.num_instances

    #     h_rate_curves = {k: [] for k in range(d + 1)}
    #     for idx in range(num):
    #         rate_at_idx = [cum_rates[j][idx] for j in range(self.max_dimension + 1)]
    #         h_vec = self.compute_h_vector(rate_at_idx)
    #         for k in range(d + 1):
    #             h_rate_curves[k].append(h_vec[k])
    #     return h_rate_curves




    # Facet betti numbers


    def vietoris_rips_simplices_birth_death(self):
        n = len(self.points)
        # rips_complex = gd.RipsComplex(points=self.points, max_edge_length=self.max_edge_length)
        # simplex_tree = rips_complex.create_simplex_tree(max_dimension=self.max_dimension + 1)
        # simplex_entries = list(simplex_tree.get_simplices())

        if self._st_current["tree"] is None:
            simplex_tree = self.simplex_tree()
        # rips_complex = gd.RipsComplex(points=self.points, max_edge_length=self.max_edge_length)
        # simplex_tree = rips_complex.create_simplex_tree(max_dimension=self.max_dimension + 1)
        # simplices_with_filtration = list(simplex_tree.get_simplices())

        simplex_entries = self.simplices_filtration()

        subset_birth = {}
        for simplex, filt in simplex_entries:
            if (len(simplex) - 1) <= self.max_dimension:
                subset_birth[frozenset(simplex)] = filt

        subset_death = {}
        for subset, birth in subset_birth.items():
            k = len(subset)
            if k == self.max_dimension + 1:
                subset_death[subset] = math.inf
            else:
                candidate_deaths = []
                for v in range(n):
                    if v not in subset:
                        superset = subset.union({v})
                        if len(superset) == k + 1 and superset in subset_birth:
                            candidate_deaths.append(subset_birth[superset])
                subset_death[subset] = min(candidate_deaths) if candidate_deaths else math.inf

        barcodes = {dim: [] for dim in range(self.max_dimension + 1)}
        for subset, birth in subset_birth.items():
            dim = len(subset) - 1
            if dim <= self.max_dimension:
                death = subset_death[subset]
                if birth != death:
                    barcodes[dim].append((birth, death))
        return barcodes

    def prepare_intervals_from_barcodes(self, dimension=None):
        barcodes = self.vietoris_rips_simplices_birth_death()
        intervals = []
        if dimension is not None:
            intervals.extend(barcodes.get(dimension, []))
        else:
            for d in barcodes:
                intervals.extend(barcodes[d])
        births = sorted(b for b, d in intervals)
        deaths = sorted(d for b, d in intervals)
        return births, deaths

    def count_active_intervals_sorted(self, births, deaths, t):
        started = bisect.bisect_right(births, t)
        ended = bisect.bisect_right(deaths, t)
        return started - ended

    def compute_active_curve(self, dimension, t_values):
        births, deaths = self.prepare_intervals_from_barcodes(dimension=dimension)
        counts = [self.count_active_intervals_sorted(births, deaths, t) for t in t_values]
        return counts

    def compute_active_rates(self, dimension, t_values):
        births, deaths = self.prepare_intervals_from_barcodes(dimension=dimension)
        rates = [self.count_active_intervals_sorted(births, deaths, t) / t for t in t_values if t != 0]
        return rates

    def facet_curves(self):
        # t_values = np.linspace(self.min_edge_length, self.max_edge_length, self.num_instances)
        if self.specific_filtration is None:
            self.specific_filtration = np.linspace(self.min_edge_length, self.max_edge_length, self.num_instances)
        t_values = self.specific_filtration
        curves = {dimension: self.compute_active_curve(dimension, t_values) for dimension in range(self.max_dimension + 1)}
        return curves

    def facet_rates(self):
        # t_values = np.linspace(self.min_edge_length, self.max_edge_length, self.num_instances)
        if self.specific_filtration is None:
            self.specific_filtration = np.linspace(self.min_edge_length, self.max_edge_length, self.num_instances)
        t_values = self.specific_filtration
        curves = {dimension: self.compute_active_rates(dimension, t_values) for dimension in range(self.max_dimension + 1)}
        return curves
















from itertools import product
import random

def occurrence(dna, kmer):
    """
    Find all start positions where a k-mer occurs in a DNA sequence.

    Parameters:
        dna (str): The DNA sequence.
        kmer (str): The k-mer to search for.

    Returns:
        list of int: Starting indices where the k-mer is found in the DNA sequence.
    """
    k = len(kmer)
    return [i for i in range(len(dna) - k + 1) if dna[i:i + k] == kmer]


def average_difference(positions):
    """
    Compute the average difference between consecutive positions in a list.

    Parameters:
        positions (list of int): A list of sorted integers.

    Returns:
        float: The average gap between consecutive elements, or 0 if not enough data.
    """
    return np.mean(np.diff(positions)) if len(positions) > 1 else 0.0

 
def generate_all_kmers(alphabet, k):
    """
    Generate all possible k-mers from the given alphabet.

    Parameters:
        alphabet (str or list): The set of symbols, e.g., 'ACGT' or ['A', 'C', 'G', 'T']
        k (int): The length of each k-mer

    Returns:
        List[str]: All k-mers of length k
    """
    return [''.join(p) for p in product(alphabet, repeat=k)]


def average_length(dna_list):
    """
    Compute the average length of DNA sequences in a list.

    Parameters:
        dna_list (list of str): A list where each element is a DNA sequence string.

    Returns:
        float: The average sequence length.
    """
    if not dna_list:
        return 0.0
    return np.mean([len(seq) for seq in dna_list])

def overall_kmer_avg_occurrence(dna_list, kmer):
    """
    Compute the average spacing between consecutive occurrences of a k-mer across all sequences.

    Parameters:
        dna_list (list of str): List of DNA sequences.
        kmer (str): The k-mer to evaluate.

    Returns:
        float: Mean of average spacings for the k-mer in sequences where it occurs more than once.
    """
    avg_diffs = [
        avg for seq in dna_list 
        if (avg := average_difference(occurrence(seq, kmer))) > 0
    ]
    return np.mean(avg_diffs) if avg_diffs else 0.0

def all_kmers_avg_occurrence(dna_list, alphabet, k):
    """
    Compute the mean average spacing across all k-mers generated from the given alphabet.

    Parameters:
        dna_list (list of str): List of DNA sequences.
        alphabet (str or list): Alphabet used to generate k-mers.
        k (int): Length of the k-mers.

    Returns:
        float: Mean of the average spacings for all possible k-mers.
    """
    kmers = generate_all_kmers(alphabet, k)
    all_averages = [overall_kmer_avg_occurrence(dna_list, kmer) for kmer in kmers]
    return np.mean(all_averages)


def generate_random_dna(length, alphabet="ACGT"):
    """
    Generate a random DNA sequence of the given length.

    Parameters:
    - length (int): Length of the DNA sequence.
    - alphabet (str): Allowed nucleotide characters (default: 'ACGT').

    Returns:
    - str: Random DNA sequence.
    """
    return ''.join(random.choices(alphabet, k=length))
