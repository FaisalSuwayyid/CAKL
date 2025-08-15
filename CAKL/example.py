# featurize_one_dna.py
# Minimal example: compute PH-based features for ONE artificial DNA sequence.

import os
import numpy as np
from itertools import product

from psrt import PH, generate_all_kmers, occurrence

# -------------------------------
# Configuration (edit as you like)
# -------------------------------
encoder        = "PSRT"          # used for output folder naming only
data_name      = "ExampleSet"    # used for output folder naming only
accession_id   = "EXAMPLE_0001"  # filename label
alphabet       = ["A", "C", "G", "T"]

# Artificial DNA (make it decently long to see structure; tweak freely)
dna = (
    "ACGT"*50 +               # 200 bp of periodic pattern
    "AAAACCCCGGGGTTTT"*5 +    # 80 bp block pattern
    "AGTCAGTCAGTT"*20         # 240+ bp quasi-repeat
)
# Choose a single k and PH settings
k_chosen      = 3          # k-mer size
max_dimension = 2          # compute up to H_2 features
num           = 2         # number of filtration steps (creates num+1 grid points)

# Filtration grid (uniform, scaled by 4**k like in your code)
specific_filtration = np.array([i * (4 ** k_chosen) for i in range(0, num + 1)])

# Output root
out_root = f"./features/{encoder}/{data_name}/{k_chosen}"
os.makedirs(out_root, exist_ok=True)

# -------------------------------
# Prepare result collectors
# -------------------------------
betti_result = (max_dimension + 1) * [None]
f_result     = (max_dimension + 1) * [None]
h_result     = (max_dimension + 2) * [None]
facet_result = (max_dimension + 1) * [None]

# -------------------------------
# Generate all k-mers and process
# -------------------------------
kmers = generate_all_kmers(alphabet, k_chosen)

for kmer in kmers:
    pts = occurrence(dna, kmer)

    # If no occurrences, append zero curves of correct length
    if len(pts) == 0:
        zero_curve = np.zeros(num + 1, dtype=float)
        # For d = 0..max_dimension
        for d in range(max_dimension + 1):
            betti_result[d] = zero_curve if betti_result[d] is None else np.vstack([betti_result[d], zero_curve])
            f_result[d]     = zero_curve if f_result[d]     is None else np.vstack([f_result[d],     zero_curve])
            facet_result[d] = zero_curve if facet_result[d] is None else np.vstack([facet_result[d], zero_curve])
        # For d = 0..max_dimension+1 (h-vector goes one higher)
        for d in range(max_dimension + 2):
            h_result[d] = zero_curve if h_result[d] is None else np.vstack([h_result[d], zero_curve])
        continue

    # Prepare point cloud (positions as 1D points in R^1)
    points = np.array(pts)[:, np.newaxis]

    # Build PH object
    ph = PH(
        points,
        max_dimension=max_dimension,
        max_edge_length=2.0,                 # keep as in your code; adjust if needed
        specific_filtration=specific_filtration
    )

    # Compute curves
    alphas, betti_curves = ph.betti_curves()
    f_curves     = ph.compute_f_vector_curves()
    h_curves     = ph.compute_h_vector_curves()
    facet_curves = ph.facet_curves()

    # Stack per dimension, filling missing with zeros_like(alphas)
    z = np.zeros_like(alphas)

    for d in range(max_dimension + 1):
        b = betti_curves.get(d, z)
        f = f_curves.get(d, z)
        fc = facet_curves.get(d, z)

        betti_result[d] = b if betti_result[d] is None else np.vstack([betti_result[d], b])
        f_result[d]     = f if f_result[d]     is None else np.vstack([f_result[d],     f])
        facet_result[d] = fc if facet_result[d] is None else np.vstack([facet_result[d], fc])

    for d in range(max_dimension + 2):
        h = h_curves.get(d, z)
        h_result[d] = h if h_result[d] is None else np.vstack([h_result[d], h])

# -------------------------------
# Save .npy outputs (one DNA)
# -------------------------------
base = os.path.join(out_root, accession_id)
for d in range(max_dimension + 1):
    np.save(f"{base}_betti{d}.npy", betti_result[d])
    np.save(f"{base}_f{d}.npy",      f_result[d])
    np.save(f"{base}_facet{d}.npy",  facet_result[d])

for d in range(max_dimension + 2):
    np.save(f"{base}_h{d}.npy", h_result[d])

print("Done.")
print("Saved files like:", f"{base}_betti0.npy", f"{base}_f0.npy", f"{base}_h0.npy", f"{base}_facet0.npy")
