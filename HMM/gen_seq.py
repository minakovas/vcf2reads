from HMM import *
from typing import List
from collections import Counter
import numpy as np

def split2kmers(s: str, k: int) -> List[str]:
    """
    split given string to kmers of given length
    Example:
    split2kmers('QWERTY', 2) -> ['QW', 'WE', 'ER', 'RT', 'TY']
    """
    lst = []
    for i in range(0, len(s) - k + 1):
        kmer = s[i:i + k]
        lst.append(kmer)
    return lst

def assess_p_site(model, seqs):
    p_site_list = []
    for seq in seqs:
        decoding = model.viterbi(seq)
        kmers = split2kmers(decoding, 2)
        kmers_counts = Counter(kmers)
        p = kmers_counts["BS"] / (kmers_counts["BS"] + kmers_counts["BB"])
        p_site_list.append(p)
    return np.mean(p_site_list), np.std(p_site_list)

if __name__ == "__main__":
    print("Generating sequences of given length:\n")
    for p_site in [0.01, 0.1, 0.5]:
        model = create_restriction_model("ACGT", p_site=p_site)
        seqs = []
        # generate 100 sequences of length 1000
        for _ in range(100):
            seqs.append(model.generate_sequence(1000))
        # obtain viterbi deconding and counn p_site from generated sequence
        p_mean, p_std = assess_p_site(model, seqs)
        print(f"Model with p_site: {p_site}, p_site assessed from generated sequences: {p_mean} ± {p_std}")

    print("\nGenerating sequences of random length (p_end=0.01):\n")
    for p_site in [0.01, 0.1, 0.5]:
        model = create_restriction_model("ACGT", p_site=p_site)
        seqs = []
        p_site_list = []
        # generate 100 sequences of random length
        for _ in range(100):
            seqs.append(model.generate_sequence(1000))
        # obtain viterbi deconding and counn p_site from generated sequence
        p_mean, p_std = assess_p_site(model, seqs)
        print(f"Model with p_site: {p_site}, p_site assessed from generated sequences: {p_mean} ± {p_std}")

