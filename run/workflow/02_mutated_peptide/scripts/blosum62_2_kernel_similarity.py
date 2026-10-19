#!/usr/bin/env python3
import math
import os
from functools import lru_cache

import numpy as np


def read_blosum62_2(blos_path):
    """
    Reads lower-triangular BLOSUM62-2 file.
    """
    line_list = []
    amino_acids = None

    with open(blos_path, "r", encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            elif s.startswith("A"):
                amino_acids = "".join(s.split())
                continue
            else:
                vals = [float(x) for x in s.split()]
                line_list.append(vals)

    if amino_acids is None:
        raise ValueError("Amino acid header line not found in BLOSUM62-2 file.")

    d = len(amino_acids)
    if len(line_list) != d:
        raise ValueError(f"Expected {d} matrix rows, got {len(line_list)}.")

    mat = np.zeros((d, d), dtype=float)
    for i in range(d):
        if len(line_list[i]) != i + 1:
            raise ValueError(f"Row {i} should have {i + 1} values, got {len(line_list[i])}.")
        for j in range(d):
            if j <= i:
                mat[i, j] = line_list[i][j]
            else:
                mat[i, j] = line_list[j][i]

    return mat, amino_acids


def build_or_load_k1(blos_path, beta=0.125, cache_path=None):
    """
    Build K1 once and persist it on disk.
    """
    if cache_path is None:
        cache_path = f"{blos_path}.beta_{beta}.k1_cache.npz"

    if os.path.exists(cache_path):
        data = np.load(cache_path, allow_pickle=False)
        k1 = data["k1"]
        amino_acids = str(data["amino_acids"])
        return k1, amino_acids

    blos_mat, amino_acids = read_blosum62_2(blos_path)
    k1 = np.power(blos_mat, beta)
    np.savez_compressed(cache_path, k1=k1, amino_acids=amino_acids)
    return k1, amino_acids


class KernelScorer:
    """
    K2, K3_sum, K3_hat
    """

    def __init__(self, blos_path, beta=0.125, cache_path=None):
        self.k1, self.amino_acids = build_or_load_k1(
            blos_path=blos_path,
            beta=beta,
            cache_path=cache_path,
        )
        self.aa_to_idx = {aa: i for i, aa in enumerate(self.amino_acids)}

    def _has_invalid_aa(self, seq):
        bad = [c for c in seq if c not in self.aa_to_idx]
        return sorted(set(bad))
            
            
            

    @lru_cache(maxsize=500000)
    def k2_prod(self, seq1, seq2, start1, start2, k):
        """
        K2 on one substring pair:
          K2 = product_i K1
        """
        prod = 1.0
        for i in range(k):
            a1 = seq1[start1 + i]
            a2 = seq2[start2 + i]
            i1 = self.aa_to_idx[a1]
            i2 = self.aa_to_idx[a2]
            prod *= self.k1[i1, i2]
        return prod

    def k3_sum(self, seq1, seq2, p_min=1, p_max=25):
        """
        K3(f,g) = sum over all substring pairs of lengths k in [p_min, p_max]
        """

        l1 = len(seq1)
        l2 = len(seq2)
        if l1 == 0 or l2 == 0:
            return 0.0

        k_upper = min(p_max, l1, l2)
        if k_upper < p_min:
            return 0.0

        total = 0.0
        for k in range(p_min, k_upper + 1):
            for start1 in range(l1 - k + 1):
                for start2 in range(l2 - k + 1):
                    total += self.k2_prod(seq1, seq2, start1, start2, k)
        return total

    def k3_hat(self, seq1, seq2, p_min=1, p_max=25):
        """
        Normalized correlation kernel:
          Khat = K3(seq1,seq2) / sqrt(K3(seq1,seq1) * K3(seq2,seq2))
        """
        bad1 = self._has_invalid_aa(seq1)
        bad2 = self._has_invalid_aa(seq2)

        if bad1 or bad2:
            return 0.0
        
        score11 = self.k3_sum(seq1, seq1, p_min=p_min, p_max=p_max)
        score22 = self.k3_sum(seq2, seq2, p_min=p_min, p_max=p_max)
        score12 = self.k3_sum(seq1, seq2, p_min=p_min, p_max=p_max)

        if score11 <= 0 or score22 <= 0:
            return 0.0

        return score12 / math.sqrt(score11 * score22)
