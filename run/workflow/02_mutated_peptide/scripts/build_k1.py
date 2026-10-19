#!/usr/bin/env python3
import argparse
import json
import math
import os
import sqlite3
from functools import lru_cache
import numpy as np

# read BLOSUM62-2 and build k1 once

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

