#!/usr/bin/env python3
import levenshtein_similarity


def hamming_distance(a, b):
    if len(a) != len(b):
        raise ValueError("Hamming distance requires equal-length strings.")
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def hamming_similarity(a, b):
    """
    Fast similarity for equal-length strings:
      1 - (#mismatches / length)

    If lengths differ, falls back to Levenshtein similarity.
    """
    if len(a) != len(b):
        return levenshtein_similarity.edit_similarity(a, b)

    if len(a) == 0:
        return 1.0

    d = hamming_distance(a, b)
    return 1.0 - d / len(a)
