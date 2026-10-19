#!/usr/bin/env python3
import math
import numpy as np

# Levenshtein similarity
def levenshtein(a, b):
    if len(a) < len(b):
        a, b = b, a

    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, start=1):
        curr = [i]
        for j, cb in enumerate(b, start=1):
            ins = curr[j - 1] + 1
            dele = prev[j] + 1
            sub = prev[j - 1] + (ca != cb)
            curr.append(min(ins, dele, sub))
        prev = curr
    return prev[-1]


def edit_similarity(a, b):
    if not a and not b:
        return 1.0
    d = levenshtein(a, b)
    return 1.0 - d / max(len(a), len(b))

