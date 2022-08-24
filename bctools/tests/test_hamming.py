#!/usr/bin/env ipython

import bctools as bc
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

HERE = Path(__file__).parent


def test_map2d():
    x = np.array([1, 2, 3, 4])
    y = np.array([50, 60, 70, 80])

    out = bc.hamming.map2d(x, y, lambda a, b: a * b)

    correct = np.array(
        [
            [50, 60, 70, 80],
            [100, 120, 140, 160],
            [150, 180, 210, 240],
            [200, 240, 280, 320],
        ]
    )
    assert np.array_equal(out, correct)


def test_map2d_sym():
    """Like test_map2d but with the symmetry option True"""
    x = np.array([1, 2, 3, 4])

    out = bc.hamming.map2d(x, x, lambda a, b: a * b, symmetry=True)
    correct = np.array([[1, 2, 3, 4], [2, 4, 6, 8], [3, 6, 9, 12], [4, 8, 12, 16]])
    assert np.array_equal(out, correct)


def test_hamming(monkeypatch):
    x = pd.Series(dict(AAAAAAA=20, AAAABBB=10, BAAAAAB=8, CCCCCCC=10, CCDDCCC=10))

    monkeypatch.setattr(plt, "show", lambda: None)
    keep, reject, tie = bc.hamming.hamming_filter(x, min_distance=4)

    assert (keep == np.array(["AAAAAAA", "CCCCCCC", "CCDDCCC"])).all()
    assert (reject == np.array(["AAAABBB", "BAAAAAB"])).all()
    assert tie == dict(CCCCCCC=["CCDDCCC"], CCDDCCC=["CCCCCCC"])
