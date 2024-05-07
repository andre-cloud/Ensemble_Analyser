import pytest
import numpy as np
from src.rrho import calc_zpe, free_gibbs_energy


def test_calc_zpe():
    # Test case 1: frequency array is empty
    freq = np.array([])
    expected_zpe = 0
    assert np.isclose(calc_zpe(freq), expected_zpe)

    # Test case 2: frequency array has one value
    freq = np.array([1000])
    expected_zpe = 0.0022781676264559806
    assert np.isclose(calc_zpe(freq), expected_zpe)

    # Test case 3: frequency array has multiple values
    freq = np.array([1000, 2000, 3000])
    expected_zpe = 0.013669005758735883
    assert np.isclose(calc_zpe(freq), expected_zpe)


def test_free_energy():
    SCF = 10.0
    T = 298.15
    freq = np.array([100, 200, 300])
    mw = 28.97
    B = np.array([1.0, 2.0, 3.0])
    m = 1
    linear = False
    cut_off = 100
    alpha = 4
    P = 101.325

    result = free_gibbs_energy(SCF, T, freq, mw, B, m, linear, cut_off, alpha, P)

    # Perform your assertions on the result with
    assert np.isclose(result, 9.977753201998967)

    linear = True
    result = free_gibbs_energy(SCF, T, freq, mw, B, m, linear, cut_off, alpha, P)

    # Perform your assertions on the result with
    assert np.isclose(result, 9.979965792434896)
