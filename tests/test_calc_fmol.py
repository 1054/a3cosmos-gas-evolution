import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'a3cosmos_gas_evolution', 'Common_Python_Code'))

import numpy as np
import pytest
from calc_fmol import (
    calc_fmol_from_metalZ_following_Krumholz2009,
    calc_fmol_from_metalZ_following_Dave2016,
    calc_fmol_from_metalZ_following_Popping2014,
)


class TestCalcFmol:

    def test_krumholz2009_scalar(self):
        metalZOH = 8.67
        result = calc_fmol_from_metalZ_following_Krumholz2009(metalZOH)
        assert 0.0 <= result <= 1.0

    def test_krumholz2009_array(self):
        metalZOH = np.array([8.0, 8.5, 8.67, 9.0])
        result = calc_fmol_from_metalZ_following_Krumholz2009(metalZOH)
        assert len(result) == 4
        assert np.all(result >= 0.0) and np.all(result <= 1.0)

    def test_krumholz2009_low_metallicity(self):
        metalZOH = 7.5
        result = calc_fmol_from_metalZ_following_Krumholz2009(metalZOH)
        assert 0.0 <= result <= 1.0

    def test_krumholz2009_high_metallicity(self):
        metalZOH = 9.5
        result = calc_fmol_from_metalZ_following_Krumholz2009(metalZOH)
        assert 0.0 <= result <= 1.0

    def test_dave2016_scalar(self):
        metalZOH = 8.67
        result = calc_fmol_from_metalZ_following_Dave2016(metalZOH)
        assert 0.0 <= result <= 1.0

    def test_dave2016_array(self):
        metalZOH = np.array([8.0, 8.5, 8.67, 9.0])
        result = calc_fmol_from_metalZ_following_Dave2016(metalZOH)
        assert len(result) == 4

    def test_dave2016_low_metallicity(self):
        metalZOH = 7.5
        result = calc_fmol_from_metalZ_following_Dave2016(metalZOH)
        assert result > -2.0

    def test_popping2014_scalar(self):
        metalZOH = 8.67
        U_MW = 1.0
        result = calc_fmol_from_metalZ_following_Popping2014(np.array([metalZOH]), U_MW=np.array([U_MW]))
        assert result[0] >= 0.0 and result[0] <= 1.0

    def test_popping2014_array(self):
        metalZOH = np.array([8.0, 8.5, 8.67, 9.0])
        U_MW = np.ones_like(metalZOH)
        result = calc_fmol_from_metalZ_following_Popping2014(metalZOH, U_MW=U_MW)
        assert len(result) == 4

    def test_fmol_increases_with_metallicity(self):
        metalZ_low = 8.0
        metalZ_high = 9.0
        result_low = calc_fmol_from_metalZ_following_Krumholz2009(metalZ_low)
        result_high = calc_fmol_from_metalZ_following_Krumholz2009(metalZ_high)
        assert result_high > result_low
