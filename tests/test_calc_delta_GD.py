import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'a3cosmos_gas_evolution', 'Common_Python_Code'))

import numpy as np
import pytest
from calc_delta_GD import (
    calc_deltaGD_from_metalZ_following_Leroy2011,
    calc_deltaGD_from_metalZ_following_Magdis2012,
    calc_deltaGD_from_metalZ_following_RemyRuyer2014a,
    calc_deltaGD_from_metalZ_following_RemyRuyer2014b,
    calc_deltaGD_from_metalZ_following_Genzel2015,
)


class TestCalcDeltaGD:

    def test_leroy2011_scalar(self):
        result = calc_deltaGD_from_metalZ_following_Leroy2011(8.67)
        assert result > 0

    def test_leroy2011_array(self):
        metalZ = np.array([8.0, 8.5, 8.67, 9.0])
        result = calc_deltaGD_from_metalZ_following_Leroy2011(metalZ)
        assert len(result) == 4
        assert np.all(result > 0)

    def test_leroy2011_low_metallicity(self):
        result = calc_deltaGD_from_metalZ_following_Leroy2011(7.5)
        assert result > 0

    def test_leroy2011_high_metallicity(self):
        result = calc_deltaGD_from_metalZ_following_Leroy2011(9.5)
        assert result > 0

    def test_magdis2012_scalar(self):
        result = calc_deltaGD_from_metalZ_following_Magdis2012(8.67)
        assert result > 0

    def test_magdis2012_array(self):
        metalZ = np.array([8.0, 8.5, 8.67, 9.0])
        result = calc_deltaGD_from_metalZ_following_Magdis2012(metalZ)
        assert len(result) == 4
        assert np.all(result > 0)

    def test_remy_ruyer_2014a_scalar(self):
        result = calc_deltaGD_from_metalZ_following_RemyRuyer2014a(8.67)
        assert result > 0

    def test_remy_ruyer_2014a_below_knee(self):
        result = calc_deltaGD_from_metalZ_following_RemyRuyer2014a(7.5)
        assert result > 0

    def test_remy_ruyer_2014a_above_knee(self):
        result = calc_deltaGD_from_metalZ_following_RemyRuyer2014a(8.69)
        assert result > 0

    def test_remy_ruyer_2014a_array(self):
        metalZ = np.array([7.5, 8.0, 8.67, 9.0])
        result = calc_deltaGD_from_metalZ_following_RemyRuyer2014a(metalZ)
        assert len(result) == 4
        assert np.all(result > 0)

    def test_remy_ruyer_2014b_scalar(self):
        result = calc_deltaGD_from_metalZ_following_RemyRuyer2014b(8.67)
        assert result > 0

    def test_remy_ruyer_2014b_array(self):
        metalZ = np.array([7.5, 8.0, 8.67, 9.0])
        result = calc_deltaGD_from_metalZ_following_RemyRuyer2014b(metalZ)
        assert len(result) == 4
        assert np.all(result > 0)

    def test_genzel2015_scalar(self):
        result = calc_deltaGD_from_metalZ_following_Genzel2015(8.67)
        assert result > 0

    def test_genzel2015_same_as_leroy2011(self):
        metalZ = 8.67
        result_genzel = calc_deltaGD_from_metalZ_following_Genzel2015(metalZ)
        result_leroy = calc_deltaGD_from_metalZ_following_Leroy2011(metalZ)
        assert np.isclose(result_genzel, result_leroy)

    def test_deltaGD_decreases_with_metallicity(self):
        metalZ_low = 8.0
        metalZ_high = 9.0
        result_low = calc_deltaGD_from_metalZ_following_Leroy2011(metalZ_low)
        result_high = calc_deltaGD_from_metalZ_following_Leroy2011(metalZ_high)
        assert result_low > result_high
