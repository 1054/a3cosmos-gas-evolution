import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'a3cosmos_gas_evolution', 'Common_Python_Code'))

import numpy as np
import pytest
from calc_alpha_CO import (
    calc_alphaCO_from_metalZ_following_Wilson1995,
    calc_alphaCO_from_metalZ_following_Magdis2012,
    calc_alphaCO_from_metalZ_following_Genzel2012,
    calc_alphaCO_from_metalZ_following_Genzel2015a,
    calc_alphaCO_from_metalZ_following_Genzel2015b,
    calc_alphaCO_from_metalZ_following_Bolatto2013,
    calc_alphaCO_from_metalZ_following_Accurso2017,
    calc_alphaCO_from_metalZ_following_Bertemes2018,
    calc_alphaCO_from_metalZ_following_Tacconi2018,
    calc_alphaCO_from_metalZ_following_Boselli2014,
)


class TestCalcAlphaCO:

    def test_wilson1995_scalar(self):
        result = calc_alphaCO_from_metalZ_following_Wilson1995(8.67)
        assert np.isscalar(result)
        assert result > 0

    def test_wilson1995_array(self):
        metalZ = np.array([8.0, 8.5, 8.67, 9.0])
        result = calc_alphaCO_from_metalZ_following_Wilson1995(metalZ)
        assert len(result) == 4
        assert np.all(result > 0)

    def test_wilson1995_solar_metalZ(self):
        result = calc_alphaCO_from_metalZ_following_Wilson1995(8.67)
        assert result > 0

    def test_magdis2012_scalar(self):
        result = calc_alphaCO_from_metalZ_following_Magdis2012(8.67)
        assert np.isscalar(result) or result.shape == ()

    def test_magdis2012_array(self):
        metalZ = np.array([8.0, 8.5, 8.67, 9.0])
        result = calc_alphaCO_from_metalZ_following_Magdis2012(metalZ)
        assert len(result) == 4
        assert np.all(result > 0)

    def test_genzel2012_scalar(self):
        result = calc_alphaCO_from_metalZ_following_Genzel2012(8.67)
        assert np.isscalar(result) or result.shape == ()

    def test_genzel2012_array(self):
        metalZ = np.array([8.0, 8.5, 8.67, 9.0])
        result = calc_alphaCO_from_metalZ_following_Genzel2012(metalZ)
        assert len(result) == 4
        assert np.all(result > 0)

    def test_genzel2015a_default_params(self):
        result = calc_alphaCO_from_metalZ_following_Genzel2015a(8.67)
        assert result > 0

    def test_genzel2015a_explicit_params(self):
        result = calc_alphaCO_from_metalZ_following_Genzel2015a(8.67, metalZ_solar=8.67, alphaCO_MilkyWay=4.36)
        assert result > 0

    def test_genzel2015b_default_params(self):
        result = calc_alphaCO_from_metalZ_following_Genzel2015b(8.67)
        assert result > 0

    def test_genzel2015b_explicit_params(self):
        result = calc_alphaCO_from_metalZ_following_Genzel2015b(8.67, metalZ_solar=8.67, alphaCO_MilkyWay=4.36)
        assert result > 0

    def test_genzel2015_consistency(self):
        result_a = calc_alphaCO_from_metalZ_following_Genzel2015a(8.67)
        result_b = calc_alphaCO_from_metalZ_following_Genzel2015b(8.67)
        assert np.isscalar(result_a) or result_a.shape == ()
        assert np.isscalar(result_b) or result_b.shape == ()

    def test_bolatto2013_high_sigma(self):
        result = calc_alphaCO_from_metalZ_following_Bolatto2013(8.67, Sigma_total=150.0)
        assert result > 0

    def test_bolatto2013_low_sigma(self):
        result = calc_alphaCO_from_metalZ_following_Bolatto2013(8.67, Sigma_total=50.0)
        assert result > 0

    def test_bolatto2013_array_input(self):
        Sigma_total = np.array([50.0, 100.0, 150.0])
        result = calc_alphaCO_from_metalZ_following_Bolatto2013(8.67, Sigma_total)
        assert len(result) == 3
        assert np.all(result > 0)

    def test_accuro2017_scalar(self):
        result = calc_alphaCO_from_metalZ_following_Accurso2017(8.67, DeltaMS=0.0)
        assert result > 0

    def test_accuro2017_array_metalZ(self):
        metalZ = np.array([8.0, 8.5, 8.67, 9.0])
        result = calc_alphaCO_from_metalZ_following_Accurso2017(metalZ, DeltaMS=0.0)
        assert len(result) == 4
        assert np.all(result > 0)

    def test_bertemes2018_scalar(self):
        result = calc_alphaCO_from_metalZ_following_Bertemes2018(8.67)
        assert result > 0

    def test_tacconi2018_scalar(self):
        result = calc_alphaCO_from_metalZ_following_Tacconi2018(8.67)
        assert result > 0

    def test_boselli2014_lgL_Hband(self):
        result = calc_alphaCO_from_metalZ_following_Boselli2014(10.0)
        assert result > 0

    def test_boselli2014_array(self):
        lgL_Hband = np.array([9.0, 10.0, 11.0])
        result = calc_alphaCO_from_metalZ_following_Boselli2014(lgL_Hband)
        assert len(result) == 3
        assert np.all(result > 0)

    def test_alphaCO_decreases_with_metalZ(self):
        metalZ_low = 8.0
        metalZ_high = 9.5
        result_low = calc_alphaCO_from_metalZ_following_Wilson1995(metalZ_low)
        result_high = calc_alphaCO_from_metalZ_following_Wilson1995(metalZ_high)
        assert result_low > result_high
