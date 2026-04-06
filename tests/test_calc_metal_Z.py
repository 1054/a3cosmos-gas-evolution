import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'a3cosmos_gas_evolution', 'Common_Python_Code'))

import numpy as np
import pytest
from calc_metal_Z import (
    calc_metalZ_from_FMR_following_Kewley2008_PP04_O3N2,
    calc_metalZ_from_FMR_following_Kewley2008_PP04_N2,
    calc_metalZ_from_FMR_following_Kewley2008_KK04,
    calc_metalZ_from_FMR_following_Kewley2008_KD02,
    calc_metalZ_from_FMR_following_Mannucci2010,
    calc_metalZ_from_FMR_following_Mannucci2010_Eq2,
    calc_metalZ_from_FMR_following_Mannucci2010_Eq4,
    calc_metalZ_from_FMR_following_Mannucci2011,
    calc_metalZ_from_FMR_following_Maiolino2008,
    calc_metalZ_from_FMR_following_Genzel2015_Eq12a,
    calc_metalZ_from_FMR_following_Genzel2015_Eq12a_with_dzliu_limit,
    calc_metalZ_from_FMR_following_Genzel2015a,
    calc_metalZ_from_FMR_following_Genzel2015b,
    calc_metalZ_from_FMR_following_Genzel2015ab_combined_by_dzliu,
    calc_metalZ_from_FMR_following_Guo2016,
    convert_metalZ_M08_to_metalZ_PP04,
    convert_metalZ_D02_to_metalZ_PP04,
    convert_metalZ_KK04_to_metalZ_PP04,
)


class TestCalcMetalZ:

    def test_kewley2008_PP04_O3N2_scalar(self):
        result = calc_metalZ_from_FMR_following_Kewley2008_PP04_O3N2(1e10)
        assert np.isscalar(result) or result.shape == ()
        assert 7.0 < result < 10.0

    def test_kewley2008_PP04_O3N2_array(self):
        M_star = np.array([1e9, 1e10, 1e11])
        result = calc_metalZ_from_FMR_following_Kewley2008_PP04_O3N2(M_star)
        assert len(result) == 3

    def test_kewley2008_PP04_N2_scalar(self):
        result = calc_metalZ_from_FMR_following_Kewley2008_PP04_N2(1e10)
        assert np.isscalar(result) or result.shape == ()
        assert 7.0 < result < 10.0

    def test_kewley2008_KK04_scalar(self):
        result = calc_metalZ_from_FMR_following_Kewley2008_KK04(1e10)
        assert np.isscalar(result) or result.shape == ()

    def test_kewley2008_KD02_scalar(self):
        result = calc_metalZ_from_FMR_following_Kewley2008_KD02(1e10)
        assert np.isscalar(result) or result.shape == ()

    def test_mannucci2010_scalar(self):
        result = calc_metalZ_from_FMR_following_Mannucci2010(1e10, 1.0)
        assert np.isscalar(result) or result.shape == ()
        assert 7.0 < result < 10.0

    def test_mannucci2010_array_Mstar(self):
        M_star = np.array([1e9, 1e10, 1e11])
        SFR = 1.0
        result = calc_metalZ_from_FMR_following_Mannucci2010(M_star, SFR)
        assert len(result) == 3

    def test_mannucci2010_array_SFR(self):
        M_star = 1e10
        SFR = np.array([0.1, 1.0, 10.0])
        result = calc_metalZ_from_FMR_following_Mannucci2010(M_star, SFR)
        assert len(result) == 3

    def test_mannucci2010_Eq2_vs_Eq4(self):
        M_star = 1e10
        SFR = 1.0
        result_Eq2 = calc_metalZ_from_FMR_following_Mannucci2010_Eq2(M_star, SFR)
        result_Eq4 = calc_metalZ_from_FMR_following_Mannucci2010_Eq4(M_star, SFR)
        assert np.isscalar(result_Eq2) or result_Eq2.shape == ()
        assert np.isscalar(result_Eq4) or result_Eq4.shape == ()

    def test_mannucci2011_scalar(self):
        result = calc_metalZ_from_FMR_following_Mannucci2011(1e10, 1.0)
        assert np.isscalar(result) or result.shape == ()

    def test_maiolino2008_z0(self):
        result = calc_metalZ_from_FMR_following_Maiolino2008(1e10, 0.0)
        assert np.isscalar(result) or result.shape == ()

    def test_maiolino2008_high_z(self):
        result = calc_metalZ_from_FMR_following_Maiolino2008(1e10, 3.0)
        assert np.isscalar(result) or result.shape == ()

    def test_genzel2015_Eq12a_scalar(self):
        result = calc_metalZ_from_FMR_following_Genzel2015_Eq12a(1e10, 1.0)
        assert np.isscalar(result) or result.shape == ()

    def test_genzel2015_Eq12a_array(self):
        M_star = np.array([1e9, 1e10, 1e11])
        z = 1.0
        result = calc_metalZ_from_FMR_following_Genzel2015_Eq12a(M_star, z)
        assert len(result) == 3

    def test_genzel2015_Eq12a_with_dzliu_limit(self):
        result = calc_metalZ_from_FMR_following_Genzel2015_Eq12a_with_dzliu_limit(1e10, 1.0)
        assert np.isscalar(result) or result.shape == ()

    def test_genzel2015a_scalar(self):
        result = calc_metalZ_from_FMR_following_Genzel2015a(1e10, 1.0, 1.0)
        assert np.isscalar(result) or result.shape == ()

    def test_genzel2015b_scalar(self):
        result = calc_metalZ_from_FMR_following_Genzel2015b(1e10, 1.0)
        assert np.isscalar(result) or result.shape == ()

    def test_genzel2015ab_combined_scalar(self):
        result = calc_metalZ_from_FMR_following_Genzel2015ab_combined_by_dzliu(1e10, 1.0, 1.0)
        assert np.isscalar(result) or result.shape == ()

    def test_genzel2015ab_combined_array(self):
        M_star = np.array([1e9, 1e10, 1e11])
        SFR = 1.0
        z = 1.0
        result = calc_metalZ_from_FMR_following_Genzel2015ab_combined_by_dzliu(M_star, SFR, z)
        assert len(result) == 3

    def test_guo2016_scalar(self):
        result = calc_metalZ_from_FMR_following_Guo2016(1e10, 1.0)
        assert np.isscalar(result) or result.shape == ()

    def test_convert_M08_to_PP04(self):
        metalZ_M08 = 8.67
        result = convert_metalZ_M08_to_metalZ_PP04(metalZ_M08)
        assert np.isscalar(result) or result.shape == ()
        assert 7.0 < result < 10.0

    def test_convert_D02_to_PP04(self):
        metalZ_D02 = 8.67
        result = convert_metalZ_D02_to_metalZ_PP04(metalZ_D02)
        assert np.isscalar(result) or result.shape == ()

    def test_convert_KK04_to_PP04(self):
        metalZ_KK04 = 8.67
        result = convert_metalZ_KK04_to_metalZ_PP04(metalZ_KK04)
        assert np.isscalar(result) or result.shape == ()

    def test_metalZ_increases_with_stellar_mass(self):
        M_star_low = 1e9
        M_star_high = 1e11
        result_low = calc_metalZ_from_FMR_following_Kewley2008_PP04_O3N2(M_star_low)
        result_high = calc_metalZ_from_FMR_following_Kewley2008_PP04_O3N2(M_star_high)
        assert result_high > result_low
