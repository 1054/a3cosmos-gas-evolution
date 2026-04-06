import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'a3cosmos_gas_evolution', 'Common_Python_Code'))

import numpy as np
import pytest
from calc_galaxy_main_sequence import (
    calc_Sargent2014_sSFR,
    calc_Speagle2014_sSFR,
    calc_Genzel2015_sSFR,
    calc_Scoville2017_sSFR,
    calc_sSFR_MS,
    calc_SFR_MS_Scoville2017,
    calc_SFR_MS_Tacconi2018,
    calc_SFR_MS_Speagle2014,
    calc_SFR_MS_Speagle2014_with_t,
    calc_SFR_MS_Speagle2014_best_fit_with_t_rewritten_by_dzliu,
    calc_SFR_MS_Speagle2014_with_z,
    calc_SFR_MS_Speagle2014_best_fit_with_z,
    calc_SFR_MS_Sargent2014,
    calc_SFR_MS_Bethermin2015,
    calc_SFR_MS_Schreiber2015,
    calc_SFR_MS_Genzel2015,
    calc_SFR_MS_Whitaker2014,
    calc_SFR_MS_Lee2015,
    calc_SFR_MS_Tomczak2016,
    calc_SFR_MS_Pearson2018,
    calc_SFR_MS_Leslie20180901,
    calc_SFR_MS_Leslie20190111,
    calc_SFR_MS_Leslie20190515,
    calc_SFR_MS_Leslie20190710,
    calc_SFR_MS_Leslie20191212,
    calc_SFR_MS_Leslie20200306,
)


class TestCalcGalaxyMainSequence:

    def test_calc_sSFR_Sargent2014_scalar(self):
        result = calc_Sargent2014_sSFR(1.0, lgMstar=10.5, DeltaMS=0.0)
        assert result > 0

    def test_calc_sSFR_Sargent2014_array(self):
        z = np.array([0.5, 1.0, 2.0])
        result = calc_Sargent2014_sSFR(z, lgMstar=10.5, DeltaMS=0.0)
        assert len(result) == 3
        assert np.all(result > 0)

    def test_calc_sSFR_Speagle2014_scalar(self):
        result = calc_Speagle2014_sSFR(5.0, lgMstar=10.5, DeltaMS=0.0)
        assert result > 0

    def test_calc_sSFR_Genzel2015_scalar(self):
        result = calc_Genzel2015_sSFR(1.0, lgMstar=10.5, DeltaMS=0.0)
        assert result > 0

    def test_calc_sSFR_Scoville2017_scalar(self):
        result = calc_Scoville2017_sSFR(1.0, lgMstar=10.5, DeltaMS=0.0)
        assert result > 0

    def test_calc_sSFR_MS_scalar(self):
        result = calc_sSFR_MS(10.5, z=1.0)
        assert result > 0

    def test_calc_SFR_MS_Scoville2017_scalar(self):
        result = calc_SFR_MS_Scoville2017(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Tacconi2018_scalar(self):
        result = calc_SFR_MS_Tacconi2018(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Tacconi2018_with_dzliu_conversion(self):
        result = calc_SFR_MS_Tacconi2018(1.0, lgMstar=10.5, with_dzliu_t_z_conversion=True)
        assert result > 0

    def test_calc_SFR_MS_Speagle2014_scalar(self):
        result = calc_SFR_MS_Speagle2014(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Speagle2014_with_t_scalar(self):
        result = calc_SFR_MS_Speagle2014_with_t(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Speagle2014_best_fit_with_t_rewritten_scalar(self):
        result = calc_SFR_MS_Speagle2014_best_fit_with_t_rewritten_by_dzliu(z=1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Speagle2014_with_z_scalar(self):
        result = calc_SFR_MS_Speagle2014_with_z(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Speagle2014_best_fit_with_z_scalar(self):
        result = calc_SFR_MS_Speagle2014_best_fit_with_z(z=1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Sargent2014_scalar(self):
        result = calc_SFR_MS_Sargent2014(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Bethermin2015_scalar(self):
        result = calc_SFR_MS_Bethermin2015(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Bethermin2015_high_z(self):
        result = calc_SFR_MS_Bethermin2015(3.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Schreiber2015_scalar(self):
        result = calc_SFR_MS_Schreiber2015(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Genzel2015_scalar(self):
        result = calc_SFR_MS_Genzel2015(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Whitaker2014_scalar(self):
        result = calc_SFR_MS_Whitaker2014(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Whitaker2014_use_equation_3(self):
        z = np.array([1.0, 1.5, 2.0, 2.5])
        result = calc_SFR_MS_Whitaker2014(z, lgMstar=10.5, use_equation_3=True)
        assert len(result) == 4

    def test_calc_SFR_MS_Lee2015_scalar(self):
        result = calc_SFR_MS_Lee2015(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Tomczak2016_scalar(self):
        result = calc_SFR_MS_Tomczak2016(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Pearson2018_scalar(self):
        result = calc_SFR_MS_Pearson2018(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Leslie20180901(self):
        result = calc_SFR_MS_Leslie20180901(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Leslie20190111(self):
        result = calc_SFR_MS_Leslie20190111(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Leslie20190515(self):
        result = calc_SFR_MS_Leslie20190515(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Leslie20190710(self):
        result = calc_SFR_MS_Leslie20190710(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Leslie20191212(self):
        result = calc_SFR_MS_Leslie20191212(1.0, lgMstar=10.5)
        assert result > 0

    def test_calc_SFR_MS_Leslie20200306(self):
        result = calc_SFR_MS_Leslie20200306(1.0, lgMstar=10.5)
        assert result > 0

    def test_sSFR_increases_with_z(self):
        result_z0 = calc_Sargent2014_sSFR(0.5, lgMstar=10.5, DeltaMS=0.0)
        result_z1 = calc_Sargent2014_sSFR(2.0, lgMstar=10.5, DeltaMS=0.0)
        assert result_z1 > result_z0

    def test_sSFR_decreases_with_mass(self):
        result_lowM = calc_Sargent2014_sSFR(1.0, lgMstar=9.0, DeltaMS=0.0)
        result_highM = calc_Sargent2014_sSFR(1.0, lgMstar=11.0, DeltaMS=0.0)
        assert result_lowM > result_highM

    def test_sSFR_increases_with_DeltaMS(self):
        result_D0 = calc_Sargent2014_sSFR(1.0, lgMstar=10.5, DeltaMS=0.0)
        result_D1 = calc_Sargent2014_sSFR(1.0, lgMstar=10.5, DeltaMS=1.0)
        assert result_D1 > result_D0
