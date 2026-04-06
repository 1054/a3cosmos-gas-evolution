import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'a3cosmos_gas_evolution', 'Common_Python_Code'))

import numpy as np
import pytest
from calc_gas_fraction import (
    calc_gas_fraction_Scoville2017,
    calc_gas_fraction_Tacconi2018,
    calc_gas_fraction_Tacconi2018_beta_2,
    calc_gas_fraction_A3COSMOS,
    calc_gas_fraction_A3COSMOS_with_Leslie2019_MS,
    calc_gas_fraction_A3COSMOS_with_Scoville2017_MS,
)


class TestCalcGasFraction:

    def test_scoville2017_return_fgas(self):
        result = calc_gas_fraction_Scoville2017(z=1.0, lgMstar=10.5, DeltaMS=0.0, return_fgas=True)
        assert 0.0 < result < 1.0

    def test_scoville2017_return_ratio(self):
        result = calc_gas_fraction_Scoville2017(z=1.0, lgMstar=10.5, DeltaMS=0.0, return_fgas=False)
        assert result > 0

    def test_scoville2017_z0(self):
        result = calc_gas_fraction_Scoville2017(z=0.0, lgMstar=10.5, DeltaMS=0.0)
        assert 0.0 < result < 1.0

    def test_scoville2017_high_z(self):
        result = calc_gas_fraction_Scoville2017(z=3.0, lgMstar=10.5, DeltaMS=0.0)
        assert 0.0 < result < 1.0

    def test_scoville2017_different_mass(self):
        result_low = calc_gas_fraction_Scoville2017(z=1.0, lgMstar=9.0, DeltaMS=0.0)
        result_high = calc_gas_fraction_Scoville2017(z=1.0, lgMstar=11.0, DeltaMS=0.0)
        assert result_low > result_high

    def test_tacconi2018_default(self):
        result = calc_gas_fraction_Tacconi2018(z=1.0, lgMstar=10.5, DeltaMS=0.0)
        assert 0.0 < result < 1.0

    def test_tacconi2018_redshift_2p5(self):
        result = calc_gas_fraction_Tacconi2018(z=2.5, lgMstar=10.5, DeltaMS=0.0, return_fgas=False)
        assert np.abs(result - 1.38467) < 1e-5

    def test_tacconi2018_beta_2(self):
        result = calc_gas_fraction_Tacconi2018_beta_2(z=1.0, lgMstar=10.5, DeltaMS=0.0)
        assert 0.0 < result < 1.0

    def test_tacconi2018_return_ratio(self):
        result = calc_gas_fraction_Tacconi2018(z=1.0, lgMstar=10.5, DeltaMS=0.0, return_fgas=False)
        assert result > 0

    def test_a3cosmos_z(self):
        result = calc_gas_fraction_A3COSMOS(z=1.0, lgMstar=10.5, DeltaMS=0.0)
        assert 0.0 < result < 1.0

    def test_a3cosmos_z_redshift_2p5(self):
        result = calc_gas_fraction_A3COSMOS(z=2.5, lgMstar=10.5, DeltaMS=0.0, return_fgas=False)
        assert np.abs(result - 1.639419) < 1e-5

    def test_a3cosmos_cosmic_age(self):
        result = calc_gas_fraction_A3COSMOS(cosmic_age=5.0, lgMstar=10.5, DeltaMS=0.0)
        assert 0.0 < result < 1.0

    def test_a3cosmos_return_ratio(self):
        result = calc_gas_fraction_A3COSMOS(z=1.0, lgMstar=10.5, DeltaMS=0.0, return_fgas=False)
        assert result > 0

    def test_a3cosmos_with_leslie2019(self):
        result = calc_gas_fraction_A3COSMOS_with_Leslie2019_MS(z=1.0, lgMstar=10.5, DeltaMS=0.0)
        assert 0.0 < result < 1.0

    def test_a3cosmos_with_scoville2017_MS(self):
        result = calc_gas_fraction_A3COSMOS_with_Scoville2017_MS(z=1.0, lgMstar=10.5, DeltaMS=0.0)
        assert 0.0 < result < 1.0

    def test_fgas_increases_with_z(self):
        result_z0 = calc_gas_fraction_Scoville2017(z=0.5, lgMstar=10.5, DeltaMS=0.0)
        result_z2 = calc_gas_fraction_Scoville2017(z=2.0, lgMstar=10.5, DeltaMS=0.0)
        assert result_z2 > result_z0

    def test_fgas_decreases_with_mass(self):
        result_lowM = calc_gas_fraction_Scoville2017(z=1.0, lgMstar=9.0, DeltaMS=0.0)
        result_highM = calc_gas_fraction_Scoville2017(z=1.0, lgMstar=11.0, DeltaMS=0.0)
        assert result_lowM > result_highM

    def test_fgas_increases_with_DeltaMS(self):
        result_D0 = calc_gas_fraction_Scoville2017(z=1.0, lgMstar=10.5, DeltaMS=0.0)
        result_D1 = calc_gas_fraction_Scoville2017(z=1.0, lgMstar=10.5, DeltaMS=1.0)
        assert result_D1 > result_D0
