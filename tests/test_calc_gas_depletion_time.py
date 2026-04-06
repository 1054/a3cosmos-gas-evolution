import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'a3cosmos_gas_evolution', 'Common_Python_Code'))

import numpy as np
import pytest
from calc_gas_depletion_time import (
    calc_gas_depletion_time_Scoville2017,
    calc_gas_depletion_time_Tacconi2018,
    calc_gas_depletion_time_A3COSMOS,
)


class TestCalcGasDepletionTime:

    def test_scoville2017_scalar(self):
        result = calc_gas_depletion_time_Scoville2017(z=1.0, lgMstar=10.5, DeltaMS=0.0)
        assert result > 0

    def test_scoville2017_z0(self):
        result = calc_gas_depletion_time_Scoville2017(z=0.0, lgMstar=10.5, DeltaMS=0.0)
        assert result > 0

    def test_scoville2017_high_z(self):
        result = calc_gas_depletion_time_Scoville2017(z=3.0, lgMstar=10.5, DeltaMS=0.0)
        assert result > 0

    def test_scoville2017_different_mass(self):
        result_low = calc_gas_depletion_time_Scoville2017(z=1.0, lgMstar=9.0, DeltaMS=0.0)
        result_high = calc_gas_depletion_time_Scoville2017(z=1.0, lgMstar=11.0, DeltaMS=0.0)
        assert result_low > 0 and result_high > 0

    def test_scoville2017_different_DeltaMS(self):
        result_D0 = calc_gas_depletion_time_Scoville2017(z=1.0, lgMstar=10.5, DeltaMS=0.0)
        result_D1 = calc_gas_depletion_time_Scoville2017(z=1.0, lgMstar=10.5, DeltaMS=1.0)
        assert result_D1 < result_D0

    def test_tacconi2018_scalar(self):
        result = calc_gas_depletion_time_Tacconi2018(z=1.0, lgMstar=10.5, DeltaMS=0.0)
        assert result > 0

    def test_tacconi2018_z0(self):
        result = calc_gas_depletion_time_Tacconi2018(z=0.0, lgMstar=10.5, DeltaMS=0.0)
        assert result > 0

    def test_tacconi2018_high_z(self):
        result = calc_gas_depletion_time_Tacconi2018(z=3.0, lgMstar=10.5, DeltaMS=0.0)
        assert result > 0

    def test_a3cosmos_z(self):
        result = calc_gas_depletion_time_A3COSMOS(z=1.0, lgMstar=10.5, DeltaMS=0.0)
        assert result > 0

    def test_a3cosmos_cosmic_age(self):
        result = calc_gas_depletion_time_A3COSMOS(cosmic_age=5.0, lgMstar=10.5, DeltaMS=0.0)
        assert result > 0

    def test_a3cosmos_z_and_cosmic_age_consistent(self):
        import apply_cosmology
        cosmo = apply_cosmology.cosmo
        z = 1.0
        cosmic_age = cosmo.age(z).value
        result_z = calc_gas_depletion_time_A3COSMOS(z=z, lgMstar=10.5, DeltaMS=0.0)
        result_age = calc_gas_depletion_time_A3COSMOS(cosmic_age=cosmic_age, lgMstar=10.5, DeltaMS=0.0)
        assert np.isclose(result_z, result_age)

    def test_tauDepl_decreases_with_z(self):
        result_z0 = calc_gas_depletion_time_Scoville2017(z=0.5, lgMstar=10.5, DeltaMS=0.0)
        result_z2 = calc_gas_depletion_time_Scoville2017(z=2.0, lgMstar=10.5, DeltaMS=0.0)
        assert result_z2 < result_z0

    def test_tauDepl_increases_with_DeltaMS(self):
        result_D0 = calc_gas_depletion_time_Scoville2017(z=1.0, lgMstar=10.5, DeltaMS=0.0)
        result_D1 = calc_gas_depletion_time_Scoville2017(z=1.0, lgMstar=10.5, DeltaMS=1.0)
        assert result_D1 < result_D0
