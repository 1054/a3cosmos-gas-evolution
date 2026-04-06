import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'a3cosmos_gas_evolution', 'Common_Python_Code'))

import numpy as np
import pytest
from astropy import units as u


class TestCalcPlanckLaw:

    def test_blackbody_nu_at_25K(self):
        from calc_planck_law import blackbody_nu
        wavelength = 100.0 * u.um
        T = 25.0 * u.K
        result = blackbody_nu(wavelength, T)
        assert result.value >= 0

    def test_blackbody_lambda_at_25K(self):
        from calc_planck_law import blackbody_lambda
        wavelength = 100.0 * u.um
        T = 25.0 * u.K
        result = blackbody_lambda(wavelength, T)
        assert result.value >= 0

    def test_blackbody_nu_at_50K(self):
        from calc_planck_law import blackbody_nu
        wavelength = 100.0 * u.um
        T = 50.0 * u.K
        result = blackbody_nu(wavelength, T)
        assert result.value >= 0

    def test_blackbody_returns_positive(self):
        from calc_planck_law import blackbody_nu, blackbody_lambda
        wavelengths = np.array([50.0, 100.0, 200.0]) * u.um
        T = 25.0 * u.K
        for w in wavelengths:
            result_nu = blackbody_nu(w, T)
            result_lambda = blackbody_lambda(w, T)
            assert result_nu.value >= 0
            assert result_lambda.value >= 0

    def test_blackbody_high_temp_greater_than_low_temp(self):
        from calc_planck_law import blackbody_nu
        wavelength = 100.0 * u.um
        T_low = 20.0 * u.K
        T_high = 30.0 * u.K
        result_low = blackbody_nu(wavelength, T_low)
        result_high = blackbody_nu(wavelength, T_high)
        assert result_high.value > result_low.value
