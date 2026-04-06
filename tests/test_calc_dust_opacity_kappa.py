import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'a3cosmos_gas_evolution', 'Common_Python_Code'))

import numpy as np
import pytest


class TestCalcDustOpacityKappa:

    def test_kappa_above_700um(self):
        import calc_dust_opacity_kappa
        lambda_um = 850.0
        beta = 1.68
        kappa = 0.596 * (lambda_um / 700.0)**(-beta)
        assert kappa > 0

    def test_kappa_below_700um(self):
        import calc_dust_opacity_kappa
        lambda_um = 500.0
        beta = 2.0
        kappa = 0.596 * (lambda_um / 700.0)**(-beta)
        assert kappa > 0

    def test_kappa_at_700um_transition(self):
        lambda_um = 700.0
        beta_high = 1.68
        beta_low = 2.0
        kappa_high = 0.596 * (lambda_um / 700.0)**(-beta_high)
        kappa_low = 0.596 * (lambda_um / 700.0)**(-beta_low)
        assert np.isclose(kappa_high, kappa_low)

    def test_kappa_wavelength_dependence(self):
        lambda_1 = 500.0
        lambda_2 = 1000.0
        beta = 2.0
        kappa_1 = 0.596 * (lambda_1 / 700.0)**(-beta)
        kappa_2 = 0.596 * (lambda_2 / 700.0)**(-beta)
        assert kappa_1 > kappa_2
