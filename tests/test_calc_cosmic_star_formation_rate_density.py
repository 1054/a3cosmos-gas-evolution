import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'a3cosmos_gas_evolution', 'Common_Python_Code'))

import numpy as np
import pytest
from calc_cosmic_star_formation_rate_density import (
    convert_age_to_z,
    calc_cosmic_star_formation_rate_density_MadauDickinson2014,
    calc_CSFRD_Madau2014,
    calc_cosmic_star_formation_rate_density_Abramson2016,
    calc_cosmic_star_formation_rate_density_Gladders2013,
    calc_CSFRD_Gladders2013,
    calc_CSFRD_Abramson2016,
    calc_cosmic_star_formation_rate_density_Liu2018,
    calc_CSFRD_Liu2018,
)


class TestCalcCosmicSFRD:

    def test_convert_age_to_z_scalar(self):
        cosmoAge = 5.0
        result = convert_age_to_z(cosmoAge)
        assert result > 0

    def test_convert_age_to_z_array(self):
        cosmoAge = np.array([2.0, 5.0, 10.0])
        result = convert_age_to_z(cosmoAge)
        assert len(result) == 3
        assert np.all(result > 0)

    def test_madau_dickinson2014_scalar(self):
        result = calc_cosmic_star_formation_rate_density_MadauDickinson2014(1.0)
        assert result > 0

    def test_madau_dickinson2014_array(self):
        z = np.array([0.5, 1.0, 2.0, 3.0])
        result = calc_cosmic_star_formation_rate_density_MadauDickinson2014(z)
        assert len(result) == 4
        assert np.all(result > 0)

    def test_madau_dickinson2014_z0(self):
        result = calc_cosmic_star_formation_rate_density_MadauDickinson2014(0.0)
        assert result > 0

    def test_madau_dickinson2014_high_z(self):
        result = calc_cosmic_star_formation_rate_density_MadauDickinson2014(5.0)
        assert result > 0

    def test_CSFRD_Madau2014(self):
        result = calc_CSFRD_Madau2014(1.0)
        assert result > 0

    def test_abramson2016_scalar(self):
        result = calc_cosmic_star_formation_rate_density_Abramson2016(1.0)
        assert result > 0

    def test_abramson2016_z0(self):
        result = calc_cosmic_star_formation_rate_density_Abramson2016(0.0)
        assert result > 0

    def test_abramson2016_high_z(self):
        result = calc_cosmic_star_formation_rate_density_Abramson2016(3.0)
        assert result > 0

    def test_Gladders2013_alias(self):
        result_gladders = calc_cosmic_star_formation_rate_density_Gladders2013(1.0)
        result_abramson = calc_cosmic_star_formation_rate_density_Abramson2016(1.0)
        assert np.isclose(result_gladders, result_abramson)

    def test_CSFRD_Gladders2013_alias(self):
        result = calc_CSFRD_Gladders2013(1.0)
        assert result > 0

    def test_CSFRD_Abramson2016_alias(self):
        result = calc_CSFRD_Abramson2016(1.0)
        assert result > 0

    def test_liu2018_double_powerlaw(self):
        result = calc_cosmic_star_formation_rate_density_Liu2018(1.0, shape='double-powerlaw')
        assert result > 0

    def test_liu2018_log_normal(self):
        result = calc_cosmic_star_formation_rate_density_Liu2018(1.0, shape='log-normal')
        assert result > 0

    def test_liu2018_invalid_shape(self):
        with pytest.raises(ValueError):
            calc_cosmic_star_formation_rate_density_Liu2018(1.0, shape='invalid')

    def test_CSFRD_Liu2018_alias(self):
        result = calc_CSFRD_Liu2018(1.0)
        assert result > 0

    def test_CSFRD_Liu2018_with_shape(self):
        result = calc_CSFRD_Liu2018(1.0, shape='log-normal')
        assert result > 0

    def test_rhoSFR_positive_at_all_z(self):
        z_values = np.array([0.0, 0.5, 1.0, 2.0, 3.0, 5.0])
        for z in z_values:
            result = calc_cosmic_star_formation_rate_density_MadauDickinson2014(z)
            assert result >= 0

    def test_rhoSFR_has_peak_at_intermediate_z(self):
        z = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0])
        rho = calc_cosmic_star_formation_rate_density_MadauDickinson2014(z)
        peak_idx = np.argmax(rho)
        assert 0 < peak_idx < len(z) - 1
