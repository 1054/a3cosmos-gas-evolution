import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'a3cosmos_gas_evolution', 'Common_Python_Code'))

import numpy as np
import pytest
from calc_star_formation_law import (
    calc_Mmolgas_from_SFR_using_SF_law_Sargent2014,
    calc_SFR_from_Mmolgas_using_SF_law_Sargent2014,
)


class TestCalcStarFormationLaw:

    def test_calc_Mmolgas_sargent2014_scalar(self):
        SFR = 10.0
        DeltaMS = 0.0
        result = calc_Mmolgas_from_SFR_using_SF_law_Sargent2014(SFR, DeltaMS)
        assert np.isscalar(result) or result.shape == ()
        assert result > 0

    def test_calc_Mmolgas_sargent2014_array_SFR(self):
        SFR = np.array([1.0, 10.0, 100.0])
        DeltaMS = 0.0
        result = calc_Mmolgas_from_SFR_using_SF_law_Sargent2014(SFR, DeltaMS)
        assert len(result) == 3
        assert np.all(result > 0)

    def test_calc_Mmolgas_sargent2014_DeltaMS_low(self):
        SFR = 10.0
        DeltaMS = 0.0
        result = calc_Mmolgas_from_SFR_using_SF_law_Sargent2014(SFR, DeltaMS)
        assert result > 0

    def test_calc_Mmolgas_sargent2014_DeltaMS_high(self):
        SFR = 10.0
        DeltaMS = 1.0
        result = calc_Mmolgas_from_SFR_using_SF_law_Sargent2014(SFR, DeltaMS)
        assert result > 0

    def test_calc_Mmolgas_sargent2014_DeltaMS_array(self):
        SFR = 10.0
        DeltaMS = np.array([0.0, 0.5, 1.0])
        result = calc_Mmolgas_from_SFR_using_SF_law_Sargent2014(SFR, DeltaMS)
        assert len(result) == 3
        assert np.all(result > 0)

    def test_calc_SFR_from_Mmolgas_sargent2014_scalar(self):
        Mmolgas = 1e10
        DeltaMS = 0.0
        result = calc_SFR_from_Mmolgas_using_SF_law_Sargent2014(Mmolgas, DeltaMS)
        assert np.isscalar(result) or result.shape == ()
        assert result > 0

    def test_calc_SFR_from_Mmolgas_sargent2014_array(self):
        Mmolgas = np.array([1e9, 1e10, 1e11])
        DeltaMS = 0.0
        result = calc_SFR_from_Mmolgas_using_SF_law_Sargent2014(Mmolgas, DeltaMS)
        assert len(result) == 3
        assert np.all(result > 0)

    def test_roundtrip_consistency(self):
        SFR_original = 10.0
        DeltaMS = 0.5
        Mmol = calc_Mmolgas_from_SFR_using_SF_law_Sargent2014(SFR_original, DeltaMS)
        SFR_recovered = calc_SFR_from_Mmolgas_using_SF_law_Sargent2014(Mmol, DeltaMS)
        assert np.isclose(SFR_original, SFR_recovered, rtol=1e-5)

    def test_Mmolgas_increases_with_SFR(self):
        SFR_low = 1.0
        SFR_high = 100.0
        DeltaMS = 0.0
        result_low = calc_Mmolgas_from_SFR_using_SF_law_Sargent2014(SFR_low, DeltaMS)
        result_high = calc_Mmolgas_from_SFR_using_SF_law_Sargent2014(SFR_high, DeltaMS)
        assert result_high > result_low

    def test_Mmolgas_increases_with_DeltaMS(self):
        SFR = 10.0
        DeltaMS_low = 0.0
        DeltaMS_high = 1.0
        result_low = calc_Mmolgas_from_SFR_using_SF_law_Sargent2014(SFR, DeltaMS_low)
        result_high = calc_Mmolgas_from_SFR_using_SF_law_Sargent2014(SFR, DeltaMS_high)
        assert result_low > 0 and result_high > 0
