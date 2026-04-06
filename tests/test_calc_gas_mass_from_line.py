import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'a3cosmos_gas_evolution', 'Common_Python_Code'))

import numpy as np
import pytest

from calc_gas_mass_from_line import (
    ALPHA_CII,
    ALPHA_CO_MW,
    ALPHA_CO_SB,
    calc_gas_mass_from_line_flux,
)


class TestCalcGasMassFromLineConstants:

    def test_alpha_CII_value(self):
        assert ALPHA_CII == 30.0

    def test_alpha_CO_MW_value(self):
        assert ALPHA_CO_MW == 4.3

    def test_alpha_CO_SB_value(self):
        assert ALPHA_CO_SB == 0.8


class TestCalcGasMassFromLineBoundaryCases:

    def test_invalid_z_None_returns_nan(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=None, lgMstar=10.5, freq=115.27, flux=1.0, line_name='CO(1-0)'
        )
        assert np.isnan(lgMgas)

    def test_invalid_z_empty_string_returns_nan(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z='', lgMstar=10.5, freq=115.27, flux=1.0, line_name='CO(1-0)'
        )
        assert np.isnan(lgMgas)

    def test_invalid_z_negative_returns_nan(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=-1.0, lgMstar=10.5, freq=115.27, flux=1.0, line_name='CO(1-0)'
        )
        assert np.isnan(lgMgas)

    def test_invalid_z_zero_returns_nan(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=0.0, lgMstar=10.5, freq=115.27, flux=1.0, line_name='CO(1-0)'
        )
        assert np.isnan(lgMgas)

    def test_invalid_z_nan_returns_nan(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=np.nan, lgMstar=10.5, freq=115.27, flux=1.0, line_name='CO(1-0)'
        )
        assert np.isnan(lgMgas)


class TestCalcGasMassFromLineCO:

    def test_co10_line_z2(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=2.0, lgMstar=10.5, freq=115.27, flux=1.0, line_name='CO(1-0)'
        )
        assert not np.isnan(lgMgas)
        assert not np.isnan(lgLline)
        assert line_type == 'CO'
        assert line_unit == 'K km s-1 pc2'
        assert lgMgas > 0
        assert lgLline > 0

    def test_co10_line_z1(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=1.0, lgMstar=10.5, freq=115.27 / 2.0, flux=1.0, line_name='CO(1-0)'
        )
        assert not np.isnan(lgMgas)
        assert line_type == 'CO'

    def test_co21_line_z2(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=2.0, lgMstar=10.5, freq=230.54, flux=1.0, line_name='CO(2-1)'
        )
        assert not np.isnan(lgMgas)
        assert line_type == 'CO'

    def test_co32_line_z2(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=2.0, lgMstar=10.5, freq=345.796, flux=1.0, line_name='CO(3-2)'
        )
        assert not np.isnan(lgMgas)
        assert line_type == 'CO'

    def test_co_line_high_z(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=5.0, lgMstar=10.5, freq=115.27 / 6.0, flux=1.0, line_name='CO(1-0)'
        )
        assert not np.isnan(lgMgas)
        assert lgMgas > 0

    def test_co_line_with_custom_alpha_CO(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=2.0, lgMstar=10.5, freq=115.27, flux=1.0, line_name='CO(1-0)',
            alpha_CO=3.0
        )
        assert not np.isnan(lgMgas)

    def test_co_line_without_logMstar(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=2.0, lgMstar=None, freq=115.27, flux=1.0, line_name='CO(1-0)'
        )
        assert not np.isnan(lgMgas)


class TestCalcGasMassFromLineCII:

    def test_cii_line_z6(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=6.0, lgMstar=10.5, freq=1900.0, flux=1.0, line_name='[CII]'
        )
        assert not np.isnan(lgMgas)
        assert not np.isnan(lgLline)
        assert line_type == 'CII'
        assert line_unit == 'L_sun'
        assert lgMgas > 0
        assert lgLline > 0

    def test_cii_line_alt_name(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=6.0, lgMstar=10.5, freq=1900.0, flux=1.0, line_name='CII'
        )
        assert not np.isnan(lgMgas)
        assert line_type == 'CII'

    def test_cii_line_high_z(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=7.0, lgMstar=10.5, freq=1900.0 / 8.0, flux=1.0, line_name='[CII]'
        )
        assert not np.isnan(lgMgas)
        assert lgMgas > 0

    def test_cii_line_with_custom_alpha_CII(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=6.0, lgMstar=10.5, freq=1900.0, flux=1.0, line_name='[CII]',
            alpha_CII=50.0
        )
        assert not np.isnan(lgMgas)


class TestCalcGasMassFromLineFlux:

    def test_different_flux_values(self):
        fluxes = [0.1, 1.0, 10.0]
        for flux in fluxes:
            lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
                z=2.0, lgMstar=10.5, freq=115.27, flux=flux, line_name='CO(1-0)'
            )
            assert not np.isnan(lgMgas)
            assert lgMgas > 0

    def test_flux_zero_returns_nan(self):
        lgMgas, lgLline, line_unit, line_type = calc_gas_mass_from_line_flux(
            z=2.0, lgMstar=10.5, freq=115.27, flux=0.0, line_name='CO(1-0)'
        )
        assert np.isnan(lgMgas)


class TestCalcGasMassFromLineOutput:

    def test_output_is_4_tuple(self):
        result = calc_gas_mass_from_line_flux(
            z=2.0, lgMstar=10.5, freq=115.27, flux=1.0, line_name='CO(1-0)'
        )
        assert len(result) == 4

    def test_co_and_cii_output_shapes_consistent(self):
        result_co = calc_gas_mass_from_line_flux(
            z=2.0, lgMstar=10.5, freq=115.27, flux=1.0, line_name='CO(1-0)'
        )
        result_cii = calc_gas_mass_from_line_flux(
            z=6.0, lgMstar=10.5, freq=1900.0, flux=1.0, line_name='[CII]'
        )
        assert len(result_co) == len(result_cii) == 4

    def test_output_shapes_consistent(self):
        result_co = calc_gas_mass_from_line_flux(
            z=2.0, lgMstar=10.5, freq=115.27, flux=1.0, line_name='CO(1-0)'
        )
        result_cii = calc_gas_mass_from_line_flux(
            z=6.0, lgMstar=10.5, freq=1900.0, flux=1.0, line_name='[CII]'
        )
        assert len(result_co) == len(result_cii) == 4

    def test_lgLline_increases_with_flux(self):
        lgLline_low = calc_gas_mass_from_line_flux(
            z=2.0, lgMstar=10.5, freq=115.27, flux=0.5, line_name='CO(1-0)'
        )[1]
        lgLline_high = calc_gas_mass_from_line_flux(
            z=2.0, lgMstar=10.5, freq=115.27, flux=1.0, line_name='CO(1-0)'
        )[1]
        assert lgLline_high > lgLline_low
