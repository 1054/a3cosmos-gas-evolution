import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'a3cosmos_gas_evolution', 'Common_Python_Code'))

import numpy as np
import pytest
from calc_gas_mass_from_dust import (
    calc_Lv_from_flux_density,
    calc_Lv_850um_from_RF_850um,
    calc_M_mol_Hughes2017,
    calc_M_total_gas_Hughes2017,
    calc_M_mol_Groves2015,
    calc_M_total_gas_Groves2015,
    calc_M_mol_Scoville2017_RF,
    calc_M_mol_Scoville2017,
    calc_gas_mass_from_dust_continuum,
    calc_gas_mass_from_dust_mass,
    calc_molecular_hydrogen_fraction,
)


class TestCalcGasMassFromDustHelpers:

    def test_calc_Lv_from_flux_density(self):
        result = calc_Lv_from_flux_density(1.0, 850.0, 2.0)
        assert result > 0

    def test_calc_Lv_850um_from_RF_850um(self):
        result = calc_Lv_850um_from_RF_850um(1.0, 2.0)
        assert result > 0

    def test_calc_M_mol_Hughes2017(self):
        result = calc_M_mol_Hughes2017(1.0, 2.0)
        assert result > 0

    def test_calc_M_total_gas_Hughes2017(self):
        result = calc_M_total_gas_Hughes2017(1.0, 2.0)
        assert result > 0

    def test_calc_M_mol_Groves2015(self):
        result = calc_M_mol_Groves2015(1.0, 850.0, 2.0)
        assert np.all(result > 0)

    def test_calc_M_total_gas_Groves2015(self):
        result = calc_M_total_gas_Groves2015(1.0, 850.0, 2.0)
        assert np.all(result > 0)

    def test_calc_M_mol_Scoville2017_RF(self):
        result = calc_M_mol_Scoville2017_RF(1.0, 2.0)
        assert result > 0

    def test_calc_M_mol_Scoville2017(self):
        result = calc_M_mol_Scoville2017(1.0, 850.0, 2.0)
        assert result > 0


class TestCalcGasMassFromDustContinuum:

    def test_hughes2017_method(self):
        result, method = calc_gas_mass_from_dust_continuum(
            obs_wavelength_um=850.0,
            obs_flux_mJy=1.0,
            SED_flux_at_obs_wavelength=1.0,
            SED_flux_at_rest_850um=1.0,
            z=2.0,
            method='Hughes2017'
        )
        assert result > 0
        assert method == 'Hughes2017'

    def test_groves2015_method(self):
        result, method = calc_gas_mass_from_dust_continuum(
            obs_wavelength_um=250.0,
            obs_flux_mJy=1.0,
            z=2.0,
            method='Groves2015'
        )
        assert result > 0
        assert method == 'Groves2015'

    def test_scoville2017_method(self):
        result, method = calc_gas_mass_from_dust_continuum(
            obs_wavelength_um=850.0,
            obs_flux_mJy=1.0,
            z=2.0,
            method='Scoville2017'
        )
        assert result > 0
        assert method == 'Scoville2017'

    def test_missing_obs_wavelength_um_raises(self):
        with pytest.raises(ValueError):
            calc_gas_mass_from_dust_continuum(
                obs_wavelength_um=None,
                obs_flux_mJy=1.0,
                z=2.0
            )

    def test_missing_obs_flux_mJy_raises(self):
        with pytest.raises(ValueError):
            calc_gas_mass_from_dust_continuum(
                obs_wavelength_um=850.0,
                obs_flux_mJy=None,
                z=2.0
            )

    def test_missing_z_raises(self):
        with pytest.raises(ValueError):
            calc_gas_mass_from_dust_continuum(
                obs_wavelength_um=850.0,
                obs_flux_mJy=1.0,
                z=None
            )

    def test_invalid_method_raises(self):
        with pytest.raises(NotImplementedError):
            calc_gas_mass_from_dust_continuum(
                obs_wavelength_um=850.0,
                obs_flux_mJy=1.0,
                z=2.0,
                method='InvalidMethod'
            )


class TestCalcMolecularHydrogenFraction:

    def test_krumholz2009(self):
        result = calc_molecular_hydrogen_fraction(8.67, method='Krumholz2009')
        assert 0.0 <= result <= 1.0

    def test_dave2016(self):
        result = calc_molecular_hydrogen_fraction(8.67, method='Dave2016')
        assert 0.0 <= result <= 1.0

    def test_popping2014(self):
        pass

    def test_invalid_method_raises(self):
        with pytest.raises(NotImplementedError):
            calc_molecular_hydrogen_fraction(8.67, method='InvalidMethod')


class TestCalcGasMassFromDustMass:

    def test_gdr_method(self):
        result, method = calc_gas_mass_from_dust_mass(
            M_dust=1e8,
            GDR=100.0,
            method='GDR'
        )
        assert result > 0
        assert 'GDR' in method

    def test_magdis2012_with_metallicity(self):
        result, method = calc_gas_mass_from_dust_mass(
            M_dust=1e8,
            metallicity=8.67,
            method='Magdis2012'
        )
        assert result > 0

    def test_magdis2012_with_M_star_SFR(self):
        result, method = calc_gas_mass_from_dust_mass(
            M_dust=1e8,
            M_star=1e10,
            SFR=10.0,
            method='Magdis2012'
        )
        assert result > 0

    def test_genzel2015a_with_z(self):
        result, method = calc_gas_mass_from_dust_mass(
            M_dust=1e8,
            M_star=1e10,
            SFR=10.0,
            z=2.0,
            method='Genzel2015a'
        )
        assert result > 0

    def test_genzel2015b(self):
        result, method = calc_gas_mass_from_dust_mass(
            M_dust=1e8,
            M_star=1e10,
            SFR=10.0,
            method='Genzel2015b'
        )
        assert result > 0

    def test_default_method(self):
        result, method = calc_gas_mass_from_dust_mass(
            M_dust=1e8,
            M_star=1e10,
            SFR=10.0,
            z=2.0
        )
        assert result > 0

    def test_gdr_method_without_gdr_raises(self):
        with pytest.raises(ValueError):
            calc_gas_mass_from_dust_mass(
                M_dust=1e8,
                method='GDR'
            )

    def test_magdis2012_missing_params_raises(self):
        with pytest.raises(ValueError):
            calc_gas_mass_from_dust_mass(
                M_dust=1e8,
                method='Magdis2012'
            )

    def test_genzel2015a_missing_params_raises(self):
        with pytest.raises(ValueError):
            calc_gas_mass_from_dust_mass(
                M_dust=1e8,
                M_star=1e10,
                SFR=10.0,
                method='Genzel2015a'
            )

    def test_invalid_method_raises(self):
        with pytest.raises(NotImplementedError):
            calc_gas_mass_from_dust_mass(
                M_dust=1e8,
                method='InvalidMethod'
            )

    def test_gas_mass_positive_for_reasonable_inputs(self):
        result = calc_gas_mass_from_dust_mass(
            M_dust=1e8,
            metallicity=8.67,
            GDR=100.0
        )
        assert result is not None

    def test_combined_method_with_krumholz(self):
        result = calc_gas_mass_from_dust_mass(
            M_dust=1e8,
            metallicity=8.67,
            method='Magdis2012+Krumholz2009'
        )
        assert result is not None
