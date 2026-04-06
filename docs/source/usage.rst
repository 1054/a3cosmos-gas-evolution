Usage
=====

Getting Help
------------

.. code-block:: python

    import a3cosmos_gas_evolution
    a3cosmos_gas_evolution.help()

Computing Cosmic Molecular Gas Density
--------------------------------------

Get the cosmic molecular gas density evolution curve (e.g., Fig. 15 of Liu et al. 2019b):

.. code-block:: python

    import a3cosmos_gas_evolution

    # A3COSMOS
    z, rho_mol_gas = a3cosmos_gas_evolution.get_cosmic_mol_gas_density_A3COSMOS()
    # Returns rho_mol_gas in solar mass per cubic mega parsec.

    # Tacconi et al. 2018
    z, rho_mol_gas = a3cosmos_gas_evolution.get_cosmic_mol_gas_density_Tacconi2018()

    # Scoville et al. 2017
    z, rho_mol_gas = a3cosmos_gas_evolution.get_cosmic_mol_gas_density_Scoville2017()

    # Plot the results
    a3cosmos_gas_evolution.plot_cosmic_mol_gas_density()

Computing Gas Fraction
----------------------

Compute the gas fraction, i.e., ``f_gas = M_molgas / (M_molgas + M_star)``:

.. code-block:: python

    import a3cosmos_gas_evolution

    # A3COSMOS
    f_gas = a3cosmos_gas_evolution.calc_gas_fraction_A3COSMOS(
        z=3.0, lgMstar=10.5, DeltaMS=0.5, return_fgas=True
    )

    # Tacconi et al. 2018
    f_gas = a3cosmos_gas_evolution.calc_gas_fraction_Tacconi2018(
        z=3.0, lgMstar=10.5, DeltaMS=0.5, return_fgas=True
    )

    # Scoville et al. 2017
    f_gas = a3cosmos_gas_evolution.calc_gas_fraction_Scoville2017(
        z=3.0, lgMstar=10.5, DeltaMS=0.5, return_fgas=True
    )

    # You can also use cosmic_age instead of z
    f_gas = a3cosmos_gas_evolution.calc_gas_fraction_A3COSMOS(
        cosmic_age=2.178, lgMstar=10.5, DeltaMS=0.5
    )

Computing Gas-to-Stellar Mass Ratio
-----------------------------------

Compute the gas-to-stellar mass ratio, i.e., ``mu_gas = M_molgas / M_star``:

.. code-block:: python

    import a3cosmos_gas_evolution

    mu_gas = a3cosmos_gas_evolution.calc_gas_fraction_A3COSMOS(
        z=3.0, lgMstar=10.5, DeltaMS=0.5, return_fgas=False
    )
    mu_gas = a3cosmos_gas_evolution.calc_gas_fraction_Tacconi2018(
        z=3.0, lgMstar=10.5, DeltaMS=0.5, return_fgas=False
    )
    mu_gas = a3cosmos_gas_evolution.calc_gas_fraction_Scoville2017(
        z=3.0, lgMstar=10.5, DeltaMS=0.5, return_fgas=False
    )

Computing Gas Depletion Time
----------------------------

Compute the gas depletion time, i.e., ``tau_depl = M_molgas / SFR``:

.. code-block:: python

    import a3cosmos_gas_evolution

    tau_depl = a3cosmos_gas_evolution.calc_gas_depletion_time_A3COSMOS(
        z=3.0, lgMstar=10.5, DeltaMS=0.0
    )
    tau_depl = a3cosmos_gas_evolution.calc_gas_depletion_time_Tacconi2018(
        z=3.0, lgMstar=10.5, DeltaMS=0.0
    )
    tau_depl = a3cosmos_gas_evolution.calc_gas_depletion_time_Scoville2017(
        z=3.0, lgMstar=10.5, DeltaMS=0.0
    )

Computing Galaxy Main-Sequence SFR
----------------------------------

Compute the galaxy main-sequence star formation rate:

.. code-block:: python

    import a3cosmos_gas_evolution

    SFR_MS = a3cosmos_gas_evolution.calc_SFR_MS_Speagle2014(z=3.0, lgMstar=10.5)
    SFR_MS = a3cosmos_gas_evolution.calc_SFR_MS_Scoville2017(z=3.0, lgMstar=10.5)
    SFR_MS = a3cosmos_gas_evolution.calc_SFR_MS_Tacconi2018(z=3.0, lgMstar=10.5)
    SFR_MS = a3cosmos_gas_evolution.calc_SFR_MS_Genzel2015(z=3.0, lgMstar=10.5)

Computing Metal Abundance
-------------------------

Compute metal abundance from the Fundamental Metallicity Relation (FMR):

.. code-block:: python

    import a3cosmos_gas_evolution

    # Using Mannucci et al. 2010
    metalZ = a3cosmos_gas_evolution.calc_metalZ_from_FMR_following_Mannucci2010(
        M_star=1e10, SFR=1.0
    )

Computing Alpha_CO
------------------

Compute the CO-to-H2 conversion factor (alpha_CO):

.. code-block:: python

    import a3cosmos_gas_evolution

    # Using Genzel et al. 2015
    alpha_CO = a3cosmos_gas_evolution.calc_alphaCO_from_metalZ_following_Genzel2015a(
        metalZ=8.67
    )

Computing Dust Properties
-------------------------

Compute dust opacity and blackbody radiation:

.. code-block:: python

    from a3cosmos_gas_evolution.Common_Python_Code import (
        calc_dust_opacity_kappa,
        calc_dust_planck_law,
        calc_dust_rayleigh_jeans_law
    )

    kappa = calc_dust_opacity_kappa(850.0)  # wavelength in microns
    flux = calc_dust_planck_law(wavelength_um=100.0, T_dust=25.0)


Gas Mass from Dust
------------------

These functions calculate molecular gas mass from observed dust continuum or dust mass.

**calc_gas_mass_from_dust_continuum**

Calculate molecular gas mass from observed dust continuum flux.

Parameters:
    - ``obs_wavelength_um``: Observed wavelength in microns
    - ``obs_flux_mJy``: Observed flux in mJy
    - ``z``: Redshift
    - ``method``: Calibration method ('Hughes2017', 'Groves2015', 'Scoville2017')

Returns:
    - Molecular gas mass in solar masses

**calc_gas_mass_from_dust_mass**

Calculate molecular gas mass from dust mass using a gas-to-dust ratio (GDR).

Parameters:
    - ``M_dust``: Dust mass in solar masses
    - ``M_star``: Stellar mass (optional, used to estimate metallicity)
    - ``SFR``: Star formation rate (optional, used to estimate metallicity)
    - ``metallicity``: Gas metallicity 12+log(O/H) (optional)
    - ``GDR``: Gas-to-dust ratio (optional, computed from metallicity if not provided)
    - ``method``: 'Magdis2012', 'Genzel2015a', 'Genzel2015b', or 'GDR'

Returns:
    - Molecular gas mass in solar masses

Usage example:

.. code-block:: python

    from a3cosmos_gas_evolution import calc_gas_mass_from_dust_mass

    # Method 1: Using dust mass and metallicity
    M_mol = calc_gas_mass_from_dust_mass(
        M_dust=1e8,
        metallicity=8.67,
        GDR=100.0
    )

    # Method 2: Using stellar mass, SFR and redshift (auto-computes metallicity and GDR)
    M_mol = calc_gas_mass_from_dust_mass(
        M_dust=1e8,
        M_star=1e10,
        SFR=10.0,
        z=2.0,
        method='Genzel2015a'
    )


Gas Mass from Line
------------------

Calculate molecular gas mass from spectral line observations (CO or [CII]).

**calculate_gas_mass_from_line_flux**

Calculate (molecular) gas mass from CO or [CII] line emission.

Parameters:
    - ``z``: Redshift
    - ``lgMstar``: Log stellar mass (optional, used for alpha_CO estimation)
    - ``freq``: Observed frequency in GHz
    - ``flux``: Line flux in Jy km/s
    - ``line_name``: Line name ('CO(1-0)', 'CO(2-1)', '[CII]', etc.)
    - ``alpha_CO``: CO-to-H2 conversion factor (optional, computed from mass-metallicity if not provided)
    - ``alpha_CII``: [CII]-to-gas conversion factor (default: 30)

Returns:
    - ``lgMgas``: Log gas gas mass in solar masses
    - ``lgLline``: Log line luminosity
    - ``line_unit``: Line luminosity unit ('L_sun' or 'K km s-1 pc2')
    - ``line_type``: Line type ('CII' or 'CO')

Usage example:

.. code-block:: python

    from a3cosmos_gas_evolution import calculate_gas_mass_from_line_flux

    # From CO line with default mass-metallicity-dependent alpha_CO
    lgMgas, lgLline, unit, ltype = calculate_gas_mass_from_line_flux(
        z=2.0,
        lgMstar=10.5,
        freq=345.796/(1+2.0),  # CO(3-2) observed frequency GHz
        flux=1.0,     # Jy km/s
        line_name='CO(3-2)'
    )

    # From CO line with custom alphaCO
    lgMgas, lgLline, unit, ltype = calculate_gas_mass_from_line_flux(
        z=2.0,
        lgMstar=10.5,
        freq=345.796/(1+2.0),  # CO(3-2) observed frequency GHz
        flux=1.0,     # Jy km/s
        line_name='CO(3-2)',
        alpha_CO=0.8
    )

    # From CO line with custom alphaCO and excitation ratio
    lgMgas, lgLline, unit, ltype = calculate_gas_mass_from_line_flux(
        z=2.0,
        lgMstar=10.5,
        freq=345.796/(1+2.0),  # CO(3-2) observed frequency GHz
        flux=1.0,     # Jy km/s
        line_name='CO(3-2)',
        line_excitation_ratio=0.7,
        alpha_CO=0.8,
    )

    # From [CII] line with default alpha_CII = 30
    lgMgas, lgLline, unit, ltype = calculate_gas_mass_from_line_flux(
        z=6.0,
        lgMstar=10.5,
        freq=1900.0/(1+6.0),  # [CII] observed frequency GHz
        flux=1.0,     # Jy km/s
        line_name='[CII]',
        alpha_CII=30.0
    )

