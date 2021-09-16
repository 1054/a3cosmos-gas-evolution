######################
a3cosmos-gas-evolution
######################

***********************************************************************************
A Python Package for Galaxy Cold Molecular Gas Mass and Depletion Time Calculation.
***********************************************************************************




A simple introduction:
======================

This Python package provides functions to calculate a galaxy's cold molecular gas mass to stellar mass ratio (gas fraction), cold molecular gas depletion time and galaxy main-sequence star formation rate. 

The motivation is that the evolution of star-forming galaxies' star formation rate (SFR) and cold molecular gas reservoir have now been reasonably well measured out to very high redshift (z~6), from present time up to as early as one giga-year after the Big Bang (Madau & Dickinson 2014; Genzel et al. 2015; Scoville et al. 2016, 2017; Tacconi et al. 2018; Liu et al. 2018, 2019). 
These studies have found that the majority of galaxies have a steady and parametrizable evolution in their stellar mass growth, SFR, and molecular gas mass (or molecular gas to total baryon fraction, i.e., gas fraction). 
These evolution functions have provided crucial constraints to cosmological simulations of dark matter halo evolution and the semi-analytic modeling of the simulated galaxy evolution in the dark matter halo (e.g., Popping et al. 2014ab, 2016, 2017, 2019ab). 
However, currently there are many parametrizations (or we say "equations") in the literature and each has its own limitation which is not very well aware by the generic users. 
Therefore, we provide this Python package which contains as many galaxy gas, dust, star formation and stellar mass evolution equations as possible for easier comparison and study. 




A simple usage:
===============

To get cold molecular gas density evolution curve, e.g., Fig. 15 of `D. Liu et al. (2019b) <https://arxiv.org/abs/1910.12883>`__

.. code-block:: python

    import a3cosmos_gas_evolution
    a3cosmos_gas_evolution.help()
    z, rho_mol_gas = a3cosmos_gas_evolution.get_cosmic_mol_gas_density_A3COSMOS() # return rho_mol_gas in solar mass per cubic mega parsec.
    z, rho_mol_gas = a3cosmos_gas_evolution.get_cosmic_mol_gas_density_Tacconi2018() # return rho_mol_gas in solar mass per cubic mega parsec.
    z, rho_mol_gas = a3cosmos_gas_evolution.get_cosmic_mol_gas_density_Scoville2017() # return rho_mol_gas in solar mass per cubic mega parsec.
    a3cosmos_gas_evolution.plot_cosmic_mol_gas_density() # or plot it with matplotlib



To compute gas fraction, i.e., ``M_molgas / (M_molgas + M_star)``

.. code-block:: python

    import a3cosmos_gas_evolution
    a3cosmos_gas_evolution.help()
    a3cosmos_gas_evolution.calc_gas_fraction_A3COSMOS(z = 3.0, lgMstar = 10.5, DeltaMS = 0.5) # or we can input cosmic_age = 2.178 instead of z = 3.0
    a3cosmos_gas_evolution.calc_gas_fraction_Tacconi2018(z = 3.0, lgMstar = 10.5, DeltaMS = 0.5)
    a3cosmos_gas_evolution.calc_gas_fraction_Scoville2017(z = 3.0, lgMstar = 10.5, DeltaMS = 0.5)



To compute gas depletion time, i.e., ``M_molgas / SFR``

.. code-block:: python

    import a3cosmos_gas_evolution
    a3cosmos_gas_evolution.help()
    a3cosmos_gas_evolution.calc_gas_depletion_time_A3COSMOS(z = 3.0, lgMstar = 10.5, DeltaMS = 0.0)
    a3cosmos_gas_evolution.calc_gas_depletion_time_Tacconi2018(z = 3.0, lgMstar = 10.5, DeltaMS = 0.0)
    a3cosmos_gas_evolution.calc_gas_depletion_time_Scoville2017(z = 3.0, lgMstar = 10.5, DeltaMS = 0.0)



To compute galaxy main-sequence star formation rate, i.e., ``SFR_MS``

.. code-block:: python

    import a3cosmos_gas_evolution
    a3cosmos_gas_evolution.help()
    a3cosmos_gas_evolution.calc_SFR_MS_Speagle2014(z = 3.0, lgMstar = 10.5)




Acknowledgement:
================
Please cite `Liu et al. (2019b) <https://arxiv.org/abs/1910.12883>`__ 
and/or 
`https://ascl.net/1910.003 <https://ascl.net/1910.003>`__
if you would like to use this package for your research. Thank you. 






