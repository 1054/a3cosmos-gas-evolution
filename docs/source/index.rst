a3cosmos-gas-evolution documentation
====================================

A Python Package for Galaxy Cold Molecular Gas Mass and Depletion Time Calculation.

This Python package provides functions to calculate a galaxy's cold molecular gas mass to stellar mass ratio (gas fraction), cold molecular gas depletion time and galaxy main-sequence star formation rate.

The motivation is that the evolution of star-forming galaxies' star formation rate (SFR) and cold molecular gas reservoir have now been reasonably well measured out to very high redshift (z~6), from present time up to as early as one giga-year after the Big Bang (Madau & Dickinson 2014; Genzel et al. 2015; Scoville et al. 2016, 2017; Tacconi et al. 2018; Liu et al. 2018, 2019).

These studies have found that the majority of galaxies have a steady and parametrizable evolution in their stellar mass growth, SFR, and molecular gas mass (or molecular gas to total baryon fraction, i.e., gas fraction).

These evolution functions have provided crucial constraints to cosmological simulations of dark matter halo evolution and the semi-analytic modeling of the simulated galaxy evolution in the dark matter halo (e.g., Popping et al. 2014ab, 2016, 2017, 2019ab).

However, currently there are many parametrizations (or we say "equations") in the literature and each has its own limitation which is not very well aware by the generic users. Therefore, we provide this Python package which contains as many galaxy gas, dust, star formation and stellar mass evolution equations as possible for easier comparison and study.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   api
   modules

Quick Start
-----------

Installation:

.. code-block:: bash

    pip install a3cosmos-gas-evolution

Usage:

.. code-block:: python

    import a3cosmos_gas_evolution
    a3cosmos_gas_evolution.help()

    # Compute gas fraction
    f_gas = a3cosmos_gas_evolution.calc_gas_fraction_A3COSMOS(z=3.0, lgMstar=10.5, DeltaMS=0.5)

    # Compute gas depletion time
    tau_depl = a3cosmos_gas_evolution.calc_gas_depletion_time_A3COSMOS(z=3.0, lgMstar=10.5, DeltaMS=0.0)

    # Compute galaxy main-sequence SFR
    SFR_MS = a3cosmos_gas_evolution.calc_SFR_MS_Speagle2014(z=3.0, lgMstar=10.5)

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Acknowledgement
===============

Please cite `Liu et al. (2019b) <https://arxiv.org/abs/1910.12883>`__
and/or `ASCL <https://ascl.net/1910.003>`__ if you would like to use this package for your research.
