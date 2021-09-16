#!/usr/bin/env python
# 

from __future__ import print_function

import os, sys, re, json, time, astropy
import numpy as np
from astropy.table import Table, Column, hstack
from copy import copy
from numpy import log10, power as pow

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=73, Om0=0.27, Tcmb0=2.725)

if sys.version_info.major >= 3:
    long = int
else:
    pass



# 
# def 
# 
def convert_age_to_z(cosmoAge):
    # reversedly compute z from age
    if type(cosmoAge) is list: cosmoAge = np.array(cosmoAge)
    spl_opz_log10 = np.linspace(np.log10(1.0+0.0), np.log10(1.0+9999.0), num=1000, endpoint=True) # for converting cosmoAge to opz (1+z)
    spl_z = 10**spl_opz_log10-1.0
    spl_cosmoAge = cosmo.age(spl_z).value
    spl_cosmoAge_log10 = np.log10(spl_cosmoAge)
    spl_cosmoAge_log10 = spl_cosmoAge_log10[::-1]
    spl_opz_log10 = spl_opz_log10[::-1]
    tmp_opz_log10 = np.interp(np.log10(cosmoAge), spl_cosmoAge_log10, spl_opz_log10)
    z = 10**tmp_opz_log10-1.0
    return z


# CSFRD Madau & Dickinson 2014
def calc_cosmic_star_formation_rate_density_MadauDickinson2014(z):
    # Madau & Dickinson (2014)
    # converted to Chabrier IMF from Salpeter IMF by a factor of 1.64
    rho_SFR = 0.015 * (1+z)**2.7 / (1.0 + ((1+z)/2.9)**5.6) / 1.64
    return rho_SFR

def calc_CSFRD_Madau2014(z):
    return calc_cosmic_star_formation_rate_density_MadauDickinson2014(z)


# CSFRD Liudz 2018
def calc_cosmic_star_formation_rate_density_Liu2018(z, shape = 'double-powerlaw'):
    # D. Liu et al. 2018 (bibcode:2018ApJ...853..172L)
    # Chabrier IMF
    # shape can be 'double-powerlaw' or 'log-normal'
    if shape.startswith('double'):
        rho_SFR = 0.00587 * (1+z)**3.0 / (1.0 + ((1+z)/2.9)**5.6)
    elif shape.startswith('log'):
        A0 = 0.575
        tau = 0.66
        T0 = 1.50
        t = cosmo.age(z).value
        rho_SFR = A0 / (t * np.sqrt(2.0*np.pi*(tau**2))) * np.exp((-(np.log(t)-T0)**2) / (2.0*(tau**2)))
    else:
        raise ValueError('Error! The input shape %s is not allowed by the called function calc_cosmic_star_formation_rate_density_Liu2018()!'%(shape))
    return rho_SFR

def calc_CSFRD_Liu2018(z, shape = 'double-powerlaw'):
    return calc_cosmic_star_formation_rate_density_Liu2018(z, shape = shape)






