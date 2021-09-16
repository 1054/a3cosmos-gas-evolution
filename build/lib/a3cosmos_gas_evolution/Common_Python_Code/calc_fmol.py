#!/usr/bin/env python
# 
# Usage:
#    from calc_fmol import ( calc_fmol_from_metalZ_following_Krumholz2009, 
#                            calc_fmol_from_metalZ_following_Dave2016, 
#                            calc_fmol_from_metalZ_following_Popping2014 )
# 
# 

from __future__ import print_function

import os, sys, re, json, time, astropy
import numpy as np
from astropy.table import Table, Column, hstack
from copy import copy

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.725)

if sys.version_info.major >= 3:
    long = int
else:
    pass



# 
# def 
# 
def calc_fmol_from_metalZ_following_Krumholz2009(metalZOH):
    # 
    # compute fmol = M_H2/(M_HI+M_H2) following
    # -- Krumholz 2009 ApJ 699 850
    # -- Eq.(2)
    # 
    metalZ_solar = 8.69
    metalZ_ = 10**(metalZOH - metalZ_solar) # in unit of metalZ_solar
    Sigma_total_gas = 30.0 # 70.0 #<TODO># Schruba 2011 Fig.13
    ksi_KMT09 = 0.77 * (1. + 3.1 * np.power(metalZ_, 0.365))
    s_ = np.log(1. + 0.6 * ksi_KMT09) / (0.04 * Sigma_total_gas * metalZ_)
    delta_ = 0.0712 * (0.1 * s_**(-1) + 0.675)**(-2.8)
    f_H2_to_total = 1.0 - np.power(1.0 + (3./4. * (s_)/(1.+delta_))**(-5), -1./5.)
    return f_H2_to_total


def calc_fmol_from_metalZ_following_Dave2016(metalZOH):
    # 
    # compute fmol = M_H2/(M_HI+M_H2) following
    # -- Dave 2016 - cosmological hydrodynamic simulation MUFASA - 2016MNRAS.462.3265D.pdf
    # -- Eq.(1,2)
    # 
    metalZ_solar = 8.69
    metalZ_ = 10**(metalZOH - metalZ_solar)
    Sigma_total_gas = 30.0 # 70.0 #<TODO># Schruba 2011 Fig.13
    ksi_KMT09 = 0.77 * (1. + 3.1 * np.power(metalZ_, 0.365))
    s_ = np.log(1. + 0.6 * ksi_KMT09 + 0.01 * ksi_KMT09**2) / (0.0396 * metalZ_ * Sigma_total_gas)
    f_H2_to_total = 1.0 - 0.75 * s_ / (1.0 + 0.25 * s_)
    return f_H2_to_total


def calc_fmol_from_metalZ_following_Popping2014(metalZOH, U_MW = None):
    # 
    # compute fmol = M_H2/(M_HI+M_H2) following
    # -- http://www.astro.rug.nl/~gpopp/popping.samgas.pdf
    # -- Eq.(8)
    # 
    metalZ_solar = 8.69
    if U_MW is None:
        U_MW = 10**(5.0 * (metalZOH - metalZ_solar)) # UV background in units of MW UV ISRF <TODO>
        U_MW[U_MW<1.0] = 1.0
    if np.isscalar(U_MW):
        U_MW = np.array([U_MW]*len(metalZOH))
    D_MW = 10**(metalZOH - metalZ_solar)
    D_star = 1.5e-3 * np.log(1. + np.power(3.*U_MW,1.7))
    alpha_ = 5. * (U_MW/2.) / (1. + (U_MW/2.)**2)
    s_ = 0.04 / (D_star + D_MW)
    g_ = (1. + alpha_ * s_ + s_**2) / (1. + s_)
    Gamma_ = np.log(1. + g_ * np.power(D_MW,3./7.) * np.power(U_MW/15.,4./7.))
    Sigma_tide = 20. * np.power(Gamma_,4./7.) / D_MW / np.sqrt(1.+U_MW*D_MW**2)
    Sigma_total_gas_r = 10.0 # 70.0 #<TODO># Schruba 2011 Fig.13
    f_H2_r = (1. + Sigma_tide / (Sigma_total_gas_r))**(-2)
    #f_H2_to_HI = f_H2_r
    #f_H2_to_total = 1./(1.+1./f_H2_to_HI)
    f_H2_to_total = f_H2_r
    print('-----')
    for i in range(len(f_H2_to_total)):
        print('metalZ 12+log10(O/H) %0.3f, f_H2_to_total %0.6f, D_star %0.6f, D_MW %0.6f, U_MW %0.6f'%(metalZOH[i], f_H2_to_total[i], D_star[i], D_MW[i], U_MW[i]))
    return f_H2_to_total








