#!/usr/bin/env python
# 
# 20190222
#     copied from "calc_stellar_mass_function.py", this code will superceed "calc_stellar_mass_function.py". 
# 

from __future__ import print_function

import os, sys, re, json, time, astropy
import numpy as np
from astropy.table import Table, Column, hstack
from copy import copy
from numpy import log, log10, power, sum, sqrt, pi, exp
pow = power
lg = log10
ln = log
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d

if not (os.path.dirname(os.path.abspath(__file__)) in sys.path): sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import apply_cosmology
cosmo = apply_cosmology.cosmo

if sys.version_info.major >= 3:
    long = int
else:
    pass




# 
# def 
# 
def Schechter_Function_for_LF(L, L_character, Phi_character, alpha):
    # 
    # Schechter (1976)
    # 
    # Phi(L) dL = (Phi_*) * (L/L_*)**(alpha) * exp(-L/L_*) dL/L_*
    #           = (Phi_*) * x**(alpha) * exp(-x) dx
    #           = (Phi_*) * 10**(lgx * alpha) * exp(-10**lgx) dx
    #           = (Phi_*) * 10**(lgx * alpha) * exp(-10**lgx) dlnx
    #           = (Phi_*) * 10**(lgx * alpha) * exp(-10**lgx) dlgx * ln(10)
    #           = (Phi_*) * 10**((lgL-lgL_*)*(alpha+1)) * exp(-10**(lgL-lgL_*)) * ln(10)  dlgx
    #           = (Our_Phi_Phi_Schechter)                                                 dlgx
    # 
    #lgx = lgL-lg_L0
    #Phi_Schechter = phi * (10**(lgx*(alpha+1))) * (np.exp(-10**lgx)) * ln(10) # per dex and already multiplied ln(10), so that its integral directly equals \int Phi(L) / L dL
    # 
    Phi_Schechter = Phi_character * (L/L_character)**(alpha) * np.exp(-(L/L_character)) # Mpc-3 dex-1
    #Phi_Schechter = Phi_Schechter * ln(10)
    return Phi_Schechter


def Saunders_Function_for_LF(L, L_character, Phi_character, alpha, sigma):
    # Saunders et al. (1990)
    Phi_Saunders = Phi_character * (L/L_character)**(1-alpha) * np.exp(-1.0/(2.0*sigma**2) * (np.log10(1.0+(L/L_character)))**2 )
    #print('Phi_character', Phi_character)
    #print('(L/L_character)**(1-alpha)', (L/L_character)**(1-alpha))
    #print('np.exp(-1.0/(2.0*sigma**2) * (np.log10(1.0+(L/L_character)))**2 )', np.exp(-1.0/(2.0*sigma**2) * (np.log10(1.0+(L/L_character)))**2 ))
    #print('Phi_Saunders', Phi_Saunders)
    return Phi_Saunders





# 
# def 
# 
def calc_radio_LF_Novak2017(z, lgL=None, galaxy_type = 'SFG'):
    # 
    # Novak 2017 bibcode:2017A&A...602A...5N
    # IMF: Chabrier 2003
    # Saunders et al. (1990)
    # Outputs: lgL_grid, lgPhi_grid
    # 
    # check z
    if not np.isscalar(z):
        raise ValueError('Please input a float number as the redshift!')
    # 
    # check galaxy_type
    if not (type(galaxy_type) is str):
        raise ValueError('Please input either "ALL", "SFG" or "QG" as the galaxy_type!')
    else:
        if not (galaxy_type in ['ALL', 'SFG', 'QG']):
            raise ValueError('Please input either "ALL", "SFG" or "QG" as the galaxy_type!')
    # 
    # make lgL_grid
    if lgL is None:
        lgL_grid = np.linspace(18.0, 25.0, num=1000, endpoint=True)
    else:
        lgL_grid = lgL
    # 
    L_grid = 10**lgL_grid
    # 
    # read LF parameters
    L_character = 1.85e21 # * 1.4e9 / 3.839e25 # vLv(1.4GHz,rest) = W Hz-1 --> Lsun
    Phi_character = 3.55e-3 # Mpc-3 dex-1
    alpha = 1.22
    sigma = 0.63
    # 
    #Phi_z0 = Saunders_Function(L_grid, L_character, Phi_character, alpha, sigma)
    # 
    # check z
    LF_zmin = 0.0
    LF_zmax = +np.inf
    if z < LF_zmin or z > LF_zmax:
        raise ValueError('calc_radio_LF_Novak2017: The input redshift is out of the allowed range of %s -- %s!'%(LF_zmin, LF_zmax))
    # 
    # scale to z via pure luminosity evolution
    alphaL = 3.16
    betaL = -0.32
    L_grid_z = (L_grid / ((1.0+z)**(alphaL+(z*betaL))))
    Phi = Saunders_Function_for_LF(L_grid_z, L_character, Phi_character, alpha, sigma)
    lgPhi = np.log10(Phi)
    # 
    if lgL is None:
        return lgL_grid, lgPhi
    else:
        return lgPhi



def calc_IR_250um_LF_Koprowski2017(z, lgL=None, galaxy_type = 'SFG'):
    # 
    # Koprowski 2017 bibcode:2017MNRAS.471.4155K
    # IMF: Chabrier 2003
    # Saunders et al. (1990)
    # Outputs: lgL_grid, lgPhi_grid
    # 
    # check z
    if not np.isscalar(z):
        raise ValueError('Please input a float number as the redshift!')
    # 
    # check galaxy_type
    if not (type(galaxy_type) is str):
        raise ValueError('Please input either "ALL", "SFG" or "QG" as the galaxy_type!')
    else:
        if not (galaxy_type in ['ALL', 'SFG', 'QG']):
            raise ValueError('Please input either "ALL", "SFG" or "QG" as the galaxy_type!')
    # 
    # make lgL_grid
    if lgL is None:
        lgL_grid = np.linspace(24.0, 27.0, num=1000, endpoint=True)
    else:
        lgL_grid = lgL
    # 
    L_grid = 10**lgL_grid
    # 
    # read LF parameters
    table_z_lower = [0.5, 1.5, 2.5, 3.5]
    table_z_upper = [1.5, 2.5, 3.5, 4.5]
    table_lgL_character = [25.20, 25.40, 25.63, 25.84] # W Hz-1
    table_lgPhi_character = [-2.88, -3.03, -3.73, -4.59] # Mpc-3 dex-1
    alpha = -0.4
    # 
    # check z
    LF_zmin = table_z_lower[0]
    LF_zmax = table_z_upper[-1]
    if z < LF_zmin or z > LF_zmax:
        raise ValueError('calc_IR_250um_LF_Koprowski2017: The input redshift is out of the allowed range of %s -- %s!'%(LF_zmin, LF_zmax))
    # 
    # scale to z (using step function... <TODO>)
    Phi = None
    lgPhi = None
    for i in range(len(table_z_upper)):
        if z >= table_z_lower[i] and z <= table_z_upper[i]:
            L_character = 10**(table_lgL_character[i])
            Phi_character = 10**(table_lgPhi_character[i])
            Phi = Schechter_Function_for_LF(L_grid, L_character, Phi_character, alpha)
            lgPhi = np.log10(Phi)
            break
    # 
    if lgL is None:
        return lgL_grid, lgPhi
    else:
        return lgPhi



def calc_IR_LF_Gruppioni2013(z, lgL=None, galaxy_type = 'SFG'):
    # 
    # Gruppioni 2013 bibcode:
    # IMF: Chabrier 2003
    # H0 = 71 km s−1 Mpc−1, Ωm = 0.27, and ΩΛ = 0.73.
    # Saunders et al. (1990)
    # Outputs: lgL_grid, lgPhi_grid
    # 
    # check z
    if not np.isscalar(z):
        raise ValueError('Please input a float number as the redshift!')
    # 
    # check galaxy_type
    if not (type(galaxy_type) is str):
        raise ValueError('Please input either "ALL", "SFG" or "QG" as the galaxy_type!')
    else:
        if not (galaxy_type in ['ALL', 'SFG', 'QG']):
            raise ValueError('Please input either "ALL", "SFG" or "QG" as the galaxy_type!')
    # 
    # make lgL_grid
    if lgL is None:
        lgL_grid = np.linspace(8.0, 14.0, num=1000, endpoint=True)
    else:
        lgL_grid = lgL
    # 
    L_grid = 10**lgL_grid
    # 
    # read LF parameters (their Table 7)
    table_data = [ [0.0 , 0.3 , 1.15, 0.52, 10.12, -2.29], 
                   [0.3 , 0.45, 1.2 , 0.5 , 10.41, -2.31], 
                   [0.45, 0.6 , 1.2 , 0.5 , 10.55, -2.35], 
                   [0.6 , 0.8 , 1.2 , 0.5 , 10.71, -2.35], 
                   [0.8 , 1.0 , 1.2 , 0.5 , 10.97, -2.40], 
                   [1.0 , 1.2 , 1.2 , 0.5 , 11.13, -2.43], 
                   [1.2 , 1.7 , 1.2 , 0.5 , 11.37, -2.70], 
                   [1.7 , 2.0 , 1.2 , 0.5 , 11.50, -3.00], 
                   [2.0 , 2.5 , 1.2 , 0.5 , 11.60, -3.01], 
                   [2.5 , 3.0 , 1.2 , 0.5 , 11.92, -3.27], 
                   [3.0 , 4.2 , 1.2 , 0.5 , 11.90, -3.74] ]
    table_data = np.array(table_data).T
    table_z_lower = table_data[0]
    table_z_upper = table_data[1]
    table_alpha = table_data[2]
    table_sigma = table_data[3]
    table_lgL_character = table_data[4] # Lsun
    table_lgPhi_character = table_data[5] # Mpc-3 dex-1
    # 
    # check z
    LF_zmin = table_z_lower[0]
    LF_zmax = table_z_upper[-1]
    if z < LF_zmin or z > LF_zmax:
        raise ValueError('calc_IR_LF_Gruppioni2013: The input redshift is out of the allowed range of %s -- %s!'%(LF_zmin, LF_zmax))
    # 
    # scale to z (using step function... <TODO>)
    Phi = None
    lgPhi = None
    for i in range(len(table_z_upper)):
        if z >= table_z_lower[i] and z <= table_z_upper[i]:
            L_character = 10**(table_lgL_character[i])
            Phi_character = 10**(table_lgPhi_character[i])
            alpha = table_alpha[i]
            sigma = table_sigma[i]
            Phi = Saunders_Function_for_LF(L_grid, L_character, Phi_character, alpha, sigma)
            lgPhi = np.log10(Phi)
            break
    # 
    if lgL is None:
        return lgL_grid, lgPhi
    else:
        return lgPhi























