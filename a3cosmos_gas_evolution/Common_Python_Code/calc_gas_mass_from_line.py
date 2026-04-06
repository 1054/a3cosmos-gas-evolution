#!/usr/bin/env python
# 
# Usage:
#    from calc_gas_mass_from_line import calc_gas_mass_from_line
# 
# 

from __future__ import print_function

import os, sys, re, json, time, astropy
import numpy as np
from astropy import units as u
from astropy.table import Table, Column, hstack
from copy import copy

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.725)

if sys.version_info.major >= 3:
    long = int
else:
    pass

from calc_metal_Z import calc_metalZ_from_FMR_following_Genzel2015a
from calc_alpha_CO import calc_alphaCO_from_metalZ_following_Genzel2015a


# ============================================================================
# MOLECULAR GAS MASS CALCULATIONS
# ============================================================================

# 转换因子
ALPHA_CII = 30.0  # L_sun/M_sun ([CII] 光度 - 气体质量转换因子，文献范围 30-100)
ALPHA_CO_MW = 4.3  # M_sun/(K km/s pc²) (银河系)
ALPHA_CO_SB = 0.8  # M_sun/(K km/s pc²) (星暴/高红移)

def calc_gas_mass_from_line_flux(
        z, 
        lgMstar, 
        freq, 
        flux, 
        line_name, 
        cosmo = cosmo,
        alpha_CO = None,
        alpha_CII = None, 
        line_excitation_ratio = None,
        verbose = False,
    ):
    """Calculate (molecular) gas mass

    Return: (log_M_gas, log_L_line, line_unit, line_name)
      - log_M_gas: (molecular) gas mass (M_sun)
      - log_L_line: line luminosity ('L_sun' or 'K km s-1 pc2')
      - line_unit: line luminosity unit ('L_sun' or 'K km s-1 pc2')
      - line_type: line type ('CII' or 'CO')
    """
    lgMgas = float(np.nan)
    lgLline = float(np.nan)
    line_unit = 'unknown'
    line_type = 'unknown'
    
    if z is None or str(z) == '' or np.isnan(float(z)) or float(z) <= 0 or np.isnan(float(flux)) or float(flux) <= 0:
        return lgMgas, lgLline, line_unit, line_type
    
    # 光度距离
    D_L = cosmo.luminosity_distance(z).to(u.Mpc).value  # Mpc
    
    if line_name.startswith('[CII]') or line_name.startswith('CII'):
        # [CII] 光度 (L_sun)
        # L = 4π × D_L² × F × hν / (1+z)
        # F (Jy km/s) → F (erg/s/cm²)
        # 1 Jy = 10^{-23} erg/s/cm²/Hz
        # F_nu = F / Δv (假设线宽)
        # 简化：L_[CII] (L_sun) ≈ 1.04×10^-3 × (D_L/Mpc)² × (F/Jy km/s) × (ν_rest/GHz) / (1+z)
        # 简化：L_[CII] (L_sun) ≈ 1.04×10^-3 × (D_L/Mpc)² × (F/Jy km/s) × (ν_obs/GHz)
        
        L_CII_Lsun = 1.04e-3 * (D_L)**2 * flux * freq
        
        # 气体质量
        if alpha_CII is None:
            alpha_CII = ALPHA_CII
        M_gas = L_CII_Lsun / alpha_CII  # M_sun
        
        lgMgas = np.log10(M_gas)
        lgLline = np.log10(L_CII_Lsun)
        line_unit = 'L_sun'
        line_type = 'CII'
    
    elif line_name.startswith('CO'):
        # CO 光度 (K km/s pc²)
        # L'_CO = 3.25×10^7 × S_CO × ν_obs^{-2} × D_L² × (1+z)^{-3}
        # S_CO in Jy km/s, ν_obs in GHz, D_L in Mpc
        
        L_prime_CO = 3.25e7 * flux * (freq)**(-2) * (D_L)**2 * (1+z)**(-3)  # K km/s pc²

        # 激发改正
        r_J1_map = {
            'CO(1-0)': 1.0,
            'CO(2-1)': 0.8,
            'CO(3-2)': 0.5,
            'CO(4-3)': 0.4,
            'CO(5-4)': 0.3,
            'CO(6-5)': 0.25,
            'CO(7-6)': 0.2,
            'CO(8-7)': 0.15,
        }
        if line_excitation_ratio is None:
            r_J1 = r_J1_map.get(line_name, 0.5)
        else:
            if line_excitation_ratio > 1.0:
                raise ValueError(f'line_excitation_ratio ({line_excitation_ratio}) must be <= 1.0')
            r_J1 = line_excitation_ratio
        L_prime_CO10 = L_prime_CO / r_J1
        
        # 气体质量 (使用高红移质量依赖的 α_CO)
        if alpha_CO is None: 
            if lgMstar is None:
                alpha_CO = ALPHA_CO_SB
            else:
                SFR = None # not needed
                metalZ = calc_metalZ_from_FMR_following_Genzel2015a(10**lgMstar, SFR, z)
                alpha_CO = calc_alphaCO_from_metalZ_following_Genzel2015a(metalZ)
                if verbose:
                    print(f'metalZ: {metalZ:.2f}, Z/Zsun: {10**(metalZ-8.69):.2f}, alpha_CO: {alpha_CO:.2f}')
        M_gas = alpha_CO * L_prime_CO10  # M_sun
        if verbose:
            print(f'alpha_CO: {alpha_CO:.2f}, L_prime_CO: {L_prime_CO:.3e}, L_prime_CO10: {L_prime_CO10:.3e}, M_gas: {M_gas:.3e}')
        
        lgMgas = float(np.log10(M_gas))
        lgLline = float(np.log10(L_prime_CO)) # noting that the return is L_prime_CO not L_prime_CO10
        line_unit = 'K km s-1 pc2'
        line_type = 'CO'
    
    return lgMgas, lgLline, line_unit, line_type


