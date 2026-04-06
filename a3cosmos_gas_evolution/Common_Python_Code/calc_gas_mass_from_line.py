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



# ============================================================================
# MOLECULAR GAS MASS CALCULATIONS
# ============================================================================

def solar_mass():
    """Solar mass in kg."""
    return 1.98847e30  # kg

def alpha_CO_zmass(z, M_star_log, metallicity=None):
    """
    Calculate α_CO from mass-metallicity relation.
    
    Based on:
    - Mass-metallicity relation: 12 + log(O/H) = f(M_*)
    - α_CO scaling with metallicity
    
    Args:
        z: Redshift
        M_star_log: Log stellar mass (log10(M_*/M_sun))
        metallicity: Optional explicit metallicity
    
    Returns:
        alpha_CO in M_sun/(K km/s pc^2)
    """
    # Mass-metallicity relation (Tremonti et al. 2004, modified)
    # 12 + log(O/H) = 9.0 - 0.35 * log(M_*/M_sun) + 0.10 * log^2(M_*/M_sun)
    # At high-z, metallicity may be lower
    
    # At z ~ 4, metallicity ~ 0.2-0.5 solar
    # α_CO ∝ Z^-1 (lower metallicity → higher α_CO)
    
    # Simplified: α_CO increases at high-z due to lower metallicity
    # α_CO(z=0) ~ 4.36 (Milky Way)
    # α_CO(z~4) ~ 3-5 (depends on galaxy type)
    
    if metallicity is None:
        # Estimate metallicity from mass (scaled by redshift)
        log_Z_solar = 0.3 - 0.3 * (1 - 1/(1+z))  # ~0.3 at z=4
    else:
        log_Z_solar = metallicity
    
    # α_CO scaling with metallicity
    Z_solar = 10**log_Z_solar
    alpha_mw = 4.36  # Milky Way value
    alpha_CO_val = alpha_mw / max(Z_solar, 0.1)  # Don't go too low
    
    # Cap at reasonable values
    alpha_CO_val = np.clip(alpha_CO_val, 0.8, 10.0)
    
    return alpha_CO_val

def calc_molecular_gas_mass(CO_luminosity, alpha_CO):
    """
    Calculate molecular gas mass from CO luminosity.
    
    M_H2 = α_CO * L'_CO
    
    Args:
        CO_luminosity: L'_CO in K km/s pc^2
        alpha_CO: α_CO in M_sun/(K km/s pc^2)
    
    Returns:
        M_H2 in solar masses
    """
    return alpha_CO * CO_luminosity


# 转换因子
ALPHA_CII = 30.0  # L_sun/M_sun ([CII] 光度 - 气体质量转换因子，文献范围 30-100)
ALPHA_CO_MW = 4.3  # M_sun/(K km/s pc²) (银河系)
ALPHA_CO_SB = 0.8  # M_sun/(K km/s pc²) (星暴/高红移)

def calculate_gas_mass_from_line(
        z, 
        logMstar, 
        freq, 
        flux, 
        line_name, 
        cosmo = cosmo,
        alpha_CO = None,
        alpha_CII = None, 
    ):
    """Calculate (molecular) gas mass

    Return: (M_gas, L_line, line_unit, line_name)
      - M_gas: (molecular) gas mass (M_sun)
      - L_line: line luminosity ('L_sun' or 'K km s-1 pc2')
      - line_unit: line luminosity unit ('L_sun' or 'K km s-1 pc2')
      - line_type: line type ('CII' or 'CO')
    """
    logMgas = np.nan
    logLline = np.nan
    line_unit = 'unknown'
    line_type = 'unknown'
    
    if z is None or str(z) == '' or np.isnan(float(z)) or float(z) <= 0:
        return logMgas, logLline, line_unit, line_type
    
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
        
        logMgas = np.log10(M_gas)
        logLline = np.log10(L_CII_Lsun)
        line_unit = 'L_sun'
        line_type = 'CII'
    
    elif line_name.startswith('CO'):
        # CO 光度 (K km/s pc²)
        # L'_CO = 3.25×10^7 × S_CO × ν_obs^{-2} × D_L² × (1+z)^{-3}
        # S_CO in Jy km/s, ν_obs in GHz, D_L in Mpc
        
        L_prime_CO = 3.25e7 * flux * (freq)**(-2) * (D_L)**2 * (1+z)**(-3)  # K km/s pc²
        
        # 气体质量 (使用高红移质量依赖的 α_CO)
        if alpha_CO is None: 
            if logMstar is None:
                alpha_CO = ALPHA_CO_SB
            else:
                alpha_CO = alpha_CO_zmass(z, logMstar)
        M_gas = alpha_CO * L_prime_CO  # M_sun
        
        logMgas = np.log10(M_gas)
        logLline = np.log10(L_prime_CO)
        line_unit = 'K km s-1 pc2'
        line_type = 'CO'
    
    return logMgas, logLline, line_unit, line_type







