#!/usr/bin/env python
# 
# 20190222
#     copied from "calc_main_sequence.py", this code will superceed "calc_main_sequence.py". 
# 

from __future__ import print_function

import os, sys, re, json, time, astropy
import numpy as np
from astropy.table import Table, Column, hstack
from copy import copy
from numpy import log10, power as pow

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
def calc_Sargent2014_sSFR(z, lgMstar=10.5, DeltaMS=0.0):
    return 0.095 * 10**(-0.21*(lgMstar-np.log10(5e10))) * np.exp(2.05*z/(1.0+0.16*z**1.54)) * 10**(DeltaMS)
def calc_Speagle2014_sSFR(cosmoAge, lgMstar=10.5, DeltaMS=0.0):
    return 10**((0.84 - 0.026*cosmoAge) * lgMstar - (6.51 - 0.11*cosmoAge)) / 10**(lgMstar) * 1e9 * 10**(DeltaMS)
def calc_Genzel2015_sSFR(z, lgMstar=10.5, DeltaMS=0.0):
    return 10**(-1.12 + 1.14 * z - 0.19 * z**2 - (0.3 + 0.13 * z) * (lgMstar - 10.5))
def calc_Scoville2017_sSFR(z, lgMstar=10.5, DeltaMS=0.0):
    lgMstar_ref = 10.5
    SFR_MS_ref = 10**(0.59*lgMstar_ref-5.77)* np.power((1.0+z),(0.22*lgMstar_ref+0.59))
    SFR_MS = SFR_MS_ref * 10**(1.72-np.log10(1+np.power(10**(lgMstar-10.31),-1.07))) / 10**(1.72-np.log10(1+np.power(10**(lgMstar_ref-10.31),-1.07)))
    sSFR_MS = SFR_MS / 10**(lgMstar) * 1e9 # Gyr
    return sSFR_MS
def calc_sSFR_MS(lgMstar, z, cosmoAge=None):
    if cosmoAge is None:
        cosmoAge = cosmo.age(z).value
    sSFR_MS = calc_Speagle2014_sSFR(cosmoAge, lgMstar)
    # 
    #sSFR_MS = calc_Sargent2014_sSFR(z, lgMstar)
    # 
    return sSFR_MS



def calc_SFR_MS_Scoville2017(z, lgMstar=10.5):
    # Scoville2017 used a composite model
    # They use the shape of MS from Lee2015, and the evolution of the MS from Speagle2014. 
    # They said that they use the Speagle2014 MS #49 fit with z as the variable, and explicitly described the equation in their Sect. 2.1:
    #   SFR_MS_Speagle2014 = 10**(0.59*lgMstar-5.77) * (1+z)**(0.22*lgMstar+0.59)
    # They use this function evaluated at log(M*) = 10.5 Msun, to describe the redshift evolution of the MS:
    #   SFR_MS_Speagle2014(z) = SFR_MS_Speagle2014(z=0) * (1+z)**(2.9)
    # Then at each redshift, they scale the Lee2015 MS shape which is explicitly described in their Sect. 2.2:
    #   SFR_MS_Lee2015(lgMstar) = 10**(1.72-log10(1+ (10**(lgMstar-10.31))**(-1.07) ))
    #   SFR_MS(z,lgMstar) = SFR_MS_Lee2015(lgMstar) / SFR_MS_Lee2015(lgMstar=10.5) * SFR_MS_Speagle2014(z) / SFR_MS_Speagle2014(z=0)
    lgMstar_ref = 10.5
    SFR_MS_ref = 10**(0.59*lgMstar_ref-5.77)* np.power((1.0+z),(0.22*lgMstar_ref+0.59))
    SFR_MS = SFR_MS_ref * 10**(1.72-np.log10(1+np.power(10**(lgMstar-10.31),-1.07))) / 10**(1.72-np.log10(1+np.power(10**(lgMstar_ref-10.31),-1.07)))
    sSFR_MS = SFR_MS / 10**(lgMstar) * 1e9 # Gyr
    return SFR_MS


def calc_SFR_MS_Tacconi2018(z, lgMstar=10.5, with_dzliu_t_z_conversion = False):
    # They used Speagle2014 time evolution equation, but applied a cosmic time to redshift conversion of:
    # t_cosmic_age = pow(10.0, 1.143 - 1.026 * log10(1+z) - 0.599 * (log10(1+z))**2 + 0.528 * (log10(1+z))**3)
    t_cosmic_age = pow(10.0, 1.143 - 1.026 * log10(1+z) - 0.599 * (log10(1+z))**2 + 0.528 * (log10(1+z))**3)
    if with_dzliu_t_z_conversion:
        # only valid for cosmo = FlatLambdaCDM(H0=73, Om0=0.27, Tcmb0=2.725)
        x = log10(1+z)
        y = (1.1295867233) + (-1.0876769081) * x + (-0.5392282897) * x**2 + (0.2310061623) * x**3
        t_cosmic_age = pow(10.0, y)
    # 
    sSFR_MS = calc_Speagle2014_sSFR(t_cosmic_age, lgMstar)
    SFR_MS = sSFR_MS / 1e9 * 10**(lgMstar)
    return SFR_MS


def calc_SFR_MS_Speagle2014(z, lgMstar=10.5, cosmic_age=None):
    # Speagle 2014 Table 7 #49, same as calc_SFR_MS_Speagle2014_with_t()
    if cosmic_age is None:
        cosmoAge = cosmo.age(z).value
    else:
        cosmoAge = cosmic_age
    # 
    SFR_MS = 10**((0.84 - 0.026*cosmoAge) * lgMstar - (6.51 - 0.11*cosmoAge))
    sSFR_MS = SFR_MS / 10**(lgMstar) * 1e9 # Gyr
    return SFR_MS


def calc_SFR_MS_Speagle2014_with_t(z, lgMstar=10.5, cosmic_age=None):
    # Speagle 2014 Table 7 #49
    if cosmic_age is None:
        cosmoAge = cosmo.age(z).value
    else:
        cosmoAge = cosmic_age
    # 
    cosmoAge = cosmo.age(z).value
    SFR_MS = 10**((0.84 - 0.026*cosmoAge) * lgMstar - (6.51 - 0.11*cosmoAge))
    sSFR_MS = SFR_MS / 10**(lgMstar) * 1e9 # Gyr
    return SFR_MS


def calc_SFR_MS_Speagle2014_best_fit_with_t_rewritten_by_dzliu(z=None, cosmic_age=None, lgMstar=10.5):
    # Speagle 2014 Table 7 #49
    if z is None and cosmic_age is None:
        print('Error! Please input z or cosmic_age when calling "calc_SFR_MS_Speagle2014_with_t_rewritten_by_dzliu()"!')
        sys.exit()
    # 
    if cosmic_age is None:
        cosmoAge = cosmo.age(z).value
    else:
        cosmoAge = cosmic_age
    # 
    cosmoAge = cosmo.age(z).value
    #SFR_MS = 10**((0.84 - 0.026*cosmoAge) * lgMstar - (6.51 - 0.11*cosmoAge))
    SFR_MS = 10**( 0.84 * (lgMstar-10.0) + (-0.15-0.026*(lgMstar-10.0))*cosmoAge + 1.89 )
    sSFR_MS = SFR_MS / 10**(lgMstar) * 1e9 # Gyr
    return SFR_MS


def calc_SFR_MS_Speagle2014_with_z(z, lgMstar=10.5, cosmic_age=None):
    # Speagle 2014 Table 8 #49
    # -- used by Scoville+2016,2017
    if cosmic_age is None:
        cosmoAge = cosmo.age(z).value
    else:
        cosmoAge = cosmic_age
    # 
    SFR_MS = 10**((0.22 * lgMstar + 0.59) * np.log10(1.+z) + 0.59 * lgMstar - 5.77)
    sSFR_MS = SFR_MS / 10**(lgMstar) * 1e9 # Gyr
    return SFR_MS


def calc_SFR_MS_Speagle2014_best_fit_with_z(z=None, lgMstar=10.5):
    # Speagle 2014 Table 8 #49
    # -- used by Scoville+2016,2017
    if z is None:
        print('Error! Please input z when calling "calc_SFR_MS_Speagle2014_best_fit_with_z()"!')
        sys.exit()
    # 
    SFR_MS = 10**((0.22 * lgMstar + 0.59) * np.log10(1.+z) + 0.59 * lgMstar - 5.77)
    sSFR_MS = SFR_MS / 10**(lgMstar) * 1e9 # Gyr
    return SFR_MS


def calc_SFR_MS_Sargent2014(z, lgMstar=10.5):
    # Sargent,2014ApJ...793...19S
    # Appendix, Eq.(A1)
    #   sSFR = N(M*) * exp(A*z / (1+B*z^C))
    #   where N(M⋆) = N(5 × 10^{10} M⊙) 10^{ν log(M⋆/[5×10^10 M⊙])}
    #     N(5 × 10^{10} M⊙) = 0.095^{+0.002}_{-0.003} Gyr^{-1}
    #     ν is the slope of the log(sSFR) versus log(M⋆) relation, ν = -0.21±0.04
    #     A = 2.05+0.33−0.20
    #     B = 0.16+0.15−0.07
    #     C = 1.54±0.32
    # 
    sSFR_MS = 0.095 * 10**(-0.21*(lgMstar-np.log10(5e10))) * np.exp(2.05*z/(1.0+0.16*z**1.54)) # Gyr
    SFR_MS = sSFR_MS * 10**(lgMstar) / 1e9
    return SFR_MS


def calc_SFR_MS_Bethermin2015(z, lgMstar=10.5):
    # see paper http://arxiv.org/pdf/1409.5796v2
    #     section 4.2
    #     paragraph 1 
    #     sSFR = (0.061±0.006Gyr**{−1}) * (1+z)**{2.82±0.12} at z<2 and as (1+z)**{2.2±0.3} at z>2.
    if np.isscalar(z):
        z = np.array([z])
    elif type(z) is list:
        z = np.array(z)
    mask = (z>=2.0)
    Bethermin2014_sSFR = (0.061 * (1.+z)**2.82)
    Bethermin2014_sSFR[mask] = ((1.+z[mask])**2.2)  / (1.+2)**2.2 * (0.061 * (1.+2)**2.82) # renormalize to z=2
    Bethermin2014_SFR = (Bethermin2014_sSFR / 1e9 * 10**(lgMstar))
    SFR_MS = Bethermin2014_SFR
    sSFR_MS = Bethermin2014_sSFR
    return SFR_MS


def calc_SFR_MS_Schreiber2015(z, lgMstar=10.5):
    # Schreiber 2015 - 2015A&A...575A..74S.pdf
    # see their Eq.(9)
    #   log10(SFR_MS[M⊙/yr]) = m − m0 + a0*r - a1 * max([0, m - m1 - a2*r]) ^ 2
    #   where r = log10(1 + z), m = log10(M∗/1e9 M⊙), and m0,m1,a0,a1,a2 are constants. 
    # Here we converts Chabrier IMF values to Salpeter IMF values for internal calculation
    if np.isscalar(z):
        z = np.array([z])
    elif type(z) is list:
        z = np.array(z)
    Schreiber2014_lgMstar = np.array(lgMstar)+np.log10(1.73) # converts Chabrier IMF values to Salpeter IMF values for internal calculation
    Schreiber2014_SFR = Schreiber2014_lgMstar-9.0-0.50+1.5*np.log10(1.+z)
    Schreiber2014_mmr = Schreiber2014_lgMstar-9.0-0.36-2.5*np.log10(1.+z)
    Schreiber2014_mask = (Schreiber2014_mmr>0)
    Schreiber2014_SFR[Schreiber2014_mask] = (Schreiber2014_SFR[Schreiber2014_mask]-0.30*Schreiber2014_mmr[Schreiber2014_mask]**2)
    Schreiber2014_SFR = np.power(10.0, Schreiber2014_SFR-np.log10(1.73)) # converts back to Chabrier IMF values
    SFR_MS = Schreiber2014_SFR
    sSFR_MS = SFR_MS / 10**(lgMstar) * 1e9
    return SFR_MS


def calc_SFR_MS_Genzel2015(z, lgMstar=10.5):
    # Genzel 2015ApJ...800...20G
    # The evolution of sSFR is taken from Whitaker+2012, see Genzel+2015 Eq.(1) on page 2, 
    #   log10(sSFR/[Gyr^{-1}]) = -1.12 + 1.14 * z - 0.19 * z**2 - (0.3 + 0.13 * z) * (lgMstar - 10.5)
    # assumed Chabrier2003 IMF
    # 
    # -- taken from Whitaker+2012
    sSFR_MS = 10**(-1.12 + 1.14 * z - 0.19 * z**2 - (0.3 + 0.13 * z) * (lgMstar - 10.5))
    SFR_MS = sSFR_MS * np.power(10.0, lgMstar)
    return SFR_MS


def calc_SFR_MS_Whitaker2014(z, lgMstar=10.5, use_equation_3 = False):
    # Whitaker,2014ApJ...795..104W
    # star forming main-sequence down to lower stellar mass at z<3
    # The evolution of sSFR, see their Eq.(6) on page 11, 
    #   sSFR/[Gyr^{-1}] = a * (1+z)^b
    #   where a and b are factors scaled with stellar mass. Using UV+IR SFR they obtain a b parameters as listed in Table 4. 
    # But at each redshift, alternatively see their Eq.(3) on page 7
    #   logSFR = a(z) (log(M*)-10.2) + b(z)
    # 
    # -- used by Genzel+2015 and Tacconi+2018
    # 
    spl_lgMstar_1 = np.array([9.2, 9.4, 9.6, 9.8, 10.0, 10.2, 10.4, 10.6, 10.8, 11.0 ])
    spl_lgMstar_2 = np.array([9.4, 9.6, 9.8, 10.0, 10.2, 10.4, 10.6, 10.8, 11.0, 11.2 ])
    spl_lga = np.array([-9.54, -9.50, -9.54, -9.58, -9.69, -9.93 , -10.11, -10.28, -10.53, -10.65 ])
    spl_b = np.array([1.95, 1.86, 1.90, 1.98, 2.16, 2.63, 2.88, 3.03, 3.37, 3.45 ])
    lga = np.interp(lgMstar, (spl_lgMstar_1+spl_lgMstar_2)/2.0, spl_lga, left=np.nan, right=np.nan)
    b = np.interp(lgMstar, (spl_lgMstar_1+spl_lgMstar_2)/2.0, spl_b, left=np.nan, right=np.nan)
    a = np.power(10.0, lga)
    sSFR_MS = (a * (1.+z)**b) # Gyr
    SFR_MS = sSFR_MS * np.power(10.0, lgMstar)
    # 
    if use_equation_3:
        spl_z_1 = np.array([0.5, 1.0, 1.5, 2.0 ])
        spl_z_2 = np.array([1.0, 1.5, 2.0, 2.5 ])
        spl_a_low = np.array([0.94, 0.99, 1.04, 0.91 ])
        spl_a_high = np.array([0.14, 0.51, 0.62, 0.67 ])
        spl_b = np.array([1.11, 1.31, 1.49, 1.62 ])
        a_low = np.interp(z, (spl_z_1+spl_z_2)/2.0, spl_a_low, left=np.nan, right=np.nan)
        a_high = np.interp(z, (spl_z_1+spl_z_2)/2.0, spl_a_high, left=np.nan, right=np.nan)
        b = np.interp(z, (spl_z_1+spl_z_2)/2.0, spl_b, left=np.nan, right=np.nan)
        SFR_MS_low = np.power(10.0, a_low * (lgMstar-10.2) + b)
        SFR_MS_high = np.power(10.0, a_high * (lgMstar-10.2) + b)
        SFR_MS = SFR_MS_low
        SFR_MS[(lgMstar>10.2)] = SFR_MS_high[(lgMstar>10.2)]
    # 
    return SFR_MS


def calc_SFR_MS_Lee2015(z, lgMstar=10.5):
    # Tomczak,2015ApJ...801...80L
    # star forming main-sequence with Herschel
    # see their Eq(2) on page 6
    # -- their equation can be rewritten as: SFR ~ 10^{s0} / ( 1 + (Mstar/M0)^{-gamma} )
    # -- their equation is also used in Tomczak2016
    spl_z_1 = np.array([0.25, 0.46, 0.63, 0.78, 0.93, 1.11])
    spl_z_2 = np.array([0.46, 0.63, 0.78, 0.93, 1.11, 1.30])
    spl_s0 = np.array([0.80, 0.99, 1.23, 1.35, 1.53, 1.72])
    spl_lgM0 = np.array([10.03, 9.82, 9.93, 9.96, 10.10, 10.31])
    spl_gamma = np.array([0.92, 1.13, 1.11, 1.28, 1.26, 1.07])
    lgM0 = np.interp(z, (spl_z_1+spl_z_2)/2.0, spl_lgM0, left=np.nan, right=np.nan)
    s0 = np.interp(z, (spl_z_1+spl_z_2)/2.0, spl_s0, left=np.nan, right=np.nan)
    gamma = np.interp(z, (spl_z_1+spl_z_2)/2.0, spl_gamma, left=np.nan, right=np.nan)
    log_SFR = s0 - np.log10(1.0 + np.power(np.power(10.0, lgMstar-lgM0), -gamma))
    SFR_MS = np.power(10.0, log_SFR)
    return SFR_MS


def calc_SFR_MS_Tomczak2016(z, lgMstar=10.5):
    # Tomczak,2016ApJ...817..118T
    # star forming main-sequence with ZFOURGE Survey
    # see their Eq(2) on page 9, and for star-forming galaxies Eq(4)
    # -- their equation can be rewritten as: SFR ~ 10^{s0} / ( 1 + (Mstar/M0)^{-gamma} )
    lgM0 = 9.458 + 0.865 * z - 0.132 * z**2
    s0 = 0.448 + 1.220 * z - 0.174 * z**2
    gamma = 1.091
    log_SFR = s0 - np.log10(1.0 + np.power(np.power(10.0, lgMstar-lgM0), -gamma))
    SFR_MS = np.power(10.0, log_SFR)
    return SFR_MS


def calc_SFR_MS_Pearson2018(z, lgMstar=10.5):
    # Pearson 2018 A&A
    # Herschel deblended IR SFR
    # abstract: log(SFR) = α[log(M⋆) − 10.5] + β
    #   β(z)=(1.10±0.07)+[(0.53±0.05)×log({0.03±0.11}+z)], Eq(9)
    #   α(z) = (0.38 ± 0.04) + (0.12 ± 0.02)z, Eq(10)
    Pearson2018_alpha = 0.38 + 0.12 * z
    Pearson2018_beta = 1.10 + (0.53 * np.log(0.03 + z)) # np.log is ln()
    SFR_MS = np.power(10.0, Pearson2018_alpha * (lgMstar-10.5) + Pearson2018_beta)
    sSFR_MS = SFR_MS / 10**(lgMstar) * 1e9 # Gyr
    return SFR_MS

def calc_SFR_MS_Leslie20180901(z, lgMstar, params = (0.5701,9.75469,4.1049,2.83077) ):
    # priv. comm.
    r = np.log10(1.+z)
    s0, m0, a1, a2 = params   
    S0_ = s0 + a1*r 
    M0_ = m0 + a2*r 
    gamma = 1.0
    log_SFR = S0_ - np.log10(1+10**(-gamma*(lgMstar-M0_)))
    SFR_MS = np.power(10.0, log_SFR)
    return SFR_MS

def calc_SFR_MS_Leslie20190111(z, lgMstar, params = (0.568,9.684,4.07,2.928,11.389) ):
    # priv. comm.
    # using a Chabrier IMF and a flat lambda CDM cosmology with H0=70,Om0=0.3.
    r = np.log10(1.+z)
    s0, m0, a1, a2, a3 = params   
    S0_ = s0 + a1*r 
    M0_ = m0 + a2*r 
    gamma = 1+np.exp(-a3*r)
    log_SFR = S0_ - np.log10(1+10**(-gamma*(lgMstar-M0_)))
    SFR_MS = np.power(10.0, log_SFR)
    return SFR_MS

def calc_SFR_MS_Leslie20190515(z, lgMstar, params = (0.38,9.28,4.49,3.70,7.71) ):
    # priv. comm.
    # using a Chabrier IMF and a flat lambda CDM cosmology with H0=70,Om0=0.3.
    r = np.log10(1.+z)
    s0, m0, a1, a2, a3 = params
    #S0_ = s0 + a1*r 
    #M0_ = m0 + a2*r 
    #gamma = 1+np.exp(-a3*r)
    log_SFR = s0 + a1*r - np.log10(1 + ( 10**lgMstar / 10**(m0+a2*r) )**(-(1.+np.exp(-a3*r))) )
    SFR_MS = np.power(10.0, log_SFR)
    return SFR_MS

def calc_SFR_MS_Leslie20190710(z, lgMstar, params = (0.38,9.28,4.48,3.69,7.70) ):
    # priv. comm.
    # using a Chabrier IMF and a flat lambda CDM cosmology with H0=70,Om0=0.3.
    r = np.log10(1.+z)
    s0, m0, a1, a2, a3 = params
    #S0_ = s0 + a1*r 
    #M0_ = m0 + a2*r 
    #gamma = 1+np.exp(-a3*r)
    log_SFR = s0 + a1*r - np.log10(1 + ( 10**lgMstar / 10**(m0+a2*r) )**(-(1.+np.exp(-a3*r))) )
    SFR_MS = np.power(10.0, log_SFR)
    return SFR_MS








