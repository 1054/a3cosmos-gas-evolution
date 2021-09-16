#!/usr/bin/env python
# 
# Usage:
#    from calc_alpha_CO import ( calc_metalZ_from_FMR_following_Mannucci2010, 
#                                calc_alphaCO_from_metalZ_following_Genzel2015a, 
#                                calc_alphaCO_from_metalZ_following_Genzel2015b )
# 

from __future__ import print_function

import os, sys, re, json, time, astropy
import numpy as np
from astropy.table import Table, Column, hstack
from copy import copy

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
def calc_alphaCO_from_metalZ_following_Wilson1995(metalZ):
    # Wilson, 1995ApJ...448L..97W
    # Eq. (4)
    metalZ_solar_old = 8.88
    metalZ_solar_new = 8.69
    XCO_Galactic = 3e20
    alphaCO_Galactic = 3e20/2e20*4.3
    return np.power(10.0, 5.95 - 0.67 * (metalZ-metalZ_solar_new+metalZ_solar_old)) * alphaCO_Galactic


def calc_alphaCO_from_metalZ_following_Magdis2012(metalZ):
    # Magdis, 2012ApJ...760....6M
    # Eq. (8)
    return np.power(10.0, (12.8) - (1.39) * metalZ )


def calc_alphaCO_from_metalZ_following_Genzel2012(metalZ):
    # Genzel 2012
    # Eq. (7) or Fig. 5 label
    return np.power(10.0, (12.0) - (1.3) * metalZ )
    #return np.power(10.0, (11.8) - (1.27) * metalZ )


def calc_alphaCO_from_metalZ_following_Genzel2015a(metalZ, metalZ_solar = 8.67, alphaCO_MilkyWay = 4.36):
    # Genzel 2015ApJ...800...20G
    # Eq. (6)
    # 
    # alphaCO_MilkyWay 4.36 ± 0.9 (M⊙ (K km s−1 pc2 )−1 ), 
    #   which is equivalent to XCO = N(H2)/(TRJ=1∆v) = 2 × 10^20 (cm−2 (K km s−1)−1; 
    #   Strong & Mattox 1996; Dame et al. 2001; Grenier et al. 2005; 
    #   Bolatto et al. 2008, 2013; Leroy et al. 2011; Abdo et al. 2010; 
    #   Ostriker et al. 2010).
    # 
    return 0.67 * np.exp(0.36 * np.power(10.0, -(metalZ-metalZ_solar)) ) * alphaCO_MilkyWay


def calc_alphaCO_from_metalZ_following_Genzel2015b(metalZ, metalZ_solar = 8.67, alphaCO_MilkyWay = 4.36):
    # Genzel 2015ApJ...800...20G
    # Eq. (7)
    # 
    # alphaCO_MilkyWay 4.36 ± 0.9 (M⊙ (K km s−1 pc2 )−1 ), 
    #   which is equivalent to XCO = N(H2)/(TRJ=1∆v) = 2 × 10^20 (cm−2 (K km s−1)−1; 
    #   Strong & Mattox 1996; Dame et al. 2001; Grenier et al. 2005; 
    #   Bolatto et al. 2008, 2013; Leroy et al. 2011; Abdo et al. 2010; 
    #   Ostriker et al. 2010).
    # 
    return np.power(10.0, -1.27 * (metalZ-metalZ_solar) ) * alphaCO_MilkyWay


def calc_alphaCO_from_metalZ_following_Genzel2015(metalZ, metalZ_solar = 8.67, alphaCO_MilkyWay = 4.36):
    return calc_alphaCO_from_metalZ_following_Genzel2015a(metalZ, metalZ_solar = metalZ_solar, alphaCO_MilkyWay = alphaCO_MilkyWay)


def calc_alphaCO_from_metalZ_following_Bolatto2013(metalZ, Sigma_total, Sigma_GMC_100 = 1.0, metalZ_solar = 8.67, alphaCO_MilkyWay = 4.36):
    # Bolatto, Wolfire, Leroy - 2013 - The CO-to-H2 Conversion Factor.pdf
    # Eq. (31) on PDF page 55.
    # 
    # Sigma_GMC_100 is the mean GMC gas surface density in unit of 100 Msolar pc^{-2}. 
    # Sigma_total is the total gas surface density in unit of 1 Msolar pc^{-2}. 
    # t_gamma is γ ≈ 0.5 for 􏰣total > 100 M⊙ pc−2 and γ = 0 otherwise
    if np.isscalar(Sigma_total):
        if Sigma_total >= 100.0:
            t_gamma = 0.5
        else:
            t_gamma = 0.0
    else:
        t_gamma = np.array(Sigma_total) * 0.0
        t_gamma[(Sigma_total>=100.0)] = 0.5
    return 2.9 * np.exp(0.4/(np.power(10.0,metalZ-metalZ_solar)*Sigma_GMC_100)) * np.power(Sigma_total/100.0,-t_gamma)


def calc_alphaCO_from_metalZ_following_Accurso2017(metalZ, DeltaMS = 0.0):
    # Accurso2017 calibrated alphaCO with L_CII/L_CO
    # see their Eq.(24)
    return np.power(10.0, 14.752 - 1.623 * metalZ + 0.062 * DeltaMS)


def calc_alphaCO_from_metalZ_following_Bertemes2018(metalZ):
    # Bertemes2018 used the geometric mean of Genzel2012 and Bolatto2013 alphaCO
    # see their Eq.(9),(10)
    # They assumed 
    alphaCO_G12 = np.power(10.0, -1.27 * (metalZ-8.67) )
    alphaCO_B13 = 0.67 * np.exp(0.36 * (-(metalZ-8.67)) )
    alphaCO_MW = 3.2
    alphaCO_combined = alphaCO_MW * np.sqrt(alphaCO_G12*alphaCO_B13) * 1.36
    return alphaCO_combined


def calc_alphaCO_from_metalZ_following_Tacconi2018(metalZ):
    # Tacconi2018 used the geometric mean of Genzel2012 and Bolatto2013 alphaCO (see also Equations (6) and (7) in Genzel2015)
    # see their Eq.(2)
    # They assumed 
    alphaCO_G12 = np.power(10.0, -1.27 * (metalZ-8.67) )
    alphaCO_B13 = 0.67 * np.exp(0.36 * (-(metalZ-8.67)) )
    alphaCO_MW = 4.36 # Helium included
    alphaCO_combined = alphaCO_MW * np.sqrt(alphaCO_G12*alphaCO_B13)
    return alphaCO_combined


def calc_alphaCO_from_metalZ_following_Boselli2014(lgL_Hband):
    # Boselli2014a (2014A&A...564A..65B) used a luminosity dependent alphaCO
    # see their page 2, paragraph 2.
    # 
    X_CO = np.power(10.0, -0.38 * lgL_Hband + 24.23 )
    alphaCO = X_CO / 2.0e20 * 4.3
    return alphaCO







