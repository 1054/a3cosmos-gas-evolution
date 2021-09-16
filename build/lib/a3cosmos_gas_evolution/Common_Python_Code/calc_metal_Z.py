#!/usr/bin/env python
# 
# Usage:
#    from calc_metal_Z import ( calc_metalZ_from_FMR_following_Mannucci2010, 
#                               calc_metalZ_from_FMR_following_Genzel2015a, 
#                               calc_metalZ_from_FMR_following_Genzel2015b, 
#                               calc_metalZ_from_FMR_following_Genzel2015ab_combined_by_dzliu, 
#                               convert_metalZ_M08_to_metalZ_PP04, 
#                               convert_metalZ_D02_to_metalZ_PP04, 
#                               convert_metalZ_KK04_to_metalZ_PP04, 
#                             )
# 
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
def calc_metalZ_from_FMR_following_Kewley2008_PP04_O3N2(M_star):
    # metalZ is 12+log(O/H)_{PP04(O3N2)}
    # 
    # Kewley & Ellison (2008ApJ...681.1183K)
    # Table 2
    # Note that their caption has typo: a + b*x + c*x**2 + d*x**3
    x = np.log10(M_star)
    a = 32.1488
    b = -8.51258
    c = 0.976384
    d = -0.0359763
    metalZ_PP04 = a + b*x + c*x**2 + d*x**3
    # 
    return metalZ_PP04


def calc_metalZ_from_FMR_following_Kewley2008_PP04_N2(M_star):
    # metalZ is 12+log(O/H)_{PP04(N2)}
    # 
    # Kewley & Ellison (2008ApJ...681.1183K)
    # Table 2
    # Note that their caption has typo: a + b*x + c*x**2 + d*x**3
    x = np.log10(M_star)
    a = 23.9049
    b = -5.62784
    c = 0.645142
    d = -0.0235065
    metalZ_PP04 = a + b*x + c*x**2 + d*x**3
    # 
    return metalZ_PP04


def calc_metalZ_from_FMR_following_Kewley2008_KK04(M_star):
    # metalZ is 12+log(O/H)_{KK04(R23)}
    # 
    # Kewley & Ellison (2008ApJ...681.1183K)
    # Table 2
    # Note that their caption has typo: a + b*x + c*x**2 + d*x**3
    x = np.log10(M_star)
    a = 27.7911
    b = -6.94493
    c = 0.808097
    d = -0.0301508
    metalZ_KK04 = a + b*x + c*x**2 + d*x**3
    # 
    return metalZ_KK04


def calc_metalZ_from_FMR_following_Kewley2008_KD02(M_star):
    # metalZ is 12+log(O/H)_{KD02}
    # 
    # Kewley & Ellison (2008ApJ...681.1183K)
    # Table 2
    # Note that their caption has typo: a + b*x + c*x**2 + d*x**3
    x = np.log10(M_star)
    a = 28.0974
    b = -7.23631
    c = 0.850344
    d = -0.0318315
    metalZ_KD02 = a + b*x + c*x**2 + d*x**3
    # 
    return metalZ_KD02


def calc_metalZ_from_FMR_following_Mannucci2010(M_star, SFR):
    # metalZ is 12+log(O/H)_{Maiolino2008}
    # 
    # Mannucci2010
    # Sect. 4, Eq. 2
    m = np.log10(M_star) - 10
    s = np.log10(SFR)
    metalZ_M08 = 8.90 + 0.37*m - 0.14*s - 0.19*m**2 + 0.12*m*s - 0.054*s**2
    # 
    return metalZ_M08


def calc_metalZ_from_FMR_following_Mannucci2010_Eq2(M_star, SFR):
    # metalZ is 12+log(O/H)_{Maiolino2008}
    # 
    # Mannucci2010
    # Sect. 4, Eq. 2, two variables
    m = np.log10(M_star) - 10
    s = np.log10(SFR)
    metalZ_M08 = 8.90 + 0.37*m - 0.14*s - 0.19*m**2 + 0.12*m*s - 0.054*s**2
    # 
    return metalZ_M08


def calc_metalZ_from_FMR_following_Mannucci2010_Eq4(M_star, SFR):
    # metalZ is 12+log(O/H)_{Maiolino2008}
    # 
    # Mannucci2010
    # Sect. 4, Eq. 4, a single variable
    x = np.log10(M_star) - 0.32 * np.log10(SFR) - 10
    metalZ_M08 = 8.90 + 0.39 * x - 0.20 * x**2 - 0.077 * x**3 + 0.064 * x**4
    # 
    return metalZ_M08


def calc_metalZ_from_FMR_following_Mannucci2011(M_star, SFR):
    # metalZ is 12+log(O/H)_{Maiolino2008}
    # 
    # Mannucci 2011 - 2011MNRAS.414.1263M.pdf
    # Sect. 2, Eq. 2
    m = np.log10(M_star) - 10
    s = np.log10(SFR)
    mu032 = np.log10(M_star) - 0.32 * np.log10(SFR)
    x = np.log10(M_star) - 0.32 * np.log10(SFR) - 10
    if np.isscalar(mu032):
        if mu032 >= 9.5:
            metalZ_M08 = 8.90 + 0.37*m - 0.14*s - 0.19*m**2 + 0.12*m*s - 0.054*s**2
            #metalZ_M08 = 8.90 + 0.39 * x - 0.20 * x**2 - 0.077 * x**3 + 0.064 * x**4
        else:
            metalZ_M08 = 8.93 + 0.51*(mu032-10)
    else:
        if np.isscalar(s) and ~np.isscalar(m):
            s = m*0.0 + s
        if np.isscalar(m) and ~np.isscalar(s):
            m = s*0.0 + m
        mask = (mu032 >= 9.5)
        metalZ_M08 = 8.90 + 0.37*m - 0.14*s - 0.19*m**2 + 0.12*m*s - 0.054*s**2
        #metalZ_M08 = 8.90 + 0.39 * x - 0.20 * x**2 - 0.077 * x**3 + 0.064 * x**4
        metalZ_M08[~mask] = 8.93 + 0.51*(mu032[~mask]-10)
    # 
    return metalZ_M08


def calc_metalZ_from_FMR_following_Maiolino2008(M_star, z):
    # metalZ is 12+log(O/H)_{Maiolino2008}
    # 
    # Maiolino2008 doi:10.1051/0004-6361/201220074
    # Sect. 7.6, Eq. 2
    # 12+log(O/H) = -0.0864 * (logMstar - logM0)**2 + K0
    # Table 5
    #   z     logM0  K0
    #   0.07  11.18  9.04
    #   0.7   11.57  9.04
    #   2.2   12.38  8.99
    #   3.5   12.76  8.79
    #   3.5   12.87  8.90
    # 
    m = np.log10(M_star)
    spl_z = np.array([0.07,0.7,2.2,3.5])
    spl_m0 = np.array([11.18,11.57,12.38,12.76])
    spl_K0 = np.array([9.04,9.04,8.99,8.79])
    m0 = np.interp(z, spl_z, spl_m0)
    K0 = np.interp(z, spl_z, spl_K0)
    metalZ_M08 = -0.0864 * (m - m0)**2 + K0
    print('calc_metalZ_from_FMR_following_Maiolino2008: logM0 = %s, K0 = %s'%(m0, K0))
    # 
    return metalZ_M08


def calc_metalZ_from_FMR_following_Magnelli2012(M_star, z):
    # metalZ is 12+log(O/H)_{Denicolo2002}
    # 
    # Magnelli2012 doi:10.1051/0004-6361/201220074
    # Sect. 4, Eq. 5
    # The Denicoló et al. calibration system provides the best agreement, at high stellar masses, between all different metallicity calibrators (Kewley & Ellison 2008)
    # We note that metallicities estimated using the “fundamental metallicity relation” of Mannucci et al. (2010) extended to z ∼ 2 in Genzel et al. (2012) are consistent with those derived using Eq. (5).
    m = np.log10(M_star)
    if z < 1.5:
        a = -4.45
    else:
        a = -4.51
    metalZ_D02 = (2.18 * m - 0.0896 * m**2 + a)
    # 
    return metalZ_D02


def calc_metalZ_from_FMR_following_Genzel2015_Eq12a_with_dzliu_limit(M_star, z):
    # Genzel et al. 2015:
    #   Following Maiolino et al. (2008) we combined the mass–metallicity relations 
    #   at different redshifts presented by Erb et al. (2006), 
    #   Maiolino et al. (2008), Zahid et al. (2014), and Wuyts et al. (2014) 
    #   in the following fitting function: 
    #     12+log(O/H) = a - 0.087 * (lgMstar - b)**2, 
    #     a = 8.74, and 
    #     b = 10.4(0.05) + 4.46(0.3) * lg(1+z) - 1.78(0.4) * (lg(1+z))**2  -- Eq. (12a)
    # 
    # dzliu_limit:
    #     when np.log10(M_star) > b, set (np.log10(M_star) - b) to zero. 
    # 
    a = 8.74
    b = 10.4 + 4.46 * np.log10(1+z) - 1.78 * (np.log10(1+z))**2
    #metalZ_PP04 = a - 0.087 * (np.log10(M_star) - b)**2
    diff = np.log10(M_star) - b
    if np.isscalar(diff):
        if diff<0.0:
            diff = 0.0 
    else:
        diff[(diff<0.0)] = 0.0
    # 
    metalZ_PP04 = a - 0.087 * (diff)**2
    # 
    return metalZ_PP04


def calc_metalZ_from_FMR_following_Genzel2015_Eq12a(M_star, z):
    # Genzel et al. 2015:
    #   Following Maiolino et al. (2008) we combined the mass–metallicity relations 
    #   at different redshifts presented by Erb et al. (2006), 
    #   Maiolino et al. (2008), Zahid et al. (2014), and Wuyts et al. (2014) 
    #   in the following fitting function: 
    #     12+log(O/H) = a - 0.087 * (lgMstar - b)**2, 
    #     a = 8.74, and 
    #     b = 10.4(0.05) + 4.46(0.3) * lg(1+z) - 1.78(0.4) * (lg(1+z))**2  -- Eq. (12a)
    a = 8.74
    b = 10.4 + 4.46 * np.log10(1+z) - 1.78 * (np.log10(1+z))**2
    metalZ_PP04 = a - 0.087 * (np.log10(M_star) - b)**2
    # 
    return metalZ_PP04


def calc_metalZ_from_FMR_following_Genzel2015a(M_star, SFR, z):
    # Genzel et al. 2015:
    #   Following Maiolino et al. (2008) we combined the mass–metallicity relations 
    #   at different redshifts presented by Erb et al. (2006), 
    #   Maiolino et al. (2008), Zahid et al. (2014), and Wuyts et al. (2014) 
    #   in the following fitting function: 
    #     12+log(O/H) = a - 0.087 * (lgMstar - b)**2, 
    #     a = 8.74, and 
    #     b = 10.4(0.05) + 4.46(0.3) * lg(1+z) - 1.78(0.4) * (lg(1+z))**2  -- Eq. (12a)
    a = 8.74
    b = 10.4 + 4.46 * np.log10(1+z) - 1.78 * (np.log10(1+z))**2
    metalZ_PP04 = a - 0.087 * (np.log10(M_star) - b)**2
    # 
    return metalZ_PP04


def calc_metalZ_from_FMR_following_Genzel2015_Eq12b(M_star, SFR):
    # Genzel et al. 2015 method 2:
    #   Mannucci et al. (2010) presented evidence for a dependence of metallicity 
    #   on the SFR for z ∼ 0 SDSS galaxies at a given stellar mass (the fundamental 
    #   metallicity relation), yielding an alternative version of Equation (12a):
    #     dZM08 = 12+log(O/H)_{M08} - 8.69 = 0.21 + 0.39 * x - 0.2 * x**2 - 0.077 * x**3 + 0.064 * x**4
    #     with x = logMstar - 0.32 * logSFR - 10, and
    #     12+log(O/H)_{PP04} - 8.9 = -0.4408 + 0.7044 * dZM08 - 0.1602 * dZM08**2 - 0.4105 * dZM08**3 - 0.1898 * dZM08**4
    # 
    #x = np.log10(M_star) - 0.32 * np.log10(SFR) - 10
    #dZM08 = 0.21 + 0.39 * x - 0.2 * x**2 - 0.077 * x**3 + 0.064 * x**4
    # 
    dZM08 = calc_metalZ_from_FMR_following_Mannucci2010_Eq4(M_star, SFR) - 8.69
    # 
    #dZM08 = calc_metalZ_from_FMR_following_Mannucci2010_Eq2(M_star, SFR) - 8.69
    # 
    metalZ = 8.9 + (-0.4408 + 0.7044 * dZM08 - 0.1602 * dZM08**2 - 0.4105 * dZM08**3 - 0.1898 * dZM08**4)
    # 
    return metalZ


def calc_metalZ_from_FMR_following_Genzel2015b(M_star, SFR):
    # Genzel et al. 2015 method 2:
    #   Mannucci et al. (2010) presented evidence for a dependence of metallicity 
    #   on the SFR for z ∼ 0 SDSS galaxies at a given stellar mass (the fundamental 
    #   metallicity relation), yielding an alternative version of Equation (12a):
    #     dZM08 = 12+log(O/H)_{M08} - 8.69 = 0.21 + 0.39 * x - 0.2 * x**2 - 0.077 * x**3 + 0.064 * x**4
    #     with x = logMstar - 0.32 * logSFR - 10, and
    #     12+log(O/H)_{PP04} - 8.9 = -0.4408 + 0.7044 * dZM08 - 0.1602 * dZM08**2 - 0.4105 * dZM08**3 - 0.1898 * dZM08**4
    # 
    #x = np.log10(M_star) - 0.32 * np.log10(SFR) - 10
    #dZM08 = 0.21 + 0.39 * x - 0.2 * x**2 - 0.077 * x**3 + 0.064 * x**4
    # 
    dZM08 = calc_metalZ_from_FMR_following_Mannucci2010_Eq4(M_star, SFR) - 8.69
    # 
    #dZM08 = calc_metalZ_from_FMR_following_Mannucci2010_Eq2(M_star, SFR) - 8.69
    # 
    metalZ = 8.9 + (-0.4408 + 0.7044 * dZM08 - 0.1602 * dZM08**2 - 0.4105 * dZM08**3 - 0.1898 * dZM08**4)
    # 
    return metalZ


def calc_metalZ_from_FMR_following_Genzel2015ab_combined_by_dzliu(M_star, SFR, z):
    # dzliu combined: 
    #   Genzel et al. 2015 method 1 and Genzel et al. 2015 method 2
    # when lgMstar >= 10.0, we use calc_metalZ_from_FMR_following_Genzel2015b()
    # otherwise we use calc_metalZ_from_FMR_following_Genzel2015a() but scaled to make sure the prediction is smooth at lgMstar=10.0
    # 
    if np.isscalar(M_star):
        input_M_star = np.array([M_star])
    else:
        input_M_star = np.array(M_star)
    # 
    if np.isscalar(SFR):
        input_SFR = np.array([SFR])
    else:
        input_SFR = np.array(SFR)
    # 
    if np.isscalar(z):
        input_z = np.array([z])
    else:
        input_z = np.array(z)
    # 
    if len(input_SFR) == 1 and len(input_M_star) > 1:
        input_SFR = np.array([input_SFR[0]]*len(input_M_star))
    # 
    if len(input_z) == 1 and len(input_M_star) > 1:
        input_z = np.array([input_z[0]]*len(input_M_star))
    # 
    mask = (input_M_star>=1e10)
    output_metalZ = input_M_star*0.0
    output_metalZ[mask] = calc_metalZ_from_FMR_following_Genzel2015b(input_M_star[mask], input_SFR[mask])
    output_metalZ[~mask] = calc_metalZ_from_FMR_following_Genzel2015a(input_M_star[~mask], input_SFR[~mask], input_z[~mask])
    renorm_metalZ = calc_metalZ_from_FMR_following_Genzel2015b(1e10, input_SFR[~mask]) - \
                    calc_metalZ_from_FMR_following_Genzel2015a(1e10, input_SFR[~mask], input_z[~mask])
                    
    output_metalZ[~mask] = output_metalZ[~mask] + renorm_metalZ
    # 
    return output_metalZ


def calc_metalZ_from_FMR_following_Guo2016(M_star, SFR):
    # metalZ is 12+log(O/H)_{Maiolino2008}
    # 
    # Yicheng Guo et al. 2016 (doi:10.3847/0004-637X/822/2/103)
    # Abstract: 12 + log(O/H) = (5.83 ± 0.19)+(0.30 ± 0.02)*log(Mstar)
    # The [O III]/Hβ flux ratio is then converted to metallicity through the calibration of Maiolino et al. (2008, M08).
    m = np.log10(M_star)
    metalZ_M08 = 5.83 + 0.30 * m
    # 
    return metalZ_M08









def convert_metalZ_M08_to_metalZ_PP04(metalZ_M08):
    # PP04: Pettini 2004
    #       N2 \equiv lg(NII6583/Halpha)
    #       metalZ_PP04 \equiv 12+lg(O/H) = 9.37 + 2.03 * N2 + 1.26 * N2**2 + 0.32 * N2**3      (PP04 Eq.1)
    #       or metalZ_PP04 = 8.90 + 0.57 * N2                                                   (PP04 Eq.2, valid in the range -2.5 < N2 < -0.3)
    #       they mentoned that the first polynomial form is not significantly better than the second linear form. 
    # M08:  Maiolino 2008
    #       N2 = -0.7732 + 1.2357 * x - 0.2811 * x**2 - 0.7201 * x**3 - 0.3330 * x**4
    #       where x = 12+lg(O/H) - 8.69 = metalZ_M08 - 8.69
    # so, the conversion from metalZ_M08 to metalZ_PP04 is:
    #       metalZ_PP04 = 8.90 + 0.57 * N2
    #                   = 8.90 + 0.57 * (-0.7732 + 1.2357 * x - 0.2811 * x**2 - 0.7201 * x**3 - 0.3330 * x**4)
    # same as Genzel+2015ApJ...800...20G Eq.12b.
    # 
    x = metalZ_M08 - 8.69
    return 8.90 + 0.57 * (-0.7732 + 1.2357 * x - 0.2811 * x**2 - 0.7201 * x**3 - 0.3330 * x**4)


def convert_metalZ_M08_to_metalZ_PP04_N2_polynomial(metalZ_M08):
    # PP04: Pettini 2004
    #       N2 \equiv lg(NII6583/Halpha)
    #       metalZ_PP04 \equiv 12+lg(O/H) = 9.37 + 2.03 * N2 + 1.26 * N2**2 + 0.32 * N2**3      (PP04 Eq.1)
    #       or metalZ_PP04 = 8.90 + 0.57 * N2                                                   (PP04 Eq.2, valid in the range -2.5 < N2 < -0.3)
    #       they mentoned that the first polynomial form is not significantly better than the second linear form. 
    # M08:  Maiolino 2008
    #       N2 = -0.7732 + 1.2357 * x - 0.2811 * x**2 - 0.7201 * x**3 - 0.3330 * x**4
    #       where x = 12+lg(O/H) - 8.69 = metalZ_M08 - 8.69
    # so, the conversion from metalZ_M08 to metalZ_PP04 is:
    #       metalZ_PP04 = 9.37 + 2.03 * N2 + 1.26 * N2**2 + 0.32 * N2**3
    #                N2 = -0.7732 + 1.2357 * x - 0.2811 * x**2 - 0.7201 * x**3 - 0.3330 * x**4
    # same as Genzel+2015ApJ...800...20G Eq.12b.
    # 
    x = metalZ_M08 - 8.69
    N2 = -0.7732 + 1.2357 * x - 0.2811 * x**2 - 0.7201 * x**3 - 0.3330 * x**4
    metalZ_PP04 = 9.37 + 2.03 * N2 + 1.26 * N2**2 + 0.32 * N2**3
    return metalZ_PP04


def convert_metalZ_D02_to_metalZ_PP04(metalZ_D02):
    # D02:  Denicolo et al. (2002) Eq.(2)
    #       N2 \equiv lg(NII6584/Halpha)
    #       metalZ_D02 \equiv 12+lg(O/H) = 9.12 + 0.73 * N2
    #       they mentoned that the first polynomial form is not significantly better than the second linear form. 
    # so, the conversion from metalZ_D02 to metalZ_PP04 (N2)
    #       according to Kewley & Ellison (2008ApJ...681.1183K) Table 3 is:
    # 
    a = -444.7831
    b = 165.42600
    c = -20.202000
    d = 0.8249386
    x = metalZ_D02
    return (a + b * x + c * x**2 + d * x**3)


def convert_metalZ_KK04_to_metalZ_PP04(metalZ_KK04):
    # KK04: 
    #       N2 \equiv lg(NII6584/Halpha)
    #       
    # so, the conversion from metalZ_KK04 to metalZ_PP04 (N2) 
    #       according to Kewley & Ellison (2008ApJ...681.1183K) Table 3 is:
    # 
    a = 916.7484
    b = -309.54480
    c = 35.051680
    d = -1.3188000
    x = metalZ_KK04
    return (a + b * x + c * x**2 + d * x**3)


def convert_metalZ_KK04_to_metalZ_PP04_O3N2(metalZ_KK04):
    # KK04: 
    #       N2 \equiv lg(NII6584/Halpha)
    #       
    # so, the conversion from metalZ_KK04 to metalZ_PP04 (O3N2) 
    #       according to Kewley & Ellison (2008ApJ...681.1183K) Table 3 is:
    # 
    a = 631.2766
    b = -210.02090
    c = 23.483050
    d = -0.8704286
    x = metalZ_KK04
    return (a + b * x + c * x**2 + d * x**3)


def convert_metalZ_KD02_to_metalZ_PP04(metalZ_KD02):
    # KD02: 
    #       Kewley & Dopita (2002)
    #       
    # so, the conversion from metalZ_KD02 to metalZ_PP04 (N2) 
    #       according to Kewley & Ellison (2008ApJ...681.1183K) Table 3 is:
    # 
    a = 569.4927
    b = -192.51820
    c = 21.918360
    d = -0.8278840
    x = metalZ_KD02
    return (a + b * x + c * x**2 + d * x**3)


def convert_metalZ_KD02_to_metalZ_PP04_O3N2(metalZ_KD02):
    # KD02: 
    #       Kewley & Dopita (2002)
    #       
    # so, the conversion from metalZ_KD02 to metalZ_PP04 (O3N2) 
    #       according to Kewley & Ellison (2008ApJ...681.1183K) Table 3 is:
    # 
    a = 664.8453
    b = -225.75330
    c = 25.768880
    d = -0.9761368
    x = metalZ_KD02
    return (a + b * x + c * x**2 + d * x**3)


def convert_metalZ_PP04_N2_to_metalZ_PP04_O3N2(metalZ_PP04_N2):
    # PP04 (N2): 
    #       Pettini 2004
    #       N2 \equiv lg(NII6583/Halpha)
    #       metalZ_PP04_N2_linear = 8.90 + 0.57 * N2                                        (PP04 Eq.2, valid in the range -2.5 < N2 < -0.3)
    #       metalZ_PP04_N2_polynomial = 9.37 + 2.03 * N2 + 1.26 * N2**2 + 0.32 * N2**3      (PP04 Eq.1)
    #       
    # so, the conversion from metalZ_PP04_N2 to metalZ_PP04 (O3N2) 
    #       according to Kewley & Ellison (2008ApJ...681.1183K) Table 3 is:
    # 
    a = -8.0069
    b = 2.74353
    c = -0.093680
    d = 0.0
    x = metalZ_PP04_N2
    return (a + b * x + c * x**2 + d * x**3)


def convert_metalZ_PP04_N2_linear_to_metalZ_PP04_N2_polynomial(metalZ_PP04_N2_linear):
    # PP04 (N2): 
    #       Pettini 2004
    #       N2 \equiv lg(NII6583/Halpha)
    #       metalZ_PP04_N2_linear = 8.90 + 0.57 * N2                                        (PP04 Eq.2, valid in the range -2.5 < N2 < -0.3)
    #       metalZ_PP04_N2_polynomial = 9.37 + 2.03 * N2 + 1.26 * N2**2 + 0.32 * N2**3      (PP04 Eq.1)
    #       
    # so, the conversion from metalZ_PP04_N2_linear to metalZ_PP04_N2_polynomial:
    # 
    N2 = (metalZ_PP04_N2_linear - 8.90) / 0.57
    metalZ_PP04_N2_polynomial = 9.37 + 2.03 * N2 + 1.26 * N2**2 + 0.32 * N2**3
    return metalZ_PP04_N2_polynomial



















# 
# dzliu selection of the best function
# 
def calc_metalZ_from_FMR_with_dzliu_selection(M_star, SFR, z):
    # metalZ is 12+log(O/H)_{PP04}
    # use Genzel+2015 Eq.12a (variable: z) when it is lower than Mannucci+2010 Eq.4 (SFR)
    # but at the massive end, if Genzel+2015 Eq.12a decreases with increasing stellar mass, we set it to not decrease. 
    input_M_star = np.log10(M_star)
    if np.isscalar(M_star):
        input_M_star = np.array([M_star])
    else:
        input_M_star = np.array(M_star)
    if np.isscalar(SFR):
        input_SFR = np.array([SFR])
    else:
        input_SFR = np.array(SFR)
    if np.isscalar(z):
        input_z = np.array([z])
    else:
        input_z = np.array(z)
    # 
    if len(input_z) == 1 and len(input_M_star) > 1:
        input_z = np.array([input_z[0]]*len(input_M_star))
        #print('calc_metalZ_from_FMR_with_dzliu_selection() Replicating input_z')
    if len(input_SFR) == 1 and len(input_M_star) > 1:
        input_SFR = np.array([input_SFR[0]]*len(input_M_star))
        #print('calc_metalZ_from_FMR_with_dzliu_selection() Replicating input_SFR')
    # 
    metalZ_G15a = calc_metalZ_from_FMR_following_Genzel2015_Eq12a(input_M_star, input_z)
    #metalZ_G15b = calc_metalZ_from_FMR_following_Genzel2015_Eq12b(input_M_star, input_SFR) # equals convert_metalZ_M08_to_metalZ_PP04_N2_linear(calc_metalZ_from_FMR_following_Mannucci2010_Eq4(input_M_star, input_SFR))
    metalZ_M10Eq4 = convert_metalZ_M08_to_metalZ_PP04_N2_polynomial(calc_metalZ_from_FMR_following_Mannucci2010_Eq4(input_M_star, input_SFR)) # apply convert_metalZ_M08_to_metalZ_PP04_N2_polynomial() conversion is better because it will turn down a little bit the metalZ at the low mass end!
    # 
    # choose min of (metalZ_G15a, metalZ_M10Eq4)
    mask1 = (metalZ_G15a < metalZ_M10Eq4)
    metalZ_out = metalZ_M10Eq4
    metalZ_out[mask1] = metalZ_G15a[mask1]
    #for i in range(1,len(metalZ_out)):
    #    if metalZ_out[i] < metalZ_out[i-1]:
    #        metalZ_out[i] = metalZ_out[i-1]
    # 
    # ref Mstar
    ref_M_star = 10**(10.4 + 4.46 * np.log10(1+input_z) - 1.78 * (np.log10(1+input_z))**2) # characteristic stellar mass deduced from Genzel2015_Eq12a
    ref_metalZ_G15a = calc_metalZ_from_FMR_following_Genzel2015_Eq12a(ref_M_star, input_z)
    #ref_metalZ_G15b = calc_metalZ_from_FMR_following_Genzel2015_Eq12b(ref_M_star, input_SFR) # equals convert_metalZ_M08_to_metalZ_PP04_N2_linear(calc_metalZ_from_FMR_following_Mannucci2010_Eq4(ref_M_star, input_SFR))
    ref_metalZ_M10Eq4 = convert_metalZ_M08_to_metalZ_PP04_N2_polynomial(calc_metalZ_from_FMR_following_Mannucci2010_Eq4(ref_M_star, input_SFR)) # apply convert_metalZ_M08_to_metalZ_PP04_N2_polynomial() conversion is better because it will turn down a little bit the metalZ at the low mass end!
    # 
    # choose min of (metalZ_G15a, metalZ_M10Eq4) when (M_star < ref_M_star)
    mask1 = np.logical_and(metalZ_G15a < metalZ_M10Eq4, input_M_star < ref_M_star)
    metalZ_out = metalZ_M10Eq4
    metalZ_out[mask1] = metalZ_G15a[mask1]
    # 
    # keeps constant at massive end
    #mask2 = (ref_metalZ_G15a < ref_metalZ_M10Eq4)
    #ref_metalZ_out = ref_metalZ_M10Eq4
    #ref_metalZ_out[mask2] = ref_metalZ_G15a[mask2]
    #mask3 = (input_M_star >= ref_M_star)
    #metalZ_out[mask3] = ref_metalZ_out[mask3]
    # 
    # keeps constant at massive end but use G15a
    mask3 = (input_M_star >= ref_M_star)
    metalZ_out[mask3] = ref_metalZ_G15a[mask3]
    # 
    return metalZ_out

















