#!/usr/bin/env python
# 
# Usage:
#    from calc_delta_GD import ( calc_deltaGD_from_metalZ_following_Leroy2011, 
#                                calc_deltaGD_from_metalZ_following_Magdis2012, 
#                                calc_deltaGD_from_metalZ_following_RemyRuyer2014a, 
#                                calc_deltaGD_from_metalZ_following_RemyRuyer2014b, 
#                                calc_deltaGD_from_metalZ_following_Genzel2015, 
#                              )
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
def calc_deltaGD_from_metalZ_following_Leroy2011(metalZ):
    # Leroy, 2011ApJ...737...12L
    # Sect. 5.1 Paragraph 2
    # metalZ is 12+log(O/H)_{PP04}
    # deltaGD = Mgas/Mdust
    # here Mgas inclues HI + H2 + Helium, i.e., atomic + mol. gas
    # They assume that δGDR does not vary between the atomic and molecular ISM. 
    return np.power(10.0, (9.4) - (0.85) * metalZ )


def calc_deltaGD_from_metalZ_following_Magdis2012(metalZ):
    # Magdis, 2012ApJ...760....6M
    # Sect. 4.2 Paragraph 1
    # 
    # metalZ is 12+log(O/H)_{PP04}
    # deltaGD = Mgas/Mdust, here Mgas includes HI + H2 + Helium
    return np.power(10.0, (10.54) - (0.99) * metalZ )


def calc_deltaGD_from_metalZ_following_RemyRuyer2014a(metalZ):
    # Remy-Ruyer et al. 2014 (doi:10.1051/0004-6361/201322803), broken power law, high metallicity part. See their Table 1. 
    #   catalog data: http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=J/A%2bA/582/A121/sample
    # XCO,MW case
    # includes HI + H2 + He + gaseous metal mass (= Z[O/H] * Mgas)
    #   Mgas = MHI + MH2 + MHe = 1.38 * (MHI + MH2)
    #   Z[O/H]_solar = 0.014 (Asplund et al. 2009)
    #   metalZ = 12+log(O/H)
    #   metalZ_solar = 8.69 (Asplund et al. 2009)
    # one example
    #   SBS 1533+574, metalZ = 8.05, MH2_alphaCO_MW = 3.70e7, MH2_alphaCO_Z = 2.03e9, Mdust.= 10**5.89, G/D = 47 - 2615. 
    #   NGC 1569, metalZ = 8.02, MH2_alphaCO_MW = 7.08e5, MH2_alphaCO_Z = 1.54e7, Mdust = 10**5.56, G/D = 2 - 42. 
    a = 2.21
    aH = 1.00
    b = 0.68
    aL = 3.08
    metaZ = metalZ
    metaZknee = 7.96
    metaZsolar = 8.69
    if not np.isscalar(metaZ):
        maskZ = (metaZ>metaZknee)
        GDR_ISM = metaZ * 0.0
        GDR_ISM[maskZ] = 10**(a + aH*(metaZsolar-metaZ[maskZ]))
        GDR_ISM[~maskZ] = 10**(b + aL*(metaZsolar-metaZ[~maskZ]))
    else:
        GDR_ISM = 10**(a + aH*(metaZsolar-metaZ)) if metaZ>7.96 else 10**(b + aL*(metaZsolar-metaZ))
    #GDR_mol = GDR_ISM/1.38 # considering only molecular gas (without HI, but still with Helium). See Remy-Ruyer et al. 2014 Eq.(3)
    return GDR_ISM


def calc_deltaGD_from_metalZ_following_RemyRuyer2014b(metalZ):
    # Remy-Ruyer et al. 2014 (doi:10.1051/0004-6361/201322803), broken power law, high metallicity part. See their Table 1. 
    # XCO,Z case
    # includes HI + H2 + He + gaseous metal mass (= Z[O/H] * Mgas)
    #   Mgas = MHI + MH2 + MHe = 1.38 * (MHI + MH2)
    #   Z[O/H]_solar = 0.014 (Asplund et al. 2009)
    #   metalZ = 12+log(O/H)
    #   metalZ_solar = 8.69 (Asplund et al. 2009)
    a = 2.21
    aH = 1.00
    b = 0.96
    aL = 3.10
    metaZ = metalZ
    metaZknee = 8.10
    metaZsolar = 8.69
    if not np.isscalar(metaZ):
        maskZ = (metaZ>metaZknee)
        GDR_ISM = metaZ * 0.0
        GDR_ISM[maskZ] = 10**(a + aH*(metaZsolar-metaZ[maskZ]))
        GDR_ISM[~maskZ] = 10**(b + aL*(metaZsolar-metaZ[~maskZ]))
    else:
        GDR_ISM = 10**(a + aH*(metaZsolar-metaZ)) if metaZ>7.96 else 10**(b + aL*(metaZsolar-metaZ))
    #GDR_mol = GDR_ISM/1.38 # considering only molecular gas (without HI, but still with Helium). See Remy-Ruyer et al. 2014 Eq.(3)
    return GDR_ISM


def calc_deltaGD_from_metalZ_following_Genzel2015(metalZ):
    # Genzel et al. 2015 (doi:10.1088/0004-637X/800/1/20), Eq. 10
    #   Following Magdis et al. (2012a) and Magnelli et al. (2012a) we converted the Draine & Li (2007) 
    #   model dust masses to (molecular) gas masses by applying the metallicity dependent dust-to-gas 
    #   ratio fitting function for z ∼ 0 SFGs found by Leroy et al. (2011):
    #return 1.0/10**(-2. + 0.85*(metalZ-8.67) ) # almost exactly the same as Leroy2011, note that they express this as dust/gas instead of gas/dust.
    return calc_deltaGD_from_metalZ_following_Leroy2011(metalZ)








