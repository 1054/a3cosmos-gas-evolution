#!/usr/bin/env python
# 
# usage:
#     from calc_galaxy_size_evolution import (calc_galaxy_size_Tacconi2018,
#                                             )
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
def calc_galaxy_size_Tacconi2018(z, Mstar):
    # Tacconi2018 Sect. 4.1
    #Mstar = 10**np.array(lgMstar)
    Re0 = 8.9 * (1.+z)**(-0.75) * (Mstar/5e10)**(0.23) # kpc
    return Re0







