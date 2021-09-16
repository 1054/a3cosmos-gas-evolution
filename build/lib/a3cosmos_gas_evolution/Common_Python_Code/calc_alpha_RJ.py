#!/usr/bin/env python
# 
# Usage:
#    from calc_alpha_RJ import ( calc_alphaRJmol_following_dzliu, 
#                                calc_alphaRJtot_following_dzliu, 
#                                calc_alphaRJmol_Hughes2017, 
#                                calc_alphaRJmol_Hughes2017 )
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







