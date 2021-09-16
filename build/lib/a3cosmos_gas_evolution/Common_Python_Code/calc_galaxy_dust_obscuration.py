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
def calc_IRX_Whitaker2017(lgMstar, SFR):
    # Whitaker2017 Eq. 1 on page 3
    a = 1.96e9
    b = -2.277
    f_obscured = 1.0 / (1.0 + a * np.exp(b * lgMstar)) # f_obscured = SFR(IR)/SFR(UV+IR)
    # Whitaker2017 Eq. 2 on page 5
    #zlo = np.array([0.02, 0.5, 1.0, 1.5, 2.0])
    #zhi = np.array([0.05, 1.0, 1.5, 2.0, 2.5])
    #a = np.array([1.234, 2.092, 3.917, 6.806, 10.701])
    #b = np.array([-2.858, -2.384, -2.596, -2.673, -2.516])
    #SFR_total = SFR # SFR_total is SFR(UV+IR)
    #f_obscured = 1.0 / (1.0 + a * np.exp(b * np.log10(SFR_total))) # f_obscured = SFR(IR)/SFR(UV+IR)
    # 
    SFR_total = SFR # SFR_total is SFR(UV+IR)
    SFR_IR = SFR_total * f_obscured
    SFR_UV = SFR_total * (1.0-f_obscured)
    IRX = np.log10(SFR_IR/SFR_UV)
    # 
    return IRX, f_obscured, SFR_UV, SFR_IR

def calc_IRX_Schreiber2017(z, lgMstar, SFR):
    # Schreiber2017 Eq. 13 on page 8
    IRX = (0.45 * z + 0.35) * (lgMstar-10.5) + 1.2 # log10(LIR/LUV)
    if np.isscalar(IRX):
        is_scalar = True
    if is_scalar:
        IRX = np.array([IRX])
    # 
    IRX[(z>3.0)] = (0.45 * 3.0 + 0.35) * (lgMstar[(z>3.0)]-10.5) + 1.2 # log10(LIR/LUV)
    if is_scalar:
        IRX = IRX[0]
    # 
    f_obscured = 1.0/(1.0+(1.0/(10**IRX)))
    # 
    SFR_total = SFR
    SFR_IR = SFR_total * 1.0/(1.0+(1.0/(10**IRX)))
    SFR_UV = SFR_total * 1.0/(1.0+(10**IRX))
    IRX = np.log10(SFR_IR/SFR_UV)
    # 
    return IRX, f_obscured, SFR_UV, SFR_IR







