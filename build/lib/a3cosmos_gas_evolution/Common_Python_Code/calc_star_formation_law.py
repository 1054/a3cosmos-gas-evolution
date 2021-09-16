#!/usr/bin/env python
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
def calc_Mmolgas_from_SFR_using_SF_law_Sargent2014(SFR, DeltaMS = 0.0):
    # Sargent et al. 2014, https://ui.adsabs.harvard.edu/abs/2014ApJ...793...19S/abstract
    # bimodal SF law
    arr_alpha = [9.22, 8.05]
    arr_beta = [0.81, 0.81]
    arr_DeltaMS = [0.0, 1.0]
    argmin_DeltaMS = np.argmin(arr_DeltaMS)
    argmax_DeltaMS = np.argmax(arr_DeltaMS)
    min_DeltaMS = arr_DeltaMS[argmin_DeltaMS]
    max_DeltaMS = arr_DeltaMS[argmax_DeltaMS]
    # 
    if np.isscalar(DeltaMS):
        if DeltaMS <= min_DeltaMS:
            t_DeltaMS = min_DeltaMS
            t_alpha = arr_alpha[argmin_DeltaMS]
            t_beta = arr_beta[argmin_DeltaMS]
        elif DeltaMS >= max_DeltaMS:
            t_DeltaMS = max_DeltaMS
            t_alpha = arr_alpha[argmax_DeltaMS]
            t_beta = arr_beta[argmax_DeltaMS]
        else:
            t_DeltaMS = DeltaMS
            t_alpha = np.interp(t_DeltaMS, arr_DeltaMS, arr_alpha)
            t_beta = np.interp(t_DeltaMS, arr_DeltaMS, arr_beta)
    else:
        t_DeltaMS = np.array(DeltaMS)
        t_alpha = t_DeltaMS * 0.0
        t_beta = t_DeltaMS * 0.0
        t_mask_1 = (t_DeltaMS <= min_DeltaMS)
        if np.count_nonzero(t_mask_1) > 0:
            t_DeltaMS[t_mask_1] = min_DeltaMS
            t_alpha[t_mask_1] = arr_alpha[argmin_DeltaMS]
            t_beta[t_mask_1] = arr_beta[argmin_DeltaMS]
        t_mask_2 = (t_DeltaMS >= max_DeltaMS)
        if np.count_nonzero(t_mask_2) > 0:
            t_DeltaMS[t_mask_2] = max_DeltaMS
            t_alpha[t_mask_2] = arr_alpha[argmax_DeltaMS]
            t_beta[t_mask_2] = arr_beta[argmax_DeltaMS]
        t_mask_3 = np.logical_and(t_DeltaMS > min_DeltaMS, t_DeltaMS < max_DeltaMS)
        if np.count_nonzero(t_mask_3) > 0:
            t_alpha[t_mask_3] = np.interp(t_DeltaMS[t_mask_3], arr_DeltaMS, arr_alpha)
            t_beta[t_mask_3] = np.interp(t_DeltaMS[t_mask_3], arr_DeltaMS, arr_beta)
    Mmol = 10**(t_alpha + t_beta * np.log10(SFR))
    return Mmol


def calc_SFR_from_Mmolgas_using_SF_law_Sargent2014(Mmolgas, DeltaMS = 0.0):
    # Sargent et al. 2014, https://ui.adsabs.harvard.edu/abs/2014ApJ...793...19S/abstract
    # bimodal SF law
    arr_alpha = [9.22, 8.05]
    arr_beta = [0.81, 0.81]
    arr_DeltaMS = [0.0, 1.0]
    argmin_DeltaMS = np.argmin(arr_DeltaMS)
    argmax_DeltaMS = np.argmax(arr_DeltaMS)
    min_DeltaMS = arr_DeltaMS[argmin_DeltaMS]
    max_DeltaMS = arr_DeltaMS[argmax_DeltaMS]
    # 
    if np.isscalar(DeltaMS):
        if DeltaMS <= min_DeltaMS:
            t_DeltaMS = min_DeltaMS
            t_alpha = arr_alpha[argmin_DeltaMS]
            t_beta = arr_beta[argmin_DeltaMS]
        elif DeltaMS >= max_DeltaMS:
            t_DeltaMS = max_DeltaMS
            t_alpha = arr_alpha[argmax_DeltaMS]
            t_beta = arr_beta[argmax_DeltaMS]
        else:
            t_DeltaMS = DeltaMS
            t_alpha = np.interp(t_DeltaMS, arr_DeltaMS, arr_alpha)
            t_beta = np.interp(t_DeltaMS, arr_DeltaMS, arr_beta)
    else:
        t_DeltaMS = np.array(DeltaMS)
        t_alpha = t_DeltaMS * 0.0
        t_beta = t_DeltaMS * 0.0
        t_mask_1 = (t_DeltaMS <= min_DeltaMS)
        if np.count_nonzero(t_mask_1) > 0:
            t_DeltaMS[t_mask_1] = min_DeltaMS
            t_alpha[t_mask_1] = arr_alpha[argmin_DeltaMS]
            t_beta[t_mask_1] = arr_beta[argmin_DeltaMS]
        t_mask_2 = (t_DeltaMS >= max_DeltaMS)
        if np.count_nonzero(t_mask_2) > 0:
            t_DeltaMS[t_mask_2] = max_DeltaMS
            t_alpha[t_mask_2] = arr_alpha[argmax_DeltaMS]
            t_beta[t_mask_2] = arr_beta[argmax_DeltaMS]
        t_mask_3 = np.logical_and(t_DeltaMS > min_DeltaMS, t_DeltaMS < max_DeltaMS)
        if np.count_nonzero(t_mask_3) > 0:
            t_alpha[t_mask_3] = np.interp(t_DeltaMS[t_mask_3], arr_DeltaMS, arr_alpha)
            t_beta[t_mask_3] = np.interp(t_DeltaMS[t_mask_3], arr_DeltaMS, arr_beta)
    SFR = 10**( (np.log10(Mmolgas) - t_alpha) / t_beta )
    return SFR





