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
def Schechter_Function(lgM, phi, lg_M0, alpha):
    # 
    # Schechter (1976)
    # 
    # Phi(M) dM = (Phi_*) * (M/M_*)**(alpha) * exp(-M/M_*) dM/M_*
    #           = (Phi_*) * x**(alpha) * exp(-x) dx
    #           = (Phi_*) * 10**(lgx * alpha) * exp(-10**lgx) dx
    #           = (Phi_*) * 10**(lgx * alpha) * exp(-10**lgx) dlnx
    #           = (Phi_*) * 10**(lgx * alpha) * exp(-10**lgx) dlgx * ln(10)
    #           = ln(10) * (Phi_*) * 10**((lgM-lgM_*)*(alpha+1)) * exp(-10**(lgM-lgM_*))
    # 
    lgx = lgM-lg_M0
    Phi_Schechter = phi * (10**(lgx*(alpha+1))) * (np.exp(-10**lgx))*ln(10)
    return Phi_Schechter




# 
# def 
# 
def calc_CO10_LF_Saintonge2017(lgMgas=None, input_type=1):
    # 
    # Saintonge 2017 (), Table 5 and Figure 6
    # IMF: Chabrier 2003
    # Outputs: lgMgas_grid, lgPhiMgas_grid
    # Input_type 1 means the analysis is done with detections only
    # Input_type 2 means the analysis is done with detections+nondetections
    # 
    # 
    # make lgMgas
    if lgMgas is None:
        lgMgas_grid = np.linspace(6.0, 13.0, num=1000, endpoint=True)
    else:
        lgMgas_grid = lgMgas
    # 
    # read GMF
    tb = Table.read(os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables/datatables_GMF/datatable_Saintonge2017_CO10_LF_%s.txt'%(input_type), format='ascii')
    GMF_zmin = np.min(tb['zLo'])
    GMF_zmax = np.max(tb['zHi'])
    GMF_lgMchar = tb['lgLchar'][0]
    GMF_phi_1 = tb['Phi_1'][0]
    GMF_alpha_1 = tb['alpha_1'][0]
    # 
    GMF_Phi_L_Prime_CO10 = Schechter_Function(lgMgas_grid, GMF_phi_1, GMF_lgMchar, GMF_alpha_1) # single component
    lgPhiMgas_grid = np.log10(GMF_Phi_L_Prime_CO10)
    # 
    # fix nan
    lgPhiMgas_grid[np.isnan(lgPhiMgas_grid)] = -100
    lgPhiMgas_grid[(lgPhiMgas_grid<-100)] = -100
    # 
    if lgMgas is None:
        return lgMgas_grid, lgPhiMgas_grid
    else:
        return lgPhiMgas_grid

def calc_CO10_LF_Saintonge2017_updated(lgMgas=None, input_type=1):
    # 
    # Saintonge 2017 updated CO LF from Dominik Riechers and Riccardo Pavesi (priv. comm.)
    # IMF: Chabrier 2003
    # Outputs: lgMgas_grid, lgPhiMgas_grid
    # 
    # 
    # make lgMgas
    if lgMgas is None:
        lgMgas_grid = np.linspace(6.0, 13.0, num=1000, endpoint=True)
    else:
        lgMgas_grid = lgMgas
    # 
    # read GMF
    tb = Table.read(os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables/datatables_GMF/datatable_Saintonge2017_CO10_LF_%s_updated.txt'%(input_type), format='ascii')
    GMF_zmin = np.min(tb['zLo'])
    GMF_zmax = np.max(tb['zHi'])
    GMF_lgMchar = tb['lgLchar'][0]
    GMF_phi_1 = tb['Phi_1'][0]
    GMF_alpha_1 = tb['alpha_1'][0]
    # 
    GMF_Phi_L_Prime_CO10 = Schechter_Function(lgMgas_grid, GMF_phi_1, GMF_lgMchar, GMF_alpha_1) # single component
    lgPhiMgas_grid = np.log10(GMF_Phi_L_Prime_CO10)
    # 
    # fix nan
    lgPhiMgas_grid[np.isnan(lgPhiMgas_grid)] = -100
    lgPhiMgas_grid[(lgPhiMgas_grid<-100)] = -100
    # 
    if lgMgas is None:
        return lgMgas_grid, lgPhiMgas_grid
    else:
        return lgPhiMgas_grid







