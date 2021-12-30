#!/usr/bin/env python
# 

from __future__ import print_function

import os, sys, re, json, time, astropy
import numpy as np
#from astropy.table import Table, Column, hstack
#from copy import copy
#from numpy import log10, power as pow

import astropy
from astropy import units as u
from astropy import constants as const
from packaging import version
if version.parse(astropy.__version__) < version.parse('4.0.0'):
    from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
    # funtions deprecated since v4.0.0, see -- https://astropy-astrofrog.readthedocs.io/en/latest/modeling/blackbody_deprecated.html
else:
    from astropy.modeling.models import BlackBody
    blackbody_lambda = lambda w, t: BlackBody(t, scale=1.0 * u.erg / (u.cm ** 2 * u.AA * u.s * u.sr))(w)
    blackbody_nu = lambda w, t: BlackBody(t)(w)

if sys.version_info.major >= 3:
    long = int
else:
    pass





if __name__ == '__main__':
    
    if len(sys.argv) <= 1:
        print('Usage: ')
        print('calc_dust_planck_law.py lambda_um T_dust[=25.0K]')
        print('')
        sys.exit()
    
    lambda_um = float(sys.argv[1])
    if len(sys.argv) > 2:
        T_dust = float(sys.argv[2])
    else:
        T_dust = 25.0
    
    #beta = 2.0
    #beta = 1.8
    
    print('T_dust = %s [K]'%(T_dust))
    
    # calc Planck law
    cal_Lv_obs_erg_s_Hz = blackbody_nu(lambda_um * u.um, T_dust * u.K) # * np.power(2.99792458e5/lambda_um, beta)
    print(cal_Lv_obs_erg_s_Hz)
    #print('%0.6e'%(cal_Lv_obs_erg_s_Hz))
    
    
    # calc Rayleigh-Jeans approximation
    #cal_Lv_obs_erg_s_Hz_2 = 2 * const.k_B.to('erg/K').value * (const.c.to('cm/s').value/(lambda_um/1e4))**2 / (const.c.to('cm/s').value)**2 * T_dust
    #print(cal_Lv_obs_erg_s_Hz_2)





