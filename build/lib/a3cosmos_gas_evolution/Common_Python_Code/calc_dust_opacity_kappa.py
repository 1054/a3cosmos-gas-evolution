#!/usr/bin/env python
# 

from __future__ import print_function

import os, sys, re, json, time, astropy
import numpy as np
#from astropy.table import Table, Column, hstack
#from copy import copy
#from numpy import log10, power as pow

#from astropy import units as u
#from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu

if sys.version_info.major >= 3:
    long = int
else:
    pass





if __name__ == '__main__':
    
    if len(sys.argv) <= 1:
        print('Usage: ')
        print('calc_dust_opacity_kappa.py lambda_um')
        print('')
        sys.exit()
    
    lambda_um = float(sys.argv[1])
    
    if lambda_um >= 700.0:
        beta = 1.68
    else:
        beta = 2.0
    
    kappa_ = 0.596 * (lambda_um / 700.0)**(-beta)
    
    print(kappa_, '[cm^2 g^{-1}]')
    #print('%0.6e'%(cal_Lv_obs_erg_s_Hz))





