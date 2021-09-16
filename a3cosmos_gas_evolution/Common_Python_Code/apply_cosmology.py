#!/usr/bin/env python
# 

import sys, inspect
from astropy.cosmology import FlatLambdaCDM

def apply_cosmology(Hubble_Constant_z0 = 70, Omega_Matter = 0.27, Omega_Lambda = 0.73, T_CMB_z0 = 2.725):
    # 
    #cosmo = FlatLambdaCDM(H0=73, Om0=0.27, Tcmb0=2.725) # replaced since 2019-02-21
    # 
    # Check out the inspect module, it is used by minimock to mock the caller's scope.
    stack = inspect.stack()
    # 
    # check 'cosmo' variable
    if ('cosmo' in (inspect.stack()[-1][0]).f_globals.keys()):
        cosmo = (inspect.stack()[-1][0]).f_globals['cosmo']
    elif ('cosmo' in (inspect.stack()[-1][0]).f_locals.keys()):
        cosmo = (inspect.stack()[-1][0]).f_locals['cosmo']
    else:
        if '%0.2f'%(Omega_Lambda + Omega_Matter) == '1.00':
            cosmo = FlatLambdaCDM(H0 = Hubble_Constant_z0, Om0 = Omega_Matter, Tcmb0 = T_CMB_z0)
        else:
            raise NotImplementedError('Error! Non-flat cosmology not implemented yet! <TODO>')
    # 
    #print(cosmo.H0)
    # 
    return cosmo



#print((inspect.stack()[-1][0]))
#print((inspect.stack()[-1][0]).f_locals.keys())

#if ('cosmo' in (inspect.stack()[-1][0]).f_globals.keys()):
#    cosmo = (inspect.stack()[-1][0]).f_globals['cosmo']
#elif ('cosmo' in (inspect.stack()[-1][0]).f_locals.keys()):
#    cosmo = (inspect.stack()[-1][0]).f_locals['cosmo']
#else:
#    cosmo = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.725)

cosmo = apply_cosmology()

#print(cosmo.H0)

