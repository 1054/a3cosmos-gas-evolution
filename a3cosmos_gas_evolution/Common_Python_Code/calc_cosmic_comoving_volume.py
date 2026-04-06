#!/usr/bin/env python
#

from __future__ import print_function

import os, sys, re, json, time, astropy
import numpy as np
#from astropy.table import Table, Column, hstack
from astropy import units as u
from copy import copy

if not (os.path.dirname(os.path.abspath(__file__)) in sys.path): sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import apply_cosmology
#cosmo = apply_cosmology.cosmo
cosmo = apply_cosmology.apply_cosmology(70, 0.3, 0.7)
print('cosmo', cosmo)

if sys.version_info.major >= 3:
    long = int
else:
    pass




def calc_cosmic_comoving_volume_1(z_edge_1, z_edge_2, obs_area_arcmin2):
    if type(obs_area_arcmin2) is u.quantity.Quantity:
        obs_area = obs_area_arcmin2
    else:
        obs_area = obs_area_arcmin2 * u.arcmin**2
    #comoving_z_list = np.linspace(z_edges[i], z_edges[i+1], num=100, endpoint=True)
    #comoving_volume = np.sum((cosmo.comoving_volume(comoving_z_list[1:]) - cosmo.comoving_volume(comoving_z_list[0:-1])) / (4.0*np.pi*u.steradian) * obs_area.to(u.steradian))
    #print('comoving_volume = %e [%s]'%(comoving_volume.value, comoving_volume.unit))
    comoving_volume = ((cosmo.comoving_volume(z_edge_2) - cosmo.comoving_volume(z_edge_1)) / (4.0*np.pi*u.steradian) * obs_area.to(u.steradian))
    print('comoving_volume = %e [%s]'%(comoving_volume.value, comoving_volume.unit))
    return comoving_volume.value



def calc_cosmic_comoving_volume_2(z_edge_1, z_edge_2, obs_area_arcmin2):
    if type(obs_area_arcmin2) is u.quantity.Quantity:
        obs_area = obs_area_arcmin2
    else:
        obs_area = obs_area_arcmin2 * u.arcmin**2
    differntial_z_list = np.linspace(z_edge_1, z_edge_2, num=10, endpoint=True)
    comoving_volume = np.sum((cosmo.differential_comoving_volume(differntial_z_list[1:]) * np.diff(differntial_z_list) * obs_area.to(u.steradian)))
    print('comoving_volume = %e [%s]'%(comoving_volume.value, comoving_volume.unit))
    ##print(cosmo.de_density_scale(z)) # should be 1
    ##print(cosmo._Ogamma0, cosmo._Onu0)
    ##print(cosmo.efunc(z), np.sqrt(0.27*(1.+z)**3 + 0*(1.+z)**2 + 0.73) ) # checked consistent
    ##print(cosmo._hubble_distance, 2.997902458e5/70 ) # checked consistent
    #sys_lumdist_output = subprocess.getoutput("/Users/dzliu/Cloud/Github/Crab.Toolkit.PdBI/bin/lumdist -h0 70 -verbose %s | grep 'lumdist d_L=' | sed -e 's/=/ /g'"%(z))
    ##print(cosmo.angular_diameter_distance(z), sys_lumdist_output.split()[8], '(z = %s)'%(z), sys_lumdist_output ) # 
    #dH_astropy = cosmo._hubble_distance
    #Ez_astropy = cosmo.efunc(z)
    #dA_astropy = cosmo.angular_diameter_distance(z)
    #dH_dzliu = 2.997902458e5/70
    #Ez_dzliu = np.sqrt(0.27*(1.+z)**3 + 0*(1.+z)**2 + 0.73)
    #dA_dzliu = float(sys_lumdist_output.split()[8])
    #print(dH_astropy/Ez_astropy*dA_astropy, dH_dzliu/Ez_dzliu*dA_dzliu) # chekced consistent
    #zp1 = 1.0 + z
    #print(cosmo.differential_comoving_volume(z), dH_astropy*((zp1*dA_astropy)**2)/Ez_astropy, dH_dzliu/Ez_dzliu*dA_dzliu**2*zp1**2) # chekced consistent
    # z = 2.0 - 2.5, area = 0.0006092348395183178 steradian, dVc = 43627725623.05944, 
    # 43627725623.05944 * 0.25 / 10
    return comoving_volume.value



if __name__ == '__main__':
    
    obs_area = 1.5546582999901375*u.deg*u.deg
    print('obs_area = %s [%s]'%(obs_area.to(u.arcmin*u.arcmin).value, obs_area.to(u.arcmin*u.arcmin).unit))
    print('obs_area = %s [%s]'%(obs_area.to(u.steradian).value, obs_area.to(u.steradian).unit))
    #print(type(obs_area))
    
    z_edges = [0.02, 0.25, 0.50, 0.75, 1.00, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]
    
    # 
    # loop z bin
    for i in range(len(z_edges)-1):
        # 
        print('z %s - %s, cosmic age %.2f - %.2f, time interval %.2f'%(\
                z_edges[i], 
                z_edges[i+1], 
                cosmo.age(z_edges[i]).to('Gyr').value, 
                cosmo.age(z_edges[i+1]).to('Gyr').value, 
                cosmo.age(z_edges[i]).to('Gyr').value - cosmo.age(z_edges[i+1]).to('Gyr').value
                )
             )
        # 
        z = (z_edges[i]+z_edges[i+1])/2.0
        # 
        calc_cosmic_comoving_volume_1(z_edges[i], z_edges[i+1], obs_area)
        calc_cosmic_comoving_volume_2(z_edges[i], z_edges[i+1], obs_area)
        # 
        # 
        # 
    


