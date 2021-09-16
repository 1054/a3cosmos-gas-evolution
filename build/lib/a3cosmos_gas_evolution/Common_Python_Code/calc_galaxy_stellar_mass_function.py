#!/usr/bin/env python
# 
# 20190222
#     copied from "calc_stellar_mass_function.py", this code will superceed "calc_stellar_mass_function.py". 
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
    #           = (Phi_*) * 10**((lgM-lgM_*)*(alpha+1)) * exp(-10**(lgM-lgM_*)) * ln(10)  dlgx
    #           = (Our_Phi_Phi_Schechter)                                                 dlgx
    # 
    lgx = lgM-lg_M0
    Phi_Schechter = phi * (10**(lgx*(alpha+1))) * (np.exp(-10**lgx)) * ln(10) # per dex and already multiplied ln(10), so that its integral directly equals \int Phi(M) / M dM
    return Phi_Schechter




# 
# def 
# 
def calc_SMF_Davidzon2017(z, lgMstar=None, galaxy_type = 'SFG'):
    # 
    # Davidzon 2017 arXiv:1701.02734
    # IMF: Chabrier 2003
    # Outputs: lgMstar_grid, lgPhiMstar_grid
    # 
    # check z
    if not np.isscalar(z):
        raise ValueError('Please input a float number as the redshift!')
    # 
    # check galaxy_type
    if not (type(galaxy_type) is str):
        raise ValueError('Please input either "ALL", "SFG" or "QG" as the galaxy_type!')
    else:
        if not (galaxy_type in ['ALL', 'SFG', 'QG']):
            raise ValueError('Please input either "ALL", "SFG" or "QG" as the galaxy_type!')
    # 
    # make lgMstar
    if lgMstar is None:
        lgMstar_grid = np.linspace(6.0, 13.0, num=1000, endpoint=True)
    else:
        lgMstar_grid = lgMstar
    # 
    # read SMF
    tb_SMF = Table.read(os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables/datatables_SMF/datatable_Davidzon2017_SMF_'+galaxy_type+'.txt', format='ascii') # zLo zHi lgMchar  Phi_1   alpha_1  Phi_2  alpha_2
    SMF_zmin = np.min(tb_SMF['zLo'])
    SMF_zmax = np.max(tb_SMF['zHi'])
    # 
    # check z
    if z < SMF_zmin or z > SMF_zmax:
        raise ValueError('calc_SMF_Davidzon2017: The input redshift is out of the allowed range of %s -- %s!'%(SMF_zmin, SMF_zmax))
    # 
    # spline SMF #<20190214># old method
    #<20190214># SMF_z = (tb_SMF['zLo'].data + tb_SMF['zHi'].data) / 2.0
    #<20190214># SMF_phi_1 = InterpolatedUnivariateSpline(SMF_z, tb_SMF['Phi_1'].data, k=1)(z)
    #<20190214># SMF_phi_2 = InterpolatedUnivariateSpline(SMF_z, tb_SMF['Phi_2'].data, k=1)(z)
    #<20190214># SMF_alpha_1 = InterpolatedUnivariateSpline(SMF_z, tb_SMF['alpha_1'].data, k=1)(z)
    #<20190214># SMF_alpha_2 = InterpolatedUnivariateSpline(SMF_z, tb_SMF['alpha_2'].data, k=1)(z)
    #<20190214># SMF_lgMchar = InterpolatedUnivariateSpline(SMF_z, tb_SMF['lgMchar'].data, k=1)(z)
    #<20190214># #print('z, lgMchar, alpha_1, phi_1, alpha_2, phi_2 =', z, SMF_lgMchar, SMF_alpha_1, SMF_phi_1, SMF_alpha_2, SMF_phi_2)
    #<20190214># SMF_PhiMstar = Schechter_Function(lgMstar_grid, SMF_phi_1, SMF_lgMchar, SMF_alpha_1) + \
    #<20190214>#                Schechter_Function(lgMstar_grid, SMF_phi_2, SMF_lgMchar, SMF_alpha_2) # two component
    #<20190214># lgPhiMstar_grid = np.log10(SMF_PhiMstar)
    # 
    # spline SMF
    lgPhiMstar_matrix = []
    for k in range(len(tb_SMF)):
        SMF_z = (tb_SMF['zLo'][k] + tb_SMF['zHi'][k]) / 2.0
        SMF_phi_1 = tb_SMF['Phi_1'][k]
        SMF_phi_2 = tb_SMF['Phi_2'][k]
        SMF_alpha_1 = tb_SMF['alpha_1'][k]
        SMF_alpha_2 = tb_SMF['alpha_2'][k]
        SMF_lgMchar = tb_SMF['lgMchar'][k]
        #print('z, lgMchar, alpha_1, phi_1, alpha_2, phi_2 =', z, SMF_lgMchar, SMF_alpha_1, SMF_phi_1, SMF_alpha_2, SMF_phi_2)
        SMF_PhiMstar = Schechter_Function(lgMstar_grid, SMF_phi_1, SMF_lgMchar, SMF_alpha_1) + \
                       Schechter_Function(lgMstar_grid, SMF_phi_2, SMF_lgMchar, SMF_alpha_2) # two component
        lgPhiMstar_grid = np.log10(SMF_PhiMstar)
        lgPhiMstar_matrix.append(copy(lgPhiMstar_grid))
    # 
    SMF_z = (tb_SMF['zLo'].data + tb_SMF['zHi'].data) / 2.0
    lgPhiMstar_matrix = np.array(lgPhiMstar_matrix) # shape == (N_SMF_z, N_SMF_lgMstar, )
    if z <= np.min(SMF_z):
        lgPhiMstar_grid = lgPhiMstar_matrix[0]
    elif z >= np.max(SMF_z):
        lgPhiMstar_grid = lgPhiMstar_matrix[-1]
    else:
        lgPhiMstar_grid = interp1d(SMF_z, lgPhiMstar_matrix, axis=0, kind='linear')(z)
    #print(lgPhiMstar_matrix.shape, SMF_z.shape, lgPhiMstar_grid.shape)
    # fix nan
    lgPhiMstar_grid[np.isnan(lgPhiMstar_grid)] = -100
    lgPhiMstar_grid[(lgPhiMstar_grid<-100)] = -100
    # 
    if lgMstar is None:
        return lgMstar_grid, lgPhiMstar_grid
    else:
        return lgPhiMstar_grid



def calc_SMF_Ilbert2013(z, lgMstar=None, galaxy_type = 'SFG'):
    # 
    # Ilbert 2013
    # IMF: Chabrier 2003
    # Outputs: lgMstar_grid, lgPhiMstar_grid
    # 
    # check z
    if not np.isscalar(z):
        raise ValueError('Please input a float number as the redshift!')
    # 
    # check galaxy_type
    if not (type(galaxy_type) is str):
        raise ValueError('Please input either "ALL" or "SFG" as the galaxy_type!')
    else:
        if not (galaxy_type in ['ALL', 'SFG', 'QG']):
            raise ValueError('Please input either "ALL" or "SFG" as the galaxy_type!')
    # 
    # make lgMstar
    if lgMstar is None:
        lgMstar_grid = np.linspace(6.0, 13.0, num=1000, endpoint=True)
    else:
        lgMstar_grid = lgMstar
    # 
    # read SMF
    tb_SMF = Table.read(os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables/datatables_SMF/datatable_Ilbert2013_SMF_'+galaxy_type+'.txt', format='ascii') # zLo zHi lgMchar  Phi_1   alpha_1  Phi_2  alpha_2
    SMF_zmin = np.min(tb_SMF['zLo'])
    SMF_zmax = np.max(tb_SMF['zHi'])
    # 
    # check z
    if z < SMF_zmin or z > SMF_zmax:
        raise ValueError('calc_SMF_Ilbert2013: The input redshift is out of the allowed range of %s -- %s!'%(SMF_zmin, SMF_zmax))
    # 
    # spline SMF
    lgPhiMstar_matrix = []
    for k in range(len(tb_SMF)):
        SMF_z = (tb_SMF['zLo'][k] + tb_SMF['zHi'][k]) / 2.0
        SMF_phi_1 = tb_SMF['Phi_1'][k]
        SMF_phi_2 = tb_SMF['Phi_2'][k]
        SMF_alpha_1 = tb_SMF['alpha_1'][k]
        SMF_alpha_2 = tb_SMF['alpha_2'][k]
        SMF_lgMchar = tb_SMF['lgMchar'][k]
        #print('z, lgMchar, alpha_1, phi_1, alpha_2, phi_2 =', z, SMF_lgMchar, SMF_alpha_1, SMF_phi_1, SMF_alpha_2, SMF_phi_2)
        SMF_PhiMstar = Schechter_Function(lgMstar_grid, SMF_phi_1, SMF_lgMchar, SMF_alpha_1) + \
                       Schechter_Function(lgMstar_grid, SMF_phi_2, SMF_lgMchar, SMF_alpha_2) # two component
        lgPhiMstar_grid = np.log10(SMF_PhiMstar)
        lgPhiMstar_matrix.append(copy(lgPhiMstar_grid))
    # 
    SMF_z = (tb_SMF['zLo'].data + tb_SMF['zHi'].data) / 2.0
    lgPhiMstar_matrix = np.array(lgPhiMstar_matrix) # shape == (N_SMF_z, N_SMF_lgMstar, )
    if z <= np.min(SMF_z):
        lgPhiMstar_grid = lgPhiMstar_matrix[0]
    elif z >= np.max(SMF_z):
        lgPhiMstar_grid = lgPhiMstar_matrix[-1]
    else:
        lgPhiMstar_grid = interp1d(SMF_z, lgPhiMstar_matrix, axis=0, kind='linear')(z)
    #print(lgPhiMstar_matrix.shape, SMF_z.shape, lgPhiMstar_grid.shape)
    # fix nan
    lgPhiMstar_grid[np.isnan(lgPhiMstar_grid)] = -100
    lgPhiMstar_grid[(lgPhiMstar_grid<-100)] = -100
    # 
    if lgMstar is None:
        return lgMstar_grid, lgPhiMstar_grid
    else:
        return lgPhiMstar_grid



def calc_SMF_Peng2010(z, lgMstar=None, galaxy_type='SFG'):
    # 
    # Peng YingJie 2010
    # IMF: Chabrier 2003
    # Outputs: lgMstar_grid, lgPhiMstar_grid
    # 
    # check z
    if not np.isscalar(z):
        raise ValueError('Please input a float number as the redshift!')
    # 
    # check galaxy_type
    if not (type(galaxy_type) is str):
        raise ValueError('Please input either "ALL", "SFG" or "QG" as the galaxy_type!')
    else:
        if not (galaxy_type in ['SFG', 'QG']):
            raise ValueError('Please input either "ALL", "SFG" or "QG" as the galaxy_type!')
    # 
    # make lgMstar
    if lgMstar is None:
        lgMstar_grid = np.linspace(6.0, 13.0, num=1000, endpoint=True)
    else:
        lgMstar_grid = lgMstar
    # 
    # read SMF
    tb_SMF = Table.read(os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables/datatables_SMF/datatable_PengYingjie2010_SMF_'+galaxy_type+'.txt', format='ascii') # zLo zHi lgMchar  Phi_1   alpha_1  Phi_2  alpha_2
    SMF_zmin = np.min(tb_SMF['zLo'])
    SMF_zmax = np.max(tb_SMF['zHi'])
    # 
    # there is only one redshift bin, but we still check the input z range
    if z < SMF_zmin or z > SMF_zmax:
        raise ValueError('calc_SMF_Peng2010: The input redshift is out of the allowed range of %s -- %s!'%(SMF_zmin, SMF_zmax))
    # 
    # just calculate SMF without interpolation
    SMF_z = (tb_SMF['zLo'].data + tb_SMF['zHi'].data) / 2.0
    SMF_phi_1 = tb_SMF['Phi_1'].data
    SMF_alpha_1 = tb_SMF['alpha_1'].data
    SMF_lgMchar = tb_SMF['lgMchar'].data
    SMF_PhiMstar = Schechter_Function(lgMstar_grid, SMF_phi_1, SMF_lgMchar, SMF_alpha_1) # one component
    if galaxy_type == 'SFG':
        SMF_PhiMstar_SFG = copy(SMF_PhiMstar)
    else:
        SMF_phi_2 = tb_SMF['Phi_2'].data
        SMF_alpha_2 = tb_SMF['alpha_2'].data
        SMF_PhiMstar_SFG = copy(SMF_PhiMstar)
        SMF_PhiMstar_QG = SMF_PhiMstar_SFG + Schechter_Function(lgMstar_grid, SMF_phi_2, SMF_lgMchar, SMF_alpha_2) # two component QG SMF
        SMF_PhiMstar_ALL = SMF_PhiMstar_QG + SMF_PhiMstar_SFG
    # 
    if galaxy_type == 'SFG':
        lgPhiMstar_grid = np.log10(SMF_PhiMstar_SFG)
    elif galaxy_type == 'QG':
        lgPhiMstar_grid = np.log10(SMF_PhiMstar_QG)
    elif galaxy_type == 'ALL':
        lgPhiMstar_grid = np.log10(SMF_PhiMstar_ALL)
    # 
    if lgMstar is None:
        return lgMstar_grid, lgPhiMstar_grid
    else:
        return lgPhiMstar_grid



def calc_SMF_Wright2018_single_component(z, lgMstar=None):
    # 
    # Wright 2018 - SMF - arXiv.1808.01754 - 2018MNRAS.480.3491W.pdf
    # IMF: Chabrier 2003
    # Outputs: lgMstar_grid, lgPhiMstar_grid
    # 
    # check z
    if not np.isscalar(z):
        raise ValueError('Please input a float number as the redshift!')
    # 
    # make lgMstar
    if lgMstar is None:
        lgMstar_grid = np.linspace(6.0, 13.0, num=1000, endpoint=True)
    else:
        lgMstar_grid = lgMstar
    # 
    # 
    SMF_lgMchar = (10.791) + (0.558) * z + (-0.102) * z**2
    SMF_alpha = (-1.160) + (-0.274) * z + (0.028) * z**2
    SMF_phi = 10**((-2.455) + (-0.883) * z + (0.093) * z**2)
    #SMF_rho = (-2.455) + (-0.883) * z + (0.093) * z**2
    #print('z, lgMchar, alpha, phi =', z, SMF_lgMchar, SMF_alpha, SMF_phi)
    SMF_PhiMstar = Schechter_Function(lgMstar_grid, SMF_phi, SMF_lgMchar, SMF_alpha)
    lgPhiMstar_grid = np.log10(SMF_PhiMstar)
    # 
    if lgMstar is None:
        return lgMstar_grid, lgPhiMstar_grid
    else:
        return lgPhiMstar_grid


def calc_SMF_Wright2018_double_component(z, lgMstar=None):
    # 
    # Wright 2018 - SMF - arXiv.1808.01754 - 2018MNRAS.480.3491W.pdf
    # IMF: Chabrier 2003
    # Outputs: lgMstar_grid, lgPhiMstar_grid
    # 
    # check z
    if not np.isscalar(z):
        raise ValueError('Please input a float number as the redshift!')
    # 
    # make lgMstar
    if lgMstar is None:
        lgMstar_grid = np.linspace(6.0, 13.0, num=1000, endpoint=True)
    else:
        lgMstar_grid = lgMstar
    # 
    # 
    SMF_lgMchar = (10.831) + (0.153) * z + (-0.033) * z**2
    SMF_alpha_1 = (-0.579) + (0.048) * z + (0.022) * z**2
    SMF_alpha_2 = (-1.489) + (-0.087) * z + (0.016) * z**2
    SMF_phi_1 = 10**((-2.312) + (-0.658) * z + (0.016) * z**2)
    SMF_phi_2 = 10**((-3.326) + (-0.158) * z + (-0.002) * z**2)
    #print('z, lgMchar, alpha_1, phi_1, alpha_2, phi_2 =', z, SMF_lgMchar, SMF_alpha_1, SMF_phi_1, SMF_alpha_2, SMF_phi_2)
    SMF_PhiMstar = Schechter_Function(lgMstar_grid, SMF_phi_1, SMF_lgMchar, SMF_alpha_1) + \
                   Schechter_Function(lgMstar_grid, SMF_phi_2, SMF_lgMchar, SMF_alpha_2)
    lgPhiMstar_grid = np.log10(SMF_PhiMstar)
    # 
    if lgMstar is None:
        return lgMstar_grid, lgPhiMstar_grid
    else:
        return lgPhiMstar_grid









def calc_Gladders2013_CSFRD(z, t0, tau):
    Gladders2013_t_age = cosmo.age(z).value
    Gladders2013_t0 = float(t0)
    Gladders2013_tau = float(tau)
    Gladders2013_SFR_1 = 1.0/(Gladders2013_t_age*sqrt(2*pi*Gladders2013_tau**2))
    Gladders2013_SFR_2 = exp(-(ln(Gladders2013_t_age)-Gladders2013_t0)**2/(2*Gladders2013_tau**2))
    Gladders2013_SFR = Gladders2013_SFR_1 * Gladders2013_SFR_2
    Gladders2013_SFR = Gladders2013_SFR / 1.64 # converted to Chabrier IMF
    return Gladders2013_SFR


def calc_MadauDickinson2014_CSFRD(z):
    # Madau & Dickinson (2014)
    # Salpeter IMF --> converted to Chabrier IMF
    if type(z) is list:
        z = np.array(z)
    rho_SFR = 0.015 * (1+z)**2.7 / (1.0 + ((1+z)/2.9)**5.6) / 1.64 # converted to Chabrier IMF
    return rho_SFR


def calc_Mstar_integrating_CSFRD_dzliu2018(z):
    # make z_list
    #z_list = np.arange(9.75, 0.75, -0.5).tolist()
    #z_list.extend([0.75, 0.50, 0.25, 0.125, 0.0625, 0.03125, 0.00000])
    opz_list = np.logspace(np.log10(1.0+0.0), np.log10(1.0+10.75), num=200, endpoint=True)
    opz_list = opz_list[::-1]
    z_list = opz_list - 1.0
    t_list = cosmo.age(z_list).value
    # 
    # integrate dzliu2018 CSFRD
    #CSFRD = calc_Gladders2013_CSFRD(z_list, 1.57, 0.67) # Gladder2013 form of CSFRD with Liu et al. (2018) parameters
    CSFRD = calc_MadauDickinson2014_CSFRD(z_list)
    Mstar_cumulated = 0.0
    CSFRD_z_list = [z_list[0]]
    CSFRD_Mstar_list = [1e-30] #<TODO># initial mass
    for i in range(len(z_list)-1):
        #z_bin = (z_list[i+1] + z_list[i]) / 2.0
        t_bin = t_list[i+1] # Gyr, time since the Big Bang
        time_bin = t_list[i+1] - t_list[i] # Gyr, time interval within the redshift bin
        #mass_loss_time_scale = 1.4 # Myr, BEHROOZI & Silk arXiv:1404.5299
        #mass_loss_time_scale = 3.0 # Myr, Ilbert et al. (2013) PDF page 11 left middle
        mass_loss_time_scale = 0.3 # Myr, Conroy & Wechsler (2009, bibcode 2009ApJ...696..620C) arxiv PDF page 5 Eq (11). 
        Mstar_formed = (CSFRD[i]+CSFRD[i+1])/2.0 * time_bin * 1e9 # Msun/yr * Gyr * 1e9 => Msun
        Mstar_loss_frac = 0.05 * ln(1.0+(t_bin)/(mass_loss_time_scale*1e-3)) 
                                # see Ilbert et al. (2013) PDF page 11 left middle; Conroy & Wechsler (2009) arxiv PDF page 5 Eq (11). 
                                # see https://arxiv.org/pdf/1404.5299.pdf PDF page 3 Eq (6); Conroy & Wechsler (2009) arxiv PDF page 5 Eq (11). 
        #Mstar_loss_frac = 0.0 # no mass loss at all <TODO>
        Mstar_cumulated += Mstar_formed * (1.0 - Mstar_loss_frac) # total new star created, locked in stellar mass
        #CSFRD_z_list.append(z_bin)
        CSFRD_z_list.append(z_list[i+1])
        CSFRD_Mstar_list.append(Mstar_cumulated)
    CSFRD_z_list = np.array(CSFRD_z_list)[::-1] # sort so that z increases
    CSFRD_Mstar_list = np.array(CSFRD_Mstar_list)[::-1]
    # spline at z for the output
    Mstar_cumulated_at_z = 10**(InterpolatedUnivariateSpline(CSFRD_z_list, np.log10(CSFRD_Mstar_list), k=1)(z))
    return Mstar_cumulated_at_z



def calc_SMF_dzliu2018(z=None, lgMstar=None, galaxy_type='SFG', z_list=None, tuning_params='', verbose=True):
    # 
    # dzliu 2018 - my own choice
    # IMF: Chabrier 2003
    # Outputs: lgMstar_grid, lgPhiMstar_grid
    # 
    # tuning_params: for example, 'D17-no-renorm;'
    # 
    # 
    # check z
    if z is not None:
        if not np.isscalar(z):
            #raise ValueError('Please input a float number as the redshift!')
            if type(z) is list:
                z = np.array(z)
            z_is_vector = True
        else:
            z_is_vector = False
    # 
    # check galaxy_type
    if not (type(galaxy_type) is str):
        raise ValueError('Please input either "ALL", "SFG" or "QG" as the galaxy_type!')
    else:
        if not (galaxy_type in ['SFG', 'QG', 'ALL']):
            raise ValueError('Please input either "ALL", "SFG" or "QG" as the galaxy_type!')
    # 
    # make lgMstar
    if lgMstar is None:
        lgMstar_grid = np.linspace(6.0, 13.0, num=1000, endpoint=True)
    else:
        lgMstar_grid = lgMstar
    # 
    # code from '/Volumes/GoogleDrive/Team Drives/DeepFields/Simulations/Cosmological_Galaxy_Modelling/a_dzliu_code_Plot_SMF_dzliu_model.sm'
    # make z_list
    if z_list is None:
        #z_list = np.arange(10.75, 0.75, -0.5).tolist()
        #z_list.extend([0.75, 0.50, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.00000])
        z_list = np.arange(10.75, 1.00, -0.5).tolist()
        z_list.extend([1.00, 0.75, 0.50, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.00000])
        #z_list = []
        #z_list.extend(np.arange(9.75, 4.00, -0.5).tolist())
        #z_list.extend(np.arange(4.00, 3.00, -0.25).tolist())
        #z_list.extend(np.arange(3.00, 2.00, -0.125).tolist())
        #z_list.extend(np.arange(2.00, 1.00, -0.0625).tolist())
        #z_list.extend(np.arange(1.00, 0.00, -0.03125).tolist())
        #z_list.extend([0.00])
    SMF_z_list = []
    SMF_Phi_list = []
    Mstar_cumulated_list = calc_Mstar_integrating_CSFRD_dzliu2018(z_list)
    for i in range(len(z_list)-1):
        z_bin = (z_list[i] + z_list[i+1]) / 2.0
        t_bin = cosmo.age(z_bin).value # Gyr since the Big Bang
        Schechter_M = 10**lgMstar_grid
        Schechter_Mstep = lgMstar_grid[1] - lgMstar_grid[0]
        do_renorm_by_CSFRD_cumulated_Mstar = True
        # 
        # first construct SMF at each redshift bin
        #if False:
        #    EEE_Mchar_MQG = lg((10**10.68)*(t_bin/cosmo.age(0).value)**(0.50)) # Quiescent galaxies' SMF's first component
        #    EEE_Mchar_NQG = lg((10**10.68)*(t_bin/cosmo.age(0).value)**(0.50)) # Quiescent galaxies' SMF's second component
        #    EEE_Mchar_SFG = lg((10**10.72)*(t_bin/cosmo.age(0).value)**(-0.50)) # Star-forming galaxies' SMF
        #    EEE_Phi_MQG = (3.400e-3)*(t_bin/cosmo.age(0).value)**(4.00)
        #    EEE_Phi_NQG = (0.126e-3)*(t_bin/cosmo.age(0).value)**(4.00) # decrease much faster with increasing z
        #    EEE_Phi_SFG = (0.900e-3)*(t_bin/cosmo.age(0).value)**(0.20)
        #    EEE_alpha_MQG = -0.39 + (z_bin)*(1.00)
        #    EEE_alpha_NQG = -1.56 + (z_bin)*(1.00)
        #    EEE_alpha_SFG = -1.40 + (z_bin)*(-0.06) # slope -- this makes too many low mass galaxiest at z>4 -- faint-end of the 24um number counts -- no, this is because stellar SED has too old age, f24um is underestimated!
        #    #EEE_alpha_SFG = -1.40 + (z)*(-0.03) # slope
        #    Schechter_P_MQG = Schechter_Function(lgMstar_grid, EEE_Phi_MQG, EEE_Mchar_MQG, EEE_alpha_MQG)
        #    Schechter_P_NQG = Schechter_Function(lgMstar_grid, EEE_Phi_NQG, EEE_Mchar_NQG, EEE_alpha_NQG)
        #    Schechter_P_SFG = Schechter_Function(lgMstar_grid, EEE_Phi_SFG, EEE_Mchar_SFG, EEE_alpha_SFG)
        #    Schechter_P_QG = Schechter_P_MQG + Schechter_P_NQG
        #    #print('z = %.04f, lgSchechter_P_SFG = %s, lgSchechter_P_QG = %s'%(z_bin, np.log10(Schechter_P_SFG), np.log10(Schechter_P_QG)))
        # 
        if z_bin < 0.02:
            Schechter_P_QG = 10**(calc_SMF_Peng2010(0.02, lgMstar=lgMstar_grid, galaxy_type='QG'))
            Schechter_P_SFG = 10**(calc_SMF_Peng2010(0.02, lgMstar=lgMstar_grid, galaxy_type='SFG'))
        elif z_bin < 0.085:
            Schechter_P_QG = 10**(calc_SMF_Peng2010(z_bin, lgMstar=lgMstar_grid, galaxy_type='QG'))
            Schechter_P_SFG = 10**(calc_SMF_Peng2010(z_bin, lgMstar=lgMstar_grid, galaxy_type='SFG'))
        elif z_bin < 0.2:
            Schechter_P_QG = 10**(calc_SMF_Peng2010(0.085, lgMstar=lgMstar_grid, galaxy_type='QG')) * (1.0-((0.2-0.085)-(0.2-z_bin))/(0.2-0.085)) + 10**(calc_SMF_Davidzon2017(0.2, lgMstar=lgMstar_grid, galaxy_type='QG')) * (0.0+((0.2-0.085)-(0.2-z_bin))/(0.2-0.085))
            Schechter_P_SFG = 10**(calc_SMF_Peng2010(0.085, lgMstar=lgMstar_grid, galaxy_type='SFG')) * (1.0-((0.2-0.085)-(0.2-z_bin))/(0.2-0.085)) + 10**(calc_SMF_Davidzon2017(0.2, lgMstar=lgMstar_grid, galaxy_type='SFG')) * (0.0+((0.2-0.085)-(0.2-z_bin))/(0.2-0.085))
        elif z_bin < 4.0:
            Schechter_P_QG = 10**(calc_SMF_Davidzon2017(z_bin, lgMstar=lgMstar_grid, galaxy_type='QG'))
            Schechter_P_SFG = 10**(calc_SMF_Davidzon2017(z_bin, lgMstar=lgMstar_grid, galaxy_type='SFG'))
            #<TODO># QG/SFG fraction seems lower at z~1-3?
            #if z_bin > 1.0 and z_bin < 3.0:
            #    Schechter_P_QG = Schechter_P_QG * np.interp(z_bin, [1.0, 1.5, 2.0, 2.5, 3.0], [1.0, 1.2, 1.2, 1.2, 1.0])
            if tuning_params is not None:
                if tuning_params.find('D17-no-renorm')>=0:
                    do_renorm_by_CSFRD_cumulated_Mstar = False
        else:
            Schechter_P_QG = 10**(calc_SMF_Davidzon2017(4.0, lgMstar=lgMstar_grid, galaxy_type='QG'))
            Schechter_P_SFG = 10**(calc_SMF_Davidzon2017(4.0, lgMstar=lgMstar_grid, galaxy_type='SFG'))
            # 
            # Note that my Schechter_Function already contains * ln(10), i.e., it is per dex!
            # 
        # 
        # 
        # then re-normalize SMF to total stellar mass integrated from CSFRD (assumed some mass loss, see the called function)
        #Mstar_cumulated = calc_Mstar_integrating_CSFRD_dzliu2018(z_bin)
        #Mstar_cumulated = calc_Mstar_integrating_CSFRD_dzliu2018(z_list[i+1]) #<20190915># BUGGY, should be z bin edge, and z_list is in descending order 
        Mstar_cumulated = Mstar_cumulated_list[i+1]
        #CSFRD_at_z_bin = calc_MadauDickinson2014_CSFRD(z_bin)
        #<TODO><20191001># mask = (Schechter_M>=1e9) #<20191001># to match the Madau & Dickinson (2014) integration limit of '0.03 * L_characteristic', i.e., 0.03 * M_characteristic = 0.03 * 10**10.5 = 
        #<TODO><20191001># Schechter_M_total = sum((Schechter_P_QG[mask]+Schechter_P_SFG[mask])*Schechter_M[mask]*Schechter_Mstep) # P is per dex, but our Schechter_Function already contains ln(10), hence we do not need to multiply ln(10) here. \int P(M) dM = \int P(M) M dln(M) = \int P(M)*ln(10) M dlg(M)
        
        if do_renorm_by_CSFRD_cumulated_Mstar:
            Schechter_M_total = sum((Schechter_P_QG+Schechter_P_SFG)*Schechter_M*Schechter_Mstep) # P is per dex, but our Schechter_Function already contains ln(10), hence we do not need to multiply ln(10) here. \int P(M) dM = \int P(M) M dln(M) = \int P(M)*ln(10) M dlg(M)
            renorm_factor = Mstar_cumulated / Schechter_M_total
            Schechter_P_SFG = Schechter_P_SFG * renorm_factor # per dex but already contains ln(10)
            Schechter_P_QG = Schechter_P_QG * renorm_factor # per dex but already contains ln(10)
        Schechter_P_ALL = Schechter_P_SFG + Schechter_P_QG
        if verbose:
            print('z = %.04f, lgMstar_CSFRD = %0.2f, lgMstar_SMF = %0.2f, renorm = %s'%(z_bin, np.log10(Mstar_cumulated), np.log10(Schechter_M_total), renorm_factor))
        #print('z = %.04f, lgCSFRD = %0.2f, lgMstar_CSFRD = %0.2f, lgMstar_SMF = %0.2f'%(z_bin, np.log10(CSFRD_at_z_bin), np.log10(Mstar_cumulated), np.log10(Schechter_M_total)))
        # 
        # --> checked OK
        # z = 0.0156, lgCSFRD = -2.04, lgMstar_CSFRD = 8.40, lgMstar_SMF = 8.35
        # 
        # 
        SMF_z_list.append(z_list[i+1]) # append the lower redshift end point of each bin
        if galaxy_type == 'SFG':
            SMF_Phi_list.append(Schechter_P_SFG)
        elif galaxy_type == 'QG':
            SMF_Phi_list.append(Schechter_P_QG)
        elif galaxy_type == 'ALL':
            SMF_Phi_list.append(Schechter_P_ALL)
        # 
    # spline at z for the output
    SMF_z_list = np.array(SMF_z_list)[::-1] # make z increasing order
    SMF_Phi_list = np.array(SMF_Phi_list)[::-1].T # make z increasing order
    #print(SMF_z_list.shape, SMF_Phi_list.shape)
    # 
    if z is None:
        lgPhiMstar_matrix = np.log10(SMF_Phi_list.T)
        return SMF_z_list, lgMstar_grid, lgPhiMstar_matrix
    elif z_is_vector:
        #print('calc_SMF_dzliu2018: np.min(SMF_z_list), np.max(SMF_z_list), z:', np.min(SMF_z_list), np.max(SMF_z_list), z)
        lgPhiMstar_grid = interp1d(SMF_z_list, np.log10(SMF_Phi_list), kind='cubic')(z) # kind='nearest' 'linear' 'quadratic' 'cubic'
    else:
        #print('calc_SMF_dzliu2018: np.min(SMF_z_list), np.max(SMF_z_list), z:', np.min(SMF_z_list), np.max(SMF_z_list), z)
        lgPhiMstar_grid = interp1d(SMF_z_list, np.log10(SMF_Phi_list), kind='cubic')(z) # kind='nearest' 'linear' 'quadratic' 'cubic'
    # 
    if lgMstar is None:
        return lgMstar_grid, lgPhiMstar_grid
    else:
        return lgPhiMstar_grid






















