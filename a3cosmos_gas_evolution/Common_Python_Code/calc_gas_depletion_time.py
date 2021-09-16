#!/usr/bin/env python
# 

from __future__ import print_function

import os, sys, re, json, time, astropy
import numpy as np
from astropy.table import Table, Column, hstack
from numpy import log10, power as pow

#from astropy.cosmology import FlatLambdaCDM
##cosmo = FlatLambdaCDM(H0=73, Om0=0.27, Tcmb0=2.725) # replaced since 2019-02-21
#cosmo = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.725)
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
def calc_gas_depletion_time_Scoville2017(z, lgMstar=10.5, DeltaMS=0.0):
    # Scoville et al. 2017 (The Astrophysical Journal, 837:150 (20pp), 2017 March 10)
    # Table 2
    Ratio_M_molgas_SFR = 3.23 * (1+z)**(-1.04) * (10**DeltaMS)**(-0.70) * (10**(lgMstar-10.0))**(-0.01)
    return Ratio_M_molgas_SFR

def calc_gas_depletion_time_Tacconi2018(z, lgMstar=10.5, DeltaMS=0.0):
    # Tacconi et al. 2018 (The Astrophysical Journal, 853:179 (22pp), 2018 February 1 https://doi.org/10.3847/1538-4357/aaa4b4)
    # Table 3
    # Best (with bootstrap errors) with S14 (Speagle+2014) MS
    A = +0.09
    B = -0.62
    C = -0.44
    D = +0.09
    E = +0.11
    DeltaRe = 0.0 # log10 (R_e / R_e_0)
    log10_Ratio_M_molgas_SFR = A + B * (np.log10(1+z)) + C * DeltaMS + D * (lgMstar-10.7) + E * (DeltaRe)
    Ratio_M_molgas_SFR = 10**log10_Ratio_M_molgas_SFR
    return Ratio_M_molgas_SFR

def calc_gas_depletion_time_A3COSMOS(z=None, cosmic_age=None, lgMstar=10.5, DeltaMS=0.0):
    # A3COSMOS 
    # Sample 20181203 
    # fitting '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plot/Plot_z_tauDepl_data_fitting/20190716_with_gas_mass_calibration_H17/fitting_dzliu_function/best_fit_via_MCMC.json'
    #popt = [-0.5820833042317259, -0.5346887707021821, -0.0038760530734522902, 0.041705991802229114, 0.1252905060544558, 0.0592184169424943]
    popt = [-0.5816185208613138, -0.5338752260098323, -0.0038228411346192814, 0.040738011277865915, 0.12494331865991803, 0.05914760529893037]
    a, b, c, d, ak, ck = popt
    if z is None and cosmic_age is None:
        print('Error! Please input either z or cosmic_age')
        sys.exit()
    else:
        if cosmic_age is not None:
            cosmoAge = cosmic_age
        else:
            cosmoAge = cosmo.age(z).value
    log10_Ratio_M_molgas_SFR = (a+ak*(lgMstar-10.0))*DeltaMS + b*(lgMstar-10.0) + (c+ck*(lgMstar-10.0))*(cosmoAge) + d
    Ratio_M_molgas_SFR = 10**log10_Ratio_M_molgas_SFR
    return Ratio_M_molgas_SFR

#def calc_gas_depletion_time_A3COSMOS(z, lgMstar=10.5, DeltaMS=0.0):
#    # A3COSMOS 
#    # Sample 20181203 
#    # fitting '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plots/Plot_z_tauDepl_data_fitting/best_fit_via_CurveFit.json'
#    #popt = [0.1, -0.4979428049052056, -0.5172364530151924, -0.5, -0.12515630072829534] #<201901>#
#    #popt = [-0.5422, -0.4824, +0.04271, -0.01896, +0.08644] #<20190205>#
#    #popt = [-0.5817, -0.5640, 0.05553, -0.01192, 0.1375] #<20190207># fixed xCOLDGASS Z_PP04_O3N2
#    #popt = [-0.5741, -0.5664, +0.05543, -0.004397, +0.1154] #<20190208># still fitting Bertemes and Lee's data
#    popt = [-0.5702, -0.5741, 0.05489, 0.008601, -0.01334] #<20190301># Eq850 this work, RemyRuyer GDR, Genzel2015_Eq12a_with_dzliu_limit MZR, KMT09 fmol, fitting all available data '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plots/Plot_z_tauDepl_data_fitting/20190301_with_gas_mass_calibration_Eq850_dzliu/best_fit_via_MCMC.json'
#    a, b, c, d, e = popt
#    #log10_Ratio_M_molgas_SFR = b*DeltaMS + c*np.log10(1+z) + e #<201901>#
#    cosmoAge = cosmo.age(z).value
#    log10_Ratio_M_molgas_SFR = a*DeltaMS + (b+c*cosmoAge)*(lgMstar-10.0) + d*(cosmoAge) + e #<20190205>#
#    
#    #<20190714># 
#    a = -0.5788835106242605
#    b = -0.5376145163745489
#    c = -0.004112269290743598
#    d = 0.04527343126448624
#    ak = 0.12203773985943833
#    ck = 0.05930642823001797
#    cosmoAge = cosmo.age(z).value
#    log10_Ratio_M_molgas_SFR = (a+ak*(lgMstar-10.0))*DeltaMS + b*(lgMstar-10.0) + (c+ck*(lgMstar-10.0))*cosmoAge + d
#    
#    
#    Ratio_M_molgas_SFR = 10**log10_Ratio_M_molgas_SFR
#    return Ratio_M_molgas_SFR

#def calc_gas_depletion_time_A3COSMOS_with_dzliu_850_gas_mass_calibration(z, lgMstar=10.5, DeltaMS=0.0):
#    # A3COSMOS 
#    # Sample 20181203 
#    # fitting '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plots/Plot_z_tauDepl_data_fitting/20190122_with_Leroy_GDR_with_Genzel_Eq12a_MZR_with_dzliu_gas_mass_calibration_with_KMT09_fmol/best_fit_via_CurveFit.json'
#    #popt = [0.1, -0.5191087642422786, -1.5519298330096012, -0.5, 0.08064040034418785] #<20190122_with_Leroy_GDR_with_Genzel_Eq12a_MZR_with_KMT09_fmol>#
#    popt = [-0.5702, -0.5741, 0.05489, 0.008601, -0.01334] #<20190301># Eq850 this work, RemyRuyer GDR, Genzel2015_Eq12a_with_dzliu_limit MZR, KMT09 fmol, fitting all available data '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plots/Plot_z_tauDepl_data_fitting/20190301_with_gas_mass_calibration_Eq850_dzliu/best_fit_via_MCMC.json'
#    a, b, c, d, e = popt
#    log10_Ratio_M_molgas_SFR = b*DeltaMS + c*np.log10(1+z) + e
#    Ratio_M_molgas_SFR = 10**log10_Ratio_M_molgas_SFR
#    return Ratio_M_molgas_SFR

#def calc_gas_depletion_time_A3COSMOS_with_Hughes2017_GasMassCalibration_Leslie2019_MS(z, lgMstar=10.5, DeltaMS=0.0):
#    # A3COSMOS 
#    # Sample 20181203 
#    # MS: Leslie20190515
#    # GasMassCalibration: Hughes2017
#    # fitting '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plots/Plot_z_tauDepl_data_fitting/20190301_with_gas_mass_calibration_H17_with_MS_Leslie20190515/best_fit_via_MCMC.json'
#    popt = [-0.4367323964697807, -0.31995589347193487, 0.0363679934437835, 0.011782986854565358, -0.19084213551539264]
#    a, b, c, d, e = popt
#    #log10_Ratio_M_molgas_SFR = b*DeltaMS + c*np.log10(1+z) + e #<201901>#
#    cosmoAge = cosmo.age(z).value
#    log10_Ratio_M_molgas_SFR = a*DeltaMS + (b+c*cosmoAge)*(lgMstar-10.0) + d*(cosmoAge) + e #<20190205>#
#    Ratio_M_molgas_SFR = 10**log10_Ratio_M_molgas_SFR
#    return Ratio_M_molgas_SFR






# 
# def core functions (adapted from "a_dzliu_code_fit_z_tauDepl_v6.py")
# 
def func_tauDepl_dzliu_log(pars, a, b, c, d, e): 
    deltaGas,DeltaMS,lgMstar,lgSFR,cosmoAge,z = pars
    # 20181030c
    #tauDepl_model_log = a*DeltaMS + b*(deltaGas) + c*(np.exp(-cosmoAge)) + d
    #tauDepl_model_log = a*(np.log10(deltaGas))*0.0 + b*DeltaMS + c*10**lgSFR*0.0 + lgMstar*0.0 + d*(cosmoAge)*0.0 + e
    tauDepl_model_log = a*(np.log10(deltaGas))*0.0 + b*DeltaMS + 10**lgSFR*0.0 + c*np.log10(1+z) + lgMstar*0.0 + d*(cosmoAge)*0.0 + e
    # 
    return tauDepl_model_log

def func_tauDepl_dzliu(pars, a, b, c, d, e): 
    deltaGas,DeltaMS,lgMstar,lgSFR,cosmoAge,z = pars
    tauDepl_model_log = func_tauDepl_dzliu_log((deltaGas,DeltaMS,lgMstar,lgSFR,cosmoAge,z), a, b, c, d, e)
    tauDepl_model = np.power(10, tauDepl_model_log) # * np.power((1+z),c)
    return tauDepl_model

def func_tauDepl_Tacconi2018_log(pars, a, b, c, d): 
    DeltaMS,lgMstar,z = pars
    tauDepl_model_log = a*DeltaMS + b*(lgMstar-10.7) + c*(np.log10(1+z)) + d
    return tauDepl_model_log

def func_tauDepl_Tacconi2018(pars):
    DeltaMS,lgMstar,z = pars
    tauDepl_model_log = func_tauDepl_Tacconi2018_log((DeltaMS,lgMstar,z), -0.44, 0.09, -0.62, 0.09)
    tauDepl_model = np.power(10, tauDepl_model_log) # * np.power((1+z),c)
    return tauDepl_model

def func_tauDepl_Scoville2017_log(pars, a, b, c, d): 
    DeltaMS,lgMstar,z = pars
    tauDepl_model_log = a*DeltaMS + b*(lgMstar-10.0) + c*(np.log10(1+z)) + d
    return tauDepl_model_log

def func_tauDepl_Scoville2017(pars):
    DeltaMS,lgMstar,z = pars
    tauDepl_model_log = func_tauDepl_Scoville2017_log((DeltaMS,lgMstar,z), -0.70, -0.01, -1.04, np.log10(3.23))
    tauDepl_model = np.power(10, tauDepl_model_log) # * np.power((1+z),c)
    return tauDepl_model









