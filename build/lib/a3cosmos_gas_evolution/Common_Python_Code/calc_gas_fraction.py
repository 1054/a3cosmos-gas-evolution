#!/usr/bin/env python
# 

from __future__ import print_function

import os, sys, re, json, time, astropy, copy
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
def calc_gas_fraction_Scoville2017(z, lgMstar=10.5, DeltaMS=0.0):
    # Scoville et al. 2017 (The Astrophysical Journal, 837:150 (20pp), 2017 March 10)
    # Table 2
    Ratio_M_molgas_M_star = 0.71 * (1+z)**1.84 * (10**DeltaMS)**0.32 * (10**(lgMstar-10.0))**(-0.70)
    return Ratio_M_molgas_M_star

def calc_gas_fraction_Tacconi2018_beta_2(z, lgMstar=10.5, DeltaMS=0.0):
    # Tacconi et al. 2018 (The Astrophysical Journal, 853:179 (22pp), 2018 February 1 https://doi.org/10.3847/1538-4357/aaa4b4)
    # Table 3
    # Best^c beta=2 with S14 (Speagle+2014) MS
    A = +0.12
    B = -3.62
    F = +0.66
    beta = 2.0
    C = +0.53
    D = -0.35
    E = +0.11
    DeltaRe = 0.0 # log10 (R_e / R_e_0)
    log10_Ratio_M_molgas_M_star = A + B * (np.log10(1+z)-F)**beta + C * DeltaMS + D * (lgMstar-10.7) + E * (DeltaRe)
    Ratio_M_molgas_M_star = 10**log10_Ratio_M_molgas_M_star
    return Ratio_M_molgas_M_star

def calc_gas_fraction_Tacconi2018(z, lgMstar=10.5, DeltaMS=0.0):
    return calc_gas_fraction_Tacconi2018_beta_2(z, lgMstar, DeltaMS)

def calc_gas_fraction_A3COSMOS(cosmic_age=None, z=None, lgMstar=10.5, DeltaMS=0.0):
    # A3COSMOS 
    # Sample 20181203 
    # fitting '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plot/Plot_z_deltaGas_data_fitting/20190715_with_gas_mass_calibration_H17/fitting_dzliu_function/best_fit_via_MCMC.json'
    #popt = [0.38380575488270824, -0.8341369033600499, -0.15468695242992148, 1.0834290154064226, 0.10008727595632205, 0.040159361310305286]
    # fitting '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plot/Plot_z_deltaGas_data_fitting/20190802_with_gas_mass_calibration_H17_with_own_uncertainty/fitting_dzliu_function/best_fit_via_MCMC.json'
    #popt = [0.4993149209610621, -0.665870841347096, -0.1488184998332125, 0.8970477945369266, 0.0989348364520346, 0.02889191643165745]
    # fitting '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plot/Plot_z_deltaGas_data_fitting/20190827_with_gas_mass_calibration_H17/fitting_dzliu_function/best_fit_via_MCMC.json' # -- now with Cortzen+2019
    popt = [0.4195381171312362, -0.6906526862112656, -0.1542557899746253, 0.9339052286147904, 0.11949845461785102, 0.03204158046196248]
    a, b, c, d, ak, ck = popt
    if z is None and cosmic_age is None:
        print('Error! Please input either z or cosmic_age')
        sys.exit()
    else:
        if cosmic_age is not None:
            cosmoAge = cosmic_age
        else:
            cosmoAge = cosmo.age(z).value
    log10_Ratio_M_molgas_M_star = (a+ak*(lgMstar-10.0))*DeltaMS + b*(lgMstar-10.0) + (c+ck*(lgMstar-10.0))*(cosmoAge) + d
    Ratio_M_molgas_M_star = 10**log10_Ratio_M_molgas_M_star
    return Ratio_M_molgas_M_star

def calc_gas_fraction_A3COSMOS_with_Leslie2019_MS(cosmic_age=None, z=None, lgMstar=10.5, DeltaMS=0.0):
    # A3COSMOS 
    # Sample 20181203 
    # fitting '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plot/Plot_z_deltaGas_data_fitting/20190827_with_gas_mass_calibration_H17_and_main_sequence_Leslie2019/fitting_dzliu_function/best_fit_via_MCMC.json'
    popt = [0.3781725710882265, -0.8360054443744662, -0.1557707273000366, 1.0923711580061433, 0.09940673198827188, 0.03940816266321434]
    a, b, c, d, ak, ck = popt
    if z is None and cosmic_age is None:
        print('Error! Please input either z or cosmic_age')
        sys.exit()
    else:
        if cosmic_age is not None:
            cosmoAge = cosmic_age
        else:
            cosmoAge = cosmo.age(z).value
    log10_Ratio_M_molgas_M_star = (a+ak*(lgMstar-10.0))*DeltaMS + b*(lgMstar-10.0) + (c+ck*(lgMstar-10.0))*(cosmoAge) + d
    Ratio_M_molgas_M_star = 10**log10_Ratio_M_molgas_M_star
    return Ratio_M_molgas_M_star

def calc_gas_fraction_A3COSMOS_with_Scoville2017_MS(cosmic_age=None, z=None, lgMstar=10.5, DeltaMS=0.0):
    # A3COSMOS 
    # Sample 20181203 
    # fitting '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plot/Plot_z_deltaGas_data_fitting/20190827_with_gas_mass_calibration_H17_and_main_sequence_Scoville2017/fitting_dzliu_function/best_fit_via_MCMC.json'
    popt = [0.42803584038283427, -1.0941775869942276, -0.168209612312487, 1.207999567254241, 0.08377200755981207, 0.06373985587594078]
    a, b, c, d, ak, ck = popt
    if z is None and cosmic_age is None:
        print('Error! Please input either z or cosmic_age')
        sys.exit()
    else:
        if cosmic_age is not None:
            cosmoAge = cosmic_age
        else:
            cosmoAge = cosmo.age(z).value
    log10_Ratio_M_molgas_M_star = (a+ak*(lgMstar-10.0))*DeltaMS + b*(lgMstar-10.0) + (c+ck*(lgMstar-10.0))*(cosmoAge) + d
    Ratio_M_molgas_M_star = 10**log10_Ratio_M_molgas_M_star
    return Ratio_M_molgas_M_star

#def calc_gas_fraction_A3COSMOS(cosmoAge, lgMstar=10.5, DeltaMS=0.0):
#    # A3COSMOS 
#    # Sample 20181203 
#    # fitting '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plots/Plot_z_deltaGas_data_fitting/best_fit_via_CurveFit.json'
#    #popt = [0.7176359685433024, -0.3304845603259402, 0.008136603197619791, -0.07859177961327789, -0.026455117351337622] #<20190122_with_Leroy_GDR_with_Genzel_Eq12a_MZR_with_KMT09_fmol>#
#    #popt = [0.6829771406895456, -0.5516107864334655, 0.021790228645462928, -0.1558086567222849, 0.752820037618907]
#    #popt = [0.6569, -0.6628, 0.02804, -0.1689, +0.9282] #<20190205_with_RemyRuyer_GDR_with_Genzel_Eq12a_MZR_with_KMT09_fmol>#
#    #popt = [0.4132, -0.7011, 0.02891, -0.1637, 1.009] #<20190207># fixed xCOLDGASS Z_PP04_O3N2
#    #popt = [0.4496, -0.6859, 0.02930, -0.1602, 0.9698] #<20190208># still fitting Bertemes and Lee's data
#    popt = [0.4635, -0.7033, 0.03099, -0.1462, 0.8439] #<20190301># Eq850 this work, RemyRuyer GDR, Genzel2015_Eq12a_with_dzliu_limit MZR, KMT09 fmol, fitting all available data '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plots/Plot_z_deltaGas_data_fitting/20190301_with_gas_mass_calibration_Eq850_dzliu/20190301_fitting_dzliu_function/'
#    a, b, c, d, e = popt
#    log10_Ratio_M_molgas_M_star = a*DeltaMS + (b+c*cosmoAge)*(lgMstar-10.0) + d*(cosmoAge) + e
#    Ratio_M_molgas_M_star = 10**log10_Ratio_M_molgas_M_star
#    return Ratio_M_molgas_M_star

#def calc_gas_fraction_A3COSMOS_with_dzliu_850_gas_mass_calibration(cosmoAge, lgMstar=10.5, DeltaMS=0.0):
#    # A3COSMOS 
#    # Sample 20181203 
#    # fitting '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plots/Plot_z_deltaGas_data_fitting/20190122_with_Leroy_GDR_with_Genzel_Eq12a_MZR_with_dzliu_gas_mass_calibration_with_KMT09_fmol/best_fit_via_CurveFit.json'
#    #popt = [0.7176359685433024, -0.3304845603259402, 0.008136603197619791, -0.07859177961327789, -0.026455117351337622] #<20190122_with_Leroy_GDR_with_Genzel_Eq12a_MZR_with_KMT09_fmol>#
#    popt = [0.4635, -0.7033, 0.03099, -0.1462, 0.8439] #<20190301># Eq850 this work, RemyRuyer GDR, Genzel2015_Eq12a_with_dzliu_limit MZR, KMT09 fmol, fitting all available data '/Users/dzliu/Cloud/GitLab/AlmaCosmos/Plots/Plot_z_deltaGas_data_fitting/20190301_with_gas_mass_calibration_Eq850_dzliu/20190301_fitting_dzliu_function/'
#    a, b, c, d, e = popt
#    log10_Ratio_M_molgas_M_star = a*DeltaMS + (b+c*cosmoAge)*(lgMstar-10.0) + d*(cosmoAge) + e
#    Ratio_M_molgas_M_star = 10**log10_Ratio_M_molgas_M_star
#    return Ratio_M_molgas_M_star






# 
# def core functions (adapted from "a_dzliu_code_fit_z_tauDepl_v6.py")
# 
def func_deltaGas_dzliu_log(pars, a, b, c, d, e): 
    if type(pars) is dict:
        for t in ['tauDepl','DeltaMS','lgMstar','lgSFR','cosmoAge','z']:
            if t in pars:
                locals()[t] = copy.copy(pars[t])
    else:
        tauDepl,DeltaMS,lgMstar,lgSFR,cosmoAge,z = pars
    # 20180331a
    # deltaGas_model_log = a*DeltaMS + b*(lgMstar-10.0) + c*numpy.power(numpy.log10(1+z)) + d
    # 20180331b
    # deltaGas_model_log = a*DeltaMS + b*(lgMstar-10.0) + c*numpy.power((1+z),d) + e # NOT GOOD
    # 20180331c
    # deltaGas_model_log = a*DeltaMS + b*(lgMstar-10.0) + c*numpy.power(numpy.log10(1+z),1) + d*numpy.power(numpy.log10(1+z),2) + e
    # 20180331d
    # deltaGas_model_log = a*DeltaMS + b*(lgMstar-10.0) + c*numpy.power(numpy.log10(1+z),d) + e
    # 20180331e
    # deltaGas_model_log = a*DeltaMS + b*(lgMstar-10.0)
    # deltaGas_model = numpy.power(10, deltaGas_model_log) * e * numpy.power(z/c, d) * numpy.exp(-z/c)
    # 20180331f
    #deltaGas_model_log = a*DeltaMS + b*(lgMstar-10.0)
    #deltaGas_model = numpy.power(10, deltaGas_model_log) * e * (z - numpy.power(z/c,d))
    # 20180331g
    #cosmoAge = cosmo.age(z).value
    #deltaGas_model_log = a*DeltaMS*(1-numpy.exp(-DeltaMS/d)) + b*(lgMstar-10.0) + c*(cosmoAge) + e
    #deltaGas_model_log = a*DeltaMS + b*(lgMstar-10.0) + c*(cosmoAge) + d
    #deltaGas_model = numpy.power(10, deltaGas_model_log) # * numpy.power((1+z),c)
    # 20180401a
    #deltaGas_model_log = a*DeltaMS + b*(lgMstar-10.0) + c*numpy.log(cosmoAge) + d
    #deltaGas_model = numpy.power(10, deltaGas_model_log) # * numpy.power((1+z),c)
    # 20180401b
    #deltaGas_model_log = a*DeltaMS + b*(lgMstar-10.0) + c*(cosmoAge) + d
    # 20181030a
    #deltaGas_model_log = a*DeltaMS + b*(lgMstar-10.0) + c*(1-np.exp(cosmoAge)) + d
    # 20181102a
    #deltaGas_model_log = (np.log10(tauDepl))*0.0 + a*DeltaMS + b*(lgSFR-2.0)*0.0 + c*(lgMstar-10.0) + d*(cosmoAge) + e
    deltaGas_model_log = (np.log10(tauDepl))*0.0 + a*DeltaMS + (lgSFR-2.0)*0.0 + (b+c*cosmoAge)*(lgMstar-10.0) + d*(cosmoAge) + e
    # 
    return deltaGas_model_log

def func_deltaGas_dzliu(pars, a, b, c, d, e): 
    if type(pars) is dict:
        deltaGas_model_log = func_deltaGas_dzliu_log(pars, a, b, c, d, e)
    else:
        tauDepl,DeltaMS,lgMstar,lgSFR,cosmoAge,z = pars
        deltaGas_model_log = func_deltaGas_dzliu_log((tauDepl,DeltaMS,lgMstar,lgSFR,cosmoAge,z), a, b, c, d, e)
    deltaGas_model = numpy.power(10, deltaGas_model_log)
    return deltaGas_model

def func_deltaGas_Tacconi2018_log(pars, a, b, c, d): 
    DeltaMS,lgMstar,z = pars
    deltaGas_model_log = a*DeltaMS + b*(lgMstar-10.7) + c*(numpy.log10(1+z)-0.66)**2 + d
    return deltaGas_model_log

def func_deltaGas_Tacconi2018(pars):
    DeltaMS,lgMstar,z = pars
    deltaGas_model_log = func_deltaGas_Tacconi2018_log((DeltaMS,lgMstar,z), 0.53, -0.35, -3.62, 0.12)
    deltaGas_model = numpy.power(10, deltaGas_model_log) # * numpy.power((1+z),c)
    return deltaGas_model

def func_deltaGas_Scoville2017_log(pars, a, b, c, d): 
    DeltaMS,lgMstar,z = pars
    deltaGas_model_log = a*DeltaMS + b*(lgMstar-10.0) + c*(numpy.log10(1+z)) + d
    return deltaGas_model_log

def func_deltaGas_Scoville2017(pars):
    DeltaMS,lgMstar,z = pars
    deltaGas_model_log = func_deltaGas_Scoville2017_log((DeltaMS,lgMstar,z), 0.32, -0.70, 1.84, numpy.log10(0.71))
    deltaGas_model = numpy.power(10, deltaGas_model_log) # * numpy.power((1+z),c)
    return deltaGas_model








