#!/usr/bin/env python
# 
# Usage:
#    from calc_fmol import ( calc_fmol_from_metalZ_following_Krumholz2009, 
#                            calc_fmol_from_metalZ_following_Dave2016, 
#                            calc_fmol_from_metalZ_following_Popping2014 )
# 
# 

from __future__ import print_function

import os, sys, re, json, time, astropy
import numpy as np
from astropy.table import Table, Column, hstack
from copy import copy

from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.27, Tcmb0=2.725)

if sys.version_info.major >= 3:
    long = int
else:
    pass



# 
# def functions
# these functions are copied from 'AlmaCosmos/Pipeline/a3cosmos-Sample-Selection/a_dzliu_code_step_8_output_datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas.py'
# 
def calc_Lv_from_flux_density(obs_flux_mJy, obs_wavelength_um, redshift):
    lumdist_Mpc = cosmo.luminosity_distance(redshift).value # Mpc
    Lv = 4 * np.pi * lumdist_Mpc**2 * obs_flux_mJy / (1.+redshift) / 40.31970 / 1e9 # Lsun Hz-1 rest-frame, using 1 Lsun Mpc^-2 = 40.31970 mJy GHz.
    Lv = Lv * 3.839e33 # erg s-1 Hz-1
    return Lv

def calc_Lv_850um_from_RF_850um(RF_850um_mJy, redshift):
    lumdist_Mpc = cosmo.luminosity_distance(redshift).value # Mpc
    vLv_850um = 4 * np.pi * lumdist_Mpc**2 * RF_850um_mJy / 40.31970 * (2.99792458e5/(850.0)) # Lsun, using 1 Lsun Mpc^-2 = 40.31970 mJy GHz.
    Lv_850um = vLv_850um * 3.839e33 / ((2.99792458e5/850.0)*1e9) # erg s-1 Hz-1
    Lv_850um_2 = 1.19e27 * (RF_850um_mJy/1e3) / (1.+redshift) * lumdist_Mpc**2 # erg s-1 Hz-1, see Hughes 2017 - 1702.07350.pdf (page 3 right)
    print(Lv_850um, Lv_850um_2) # checked that `calc_Lv_850um_from_RF_850um(1., 2.0)` gives "1.0156146866563148e+32 1.0101020437871695e+32"
    return Lv_850um

def calc_M_mol_Hughes2017(S_850um_mJy, redshift):
    # Hughes+2017 
    # Here we take their molecular gas mass calibration. They also provided M_ISM calibration. 
    # For massive (log (vLv / [erg s-1 Hz-1]) ~ 32) galaxies, M_ISM - M_mol = 0.01 dex.
    # For low mass (log (vLv / [erg s-1 Hz-1]) ~ 30) galaxies, M_ISM - M_mol = 0.13 dex.
    # 
    Lv_850um_erg_s_Hz = calc_Lv_from_flux_density(S_850um_mJy, 850.0*(1.+redshift), redshift)
    M_mol_Msun = np.power(10.0, np.log10(Lv_850um_erg_s_Hz) * 0.93 - 17.74) # Hughes 2017 - 1702.07350.pdf (page 4, equation 5)
    return M_mol_Msun

def calc_M_total_gas_Hughes2017(S_850um_mJy, redshift):
    # Hughes+2017 
    # Here we take their molecular gas mass calibration. They also provided M_ISM calibration. 
    # For massive (log (vLv / [erg s-1 Hz-1]) ~ 32) galaxies, M_ISM - M_mol = 0.01 dex.
    # For low mass (log (vLv / [erg s-1 Hz-1]) ~ 30) galaxies, M_ISM - M_mol = 0.13 dex.
    # 
    Lv_850um_erg_s_Hz = calc_Lv_from_flux_density(S_850um_mJy, 850.0*(1.+redshift), redshift)
    M_total_gas_Msun = np.power(10.0, np.log10(Lv_850um_erg_s_Hz) * 0.86 - 15.38) # Hughes 2017 - 1702.07350.pdf (page 4, equation 4, for total ISM mass)
    return M_total_gas_Msun

def calc_M_mol_Groves2015(obs_flux_mJy, obs_wavelength_um, redshift):
    # Groves+2015
    # They gas mass calibration are for total atomic+mol gas, 
    # but they also mentioned that for massive galaxies M_mol/M_ISM ~ 1, see their last paragraph in Sect. 3.3. 
    # 
    #Brent_W = np.array([250., 350., 500.])
    #Brent_A = np.array([1.57, 1.49, 1.72]) # see Groves, bibcode=2015ApJ...799...96G, Table 6
    #Brent_B = np.array([0.86, 0.92, 0.96]) # see Groves, bibcode=2015ApJ...799...96G, Table 6, log(Mgas) = A + B * log(nuLnu_W)
    Brent_W = np.array([70.0, 100., 160., 250., 350., 500.])
    Brent_A = np.array([4.60, 3.27, 1.91, 1.17, 1.17, 1.44]) # see Groves, bibcode=2015ApJ...799...96G, Table 5, for log(M⋆/M⊙) > 9 sample
    Brent_B = np.array([0.50, 0.63, 0.78, 0.90, 0.90, 0.99]) # see Groves, bibcode=2015ApJ...799...96G, Table 5, for log(M⋆/M⊙) > 9 sample, log(Mgas) = A + B * log(nuLnu_W)
    rest_wavelength_um = obs_wavelength_um / (1.+redshift)
    dL = cosmo.luminosity_distance(redshift).value # Mpc
    Brent_log_M_mol_Msun = np.zeros((len(obs_flux_mJy),len(Brent_W),))
    #print(Brent_log_M_mol_Msun)
    for i in range(len(Brent_W)):
        A = Brent_A[i]
        B = Brent_B[i]
        vLv = (2.99792458e5/obs_wavelength_um)*(obs_flux_mJy)*(4*np.pi*dL**2)/40.31970 # "/40.31970" converts mJy GHz to Lsun Mpc-2 # 1.19e27 * Snu / (1+z) * dL**2 # erg s-1 Hz-1, see Hughes 2017 - 1702.07350.pdf (page 3 right)
        log_M_ISM_Msun = A + B * np.log10(vLv) # Groves et al. 2015, page 6, right middle-to-bottom
        Brent_log_M_mol_Msun[:,i] = log_M_ISM_Msun # calculate one M_mol for each Brent parameter, here we assume M_ISM is M_mol for massive galaxies.
    # 
    log_M_mol_Msun = []
    for i in range(len(obs_flux_mJy)):
        if rest_wavelength_um[i] >= np.min(Brent_W) and rest_wavelength_um[i] <= np.max(Brent_W):
            Brent_Interp = scipy.interpolate.interp1d(np.log10(Brent_W), Brent_log_M_mol_Msun[i,:])
            log_M_mol_Msun.append(Brent_Interp(np.log10(rest_wavelength_um[i])))
        else:
            log_M_mol_Msun.append(np.nan)
    M_mol_Msun = 10**np.array(log_M_mol_Msun)
    # 
    Schinnerer2016_A = 1.57 - 8e-4 * (rest_wavelength_um-250.0)
    Schinnerer2016_B = 0.86 + 6e-4 * (rest_wavelength_um-250.0)
    lgMH2_Schinnerer2016 = Schinnerer2016_A + Schinnerer2016_B * np.log10(vLv) # Schinnerer et al. 2016, equation (1)
    # 
    return M_mol_Msun

def calc_M_total_gas_Groves2015(obs_flux_mJy, obs_wavelength_um, redshift):
    # Groves+2015
    # They gas mass calibration are for total atomic+mol gas, 
    # but they also mentioned that for massive galaxies M_mol/M_ISM ~ 1, see their last paragraph in Sect. 3.3. 
    # 
    #Brent_W = np.array([250., 350., 500.])
    #Brent_A = np.array([1.57, 1.49, 1.72]) # see Groves, bibcode=2015ApJ...799...96G, Table 6
    #Brent_B = np.array([0.86, 0.92, 0.96]) # see Groves, bibcode=2015ApJ...799...96G, Table 6, log(Mgas) = A + B * log(nuLnu_W)
    Brent_W = np.array([70.0, 100., 160., 250., 350., 500.])
    Brent_A = np.array([4.15, 3.62, 3.52, 3.17, 3.08, 3.19]) # see Groves, bibcode=2015ApJ...799...96G, Table 5, for all log(M⋆/M⊙) sample
    Brent_B = np.array([0.55, 0.60, 0.61, 0.69, 0.74, 0.78]) # see Groves, bibcode=2015ApJ...799...96G, Table 5, for all log(M⋆/M⊙) sample, log(Mgas) = A + B * log(nuLnu_W)
    rest_wavelength_um = obs_wavelength_um / (1.+redshift)
    dL = cosmo.luminosity_distance(redshift).value # Mpc
    Brent_log_M_mol_Msun = np.zeros((len(obs_flux_mJy),len(Brent_W),))
    #print(Brent_log_M_mol_Msun)
    for i in range(len(Brent_W)):
        A = Brent_A[i]
        B = Brent_B[i]
        vLv = (2.99792458e5/obs_wavelength_um)*(obs_flux_mJy)*(4*np.pi*dL**2)/40.31970 # "/40.31970" converts mJy GHz to Lsun Mpc-2 # 1.19e27 * Snu / (1+z) * dL**2 # erg s-1 Hz-1, see Hughes 2017 - 1702.07350.pdf (page 3 right)
        log_M_ISM_Msun = A + B * np.log10(vLv) # Groves et al. 2015, page 6, right middle-to-bottom
        Brent_log_M_mol_Msun[:,i] = log_M_ISM_Msun # calculate one M_mol for each Brent parameter, here we assume M_ISM is M_mol for massive galaxies.
    # 
    log_M_mol_Msun = []
    for i in range(len(obs_flux_mJy)):
        if rest_wavelength_um[i] >= np.min(Brent_W) and rest_wavelength_um[i] <= np.max(Brent_W):
            Brent_Interp = scipy.interpolate.interp1d(np.log10(Brent_W), Brent_log_M_mol_Msun[i,:]) # interpolate the coefficients <TODO>
            log_M_mol_Msun.append(Brent_Interp(np.log10(rest_wavelength_um[i])))
        else:
            log_M_mol_Msun.append(np.nan)
    M_mol_Msun = 10**np.array(log_M_mol_Msun)
    # 
    return M_mol_Msun

def calc_M_mol_Scoville2017_RF(S_850um_mJy, redshift):
    # As mentioned in Hughes+2017, Scoville+2014 calibrated the atomic+mol gas mass, 
    #   but Scoville+2016 and 2017 calibrated only the mol gas mass.
    # 
    Lv_850um_erg_s_Hz = calc_Lv_from_flux_density(S_850um_mJy, 850.0*(1.+redshift), redshift)
    M_mol_Msun = Lv_850um_erg_s_Hz / 6.7e19 # Scoville et al. 2016ApJ...820...83S (appendix and also page 16 equation 16) -- molecular gas including Helium
    #lumdist_Mpc = cosmo.luminosity_distance(redshift).value # Mpc
    #M_mol_Msun_2 = 1.78 * RF_850um_mJy/(1.+redshift) * (lumdist_Mpc/1e3)**2 * 1e10 # Scoville et al. 2016ApJ...820...83S (appendix and also page 16 equation 16) -- molecular gas including Helium
    #print('%0.6e %0.6e' % (M_mol_Msun, M_mol_Msun_2) ) # checked that `calc_M_mol_Scoville2017(1., 2.0)` "1.515843e+12 1.510909e+12"
    return M_mol_Msun

def calc_M_mol_Scoville2017(obs_flux_mJy, obs_wavelength_um, redshift):
    # As mentioned in Hughes+2017, Scoville+2014 calibrated the atomic+mol gas mass, 
    #   but Scoville+2016 and 2017 calibrated only the mol gas mass.
    # 
    Lv_obs_erg_s_Hz = calc_Lv_from_flux_density(obs_flux_mJy, obs_wavelength_um, redshift)
    beta = 1.8 # do K-correction assuming planck function, beta=1.8, Tdust=25K
    Tdust = 25.0 # K
    cal_Lv_obs_erg_s_Hz = blackbody_nu(obs_wavelength_um/(1.+redshift) * u.um, Tdust * u.K) * np.power(2.99792458e5/(obs_wavelength_um/(1.+redshift)), beta)
    cal_Lv_850um_erg_s_Hz = blackbody_nu(850.0 * u.um, 25.0 * u.K) * np.power(2.99792458e5/850.0, beta)
    #print(cal_Lv_obs_erg_s_Hz)
    #print(cal_Lv_850um_erg_s_Hz)
    #print(cal_Lv_obs_erg_s_Hz / cal_Lv_850um_erg_s_Hz)
    Lv_850um_erg_s_Hz = Lv_obs_erg_s_Hz * (cal_Lv_850um_erg_s_Hz/cal_Lv_obs_erg_s_Hz)
    M_mol_Msun = Lv_850um_erg_s_Hz / 6.7e19 # Scoville et al. 2016ApJ...820...83S (appendix and also page 16 equation 16) -- molecular gas including Helium
    return M_mol_Msun



# 
# main function 1
# 
def calc_gas_mass_from_dust_continuum(obs_wavelength_um, obs_flux_mJy, SED_flux_at_obs_wavelength = None, SED_flux_at_rest_850um = None, z = None, method = 'Hughes2017'):
    # 
    # check user input
    # 
    if obs_wavelength_um is None:
        raise ValueError('Error! Please input \'obs_wavelength_um\' for calc_gas_mass_from_dust_continuum_with_band_conversion()!')
    # 
    if obs_flux_mJy is None:
        raise ValueError('Error! Please input \'obs_flux_mJy\' for calc_gas_mass_from_dust_continuum_with_band_conversion()!')
    # 
    if z is None:
        raise ValueError('Error! Please input \'z\' for calc_gas_mass_from_dust_continuum_with_band_conversion()!')
    # 
    # convert observed dust continuum at observed wavelength to (cold molecular) gas mass using the input method. 
    # 
    allowed_methods = ['Hughes2017', 'Groves2015', 'Scoville2017']
    # 
    if method == 'Hughes2017':
        # 
        if SED_flux_at_obs_wavelength is None or SED_flux_at_rest_850um is None:
            raise ValueError('Error! Please input \'SED_flux_at_obs_wavelength\' and \'SED_flux_at_rest_850um\' for calc_gas_mass_from_dust_continuum_with_band_conversion() with method \"Hughes2017\"!')
        # 
        M_mol_gas = calc_M_mol_Hughes2017(SED_flux_at_rest_850um / SED_flux_at_obs_wavelength * obs_flux_mJy, z) # band conversion used MAGPHYS SED, then alphamol850 following Hughes+2017.
        # 
    elif method == 'Groves2015':
        # 
        M_mol_gas = calc_M_mol_Groves2015(obs_flux_mJy, obs_wavelength_um, z) # no need band conversion here, it is inside the called function.
        # 
    elif method == 'Scoville2017':
        # 
        M_mol_gas = calc_M_mol_Scoville2017(obs_flux_mJy, obs_wavelength_um, z) # no need band conversion here, it is inside the called function by assuming a simple dust RJ tail slope exactly following Scoville+2017.
        # 
    else:
        # 
        raise NotImplementedError('Sorry! The gas mass calibration has not been implemented for the method "%s"! We have following methods: %s'%(method, allowed_methods))
        # 
    # 
    return M_mol_gas, method




# 
# def functions
# 
def calc_molecular_hydrogen_fraction(metalZOH, U_MW = None, method = 'Krumholz2009'):
    # 
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if not script_dir in sys.path:
        sys.path.append(script_dir)
    import inspect
    # 
    allowed_methods = ['Krumholz2009', 'Dave2016', 'Popping2014']
    # 
    if re.match(r'.*\bKrumholz2009\b.*', method, re.IGNORECASE):
        # 
        if not ('calc_fmol_from_metalZ_following_Krumholz2009' in dir()):
            from calc_fmol import calc_fmol_from_metalZ_following_Krumholz2009
        # 
        f_mol = calc_fmol_from_metalZ_following_Krumholz2009(metalZOH)
        # 
    elif re.match(r'.*\bDave2016\b.*', method, re.IGNORECASE):
        # 
        if not ('calc_fmol_from_metalZ_following_Dave2016' in dir()):
            from calc_fmol import calc_fmol_from_metalZ_following_Dave2016
        # 
        f_mol = calc_fmol_from_metalZ_following_Dave2016(metalZOH)
        # 
    elif re.match(r'.*\bPopping2014\b.*', method, re.IGNORECASE):
        # 
        if not ('calc_fmol_from_metalZ_following_Popping2014' in dir()):
            from calc_fmol import calc_fmol_from_metalZ_following_Popping2014
        # 
        f_mol = calc_fmol_from_metalZ_following_Popping2014(metalZOH, U_MW = U_MW)
        # 
    else:
        # 
        raise NotImplementedError('Sorry! The molecular hydrogen fraction calculation has not been implemented for the method "%s"! We have following methods: %s'%(method, allowed_methods))
        # 
    return f_mol
    


# 
# main function 2
# 
def calc_gas_mass_from_dust_mass(M_dust, M_star = None, SFR = None, metallicity = None, GDR = None, z = None, method = ''):
    # 
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if not script_dir in sys.path:
        sys.path.append(script_dir)
    #<TETS>#import inspect
    #<TETS>#print('inspect.isfunction(calc_deltaGD_from_metalZ_following_Genzel2015)', inspect.isfunction('calc_deltaGD_from_metalZ_following_Genzel2015'))
    #<TETS>#print('in dir()', 'calc_deltaGD_from_metalZ_following_Genzel2015' in dir())
    #<TETS>#from calc_delta_GD import calc_deltaGD_from_metalZ_following_Genzel2015
    #<TETS>#print('inspect.isfunction(calc_deltaGD_from_metalZ_following_Genzel2015)', inspect.isfunction('calc_deltaGD_from_metalZ_following_Genzel2015'))
    #<TETS>#print('in dir()', 'calc_deltaGD_from_metalZ_following_Genzel2015' in dir())
    # 
    M_total_gas = None
    M_mol_gas = None
    # 
    allowed_methods = ['Magdis2012', 'Genzel2015a', 'Genzel2015b', 'GDR']
    # 
    # if method is not given by the user
    if method == '' and GDR is not None:
        method = 'GDR'
    elif method == '':
        method = 'default'
    # 
    if re.match(r'.*\bMagdis2012\b.*', method, re.IGNORECASE):
        # 
        if metallicity is None:
            if M_star is None or SFR is None:
                raise ValueError('Error! Please input \'M_star\' and \'SFR\' or \'metallicity\' for calc_gas_mass_from_dust_mass() with the method \"Magdis2012\"!')
            else:
                if not ('calc_metalZ_from_FMR_following_Mannucci2010' in dir()):
                    from calc_metal_Z import calc_metalZ_from_FMR_following_Mannucci2010
                metallicity = calc_metalZ_from_FMR_following_Mannucci2010(M_star, SFR)
                method += ', lgMstar=%.2f, SFR=%g, [Mannucci2010] metallicity=%.2f'%(np.log10(M_star), SFR, metallicity)
        else:
            method += ', metallicity=%.2f'%(metallicity)
        # 
        if GDR is None:
            if not ('calc_deltaGD_from_metalZ_following_Magdis2012' in dir()):
                from calc_delta_GD import calc_deltaGD_from_metalZ_following_Magdis2012
            GDR = calc_deltaGD_from_metalZ_following_Magdis2012(metallicity)
            method += ', [Magdis2012] GDR=%.2f'%(GDR)
        else:
            # 
            method += ', GDR=%.2f'%(GDR)
        # 
        M_mol_gas = M_dust * GDR
        # 
        # if method is like 'Magdis2012+Krumholz2009'
        if method.find('+') >= 0:
            f_mol = calc_molecular_hydrogen_fraction(metallicity, method = method)
            M_total_gas = M_mol_gas
            M_mol_gas = M_total_gas * f_mol
        # 
    elif re.match(r'.*\bGenzel2015a\b.*', method, re.IGNORECASE):
        # 
        if metallicity is None:
            if M_star is None or SFR is None or z is None:
                raise ValueError('Error! Please input (\'M_star\' and \'SFR\' and \'z\') or \'metallicity\' for calc_gas_mass_from_dust_mass() with the method \"Genzel2015a\"!')
            else:
                if not ('calc_metalZ_from_FMR_following_Genzel2015a' in dir()):
                    from calc_metal_Z import calc_metalZ_from_FMR_following_Genzel2015a
                metallicity = calc_metalZ_from_FMR_following_Genzel2015a(M_star, SFR, z)
                method += ', lgMstar=%.2f, SFR=%g, z=%g, [Genzel2015a] metallicity=%.2f'%(np.log10(M_star), SFR, z, metallicity)
        else:
            method += ', metallicity=%.2f'%(metallicity)
        # 
        if GDR is None:
            if not ('calc_deltaGD_from_metalZ_following_Genzel2015' in dir()):
                from calc_delta_GD import calc_deltaGD_from_metalZ_following_Genzel2015
            GDR = calc_deltaGD_from_metalZ_following_Genzel2015(metallicity)
            method += ', [Genzel2015] GDR=%.2f'%(GDR)
        else:
            # 
            method += ', GDR=%.2f'%(GDR)
        # 
        M_mol_gas = M_dust * GDR
        # 
        # if method is like 'Genzel2015a+Krumholz2009'
        if method.find('+') >= 0:
            f_mol = calc_molecular_hydrogen_fraction(metallicity, method = method)
            M_total_gas = M_mol_gas
            M_mol_gas = M_total_gas * f_mol
        # 
    elif re.match(r'.*\b(Genzel2015b|Genzel2015)\b.*', method, re.IGNORECASE):
        # 
        if metallicity is None:
            if M_star is None or SFR is None:
                raise ValueError('Error! Please input (\'M_star\' and \'SFR\') or \'metallicity\' for calc_gas_mass_from_dust_mass() with the method \"Genzel2015b\"!')
            else:
                if not ('calc_metalZ_from_FMR_following_Genzel2015b' in dir()):
                    from calc_metal_Z import calc_metalZ_from_FMR_following_Genzel2015b
                metallicity = calc_metalZ_from_FMR_following_Genzel2015b(M_star, SFR)
                method += ', lgMstar=%.2f, SFR=%g, [Genzel2015b] metallicity=%.2f'%(np.log10(M_star), SFR, metallicity)
        # 
        if GDR is None:
            if not ('calc_deltaGD_from_metalZ_following_Genzel2015' in dir()):
                from calc_delta_GD import calc_deltaGD_from_metalZ_following_Genzel2015
            GDR = calc_deltaGD_from_metalZ_following_Genzel2015(metallicity)
            method += ', [Genzel2015] GDR=%.2f'%(GDR)
        else:
            # 
            method += ', GDR=%.2f'%(GDR)
        # 
        M_mol_gas = M_dust * GDR
        # 
        # if method is like 'Genzel2015a+Krumholz2009'
        if method.find('+') >= 0:
            f_mol = calc_molecular_hydrogen_fraction(metallicity, method = method)
            M_total_gas = M_mol_gas
            M_mol_gas = M_total_gas * f_mol
        # 
        # 
    elif re.match(r'.*\b(default)\b.*', method, re.IGNORECASE):
        # 
        if metallicity is None:
            if M_star is None or SFR is None or z is None:
                raise ValueError('Error! Please input (\'M_star\' and \'SFR\' and \'z\') or \'metallicity\' for calc_gas_mass_from_dust_mass() with the method \"default\"!')
            else:
                if not ('calc_metalZ_from_FMR_following_Genzel2015ab_combined_by_dzliu' in dir()):
                    from calc_metal_Z import calc_metalZ_from_FMR_following_Genzel2015ab_combined_by_dzliu
                metallicity = calc_metalZ_from_FMR_following_Genzel2015ab_combined_by_dzliu(M_star, SFR, z)
                method += ', lgMstar=%.2f, SFR=%g, [Genzel2015abcomb] metallicity=%.2f'%(np.log10(M_star), SFR, metallicity)
        else:
            method += ', metallicity=%.2f'%(metallicity)
        # 
        if GDR is None:
            if not ('calc_deltaGD_from_metalZ_following_RemyRuyer2014b' in dir()):
                from calc_delta_GD import calc_deltaGD_from_metalZ_following_RemyRuyer2014b
            GDR = calc_deltaGD_from_metalZ_following_RemyRuyer2014b(metallicity)
            method += ', [RemyRuyer2014b] GDR=%.2f'%(GDR)
        else:
            # 
            method += ', GDR=%.2f'%(GDR)
        # 
        M_mol_gas = M_dust * GDR
        # 
        # if method is like 'Genzel2015a+Krumholz2009'
        if method.find('+') >= 0:
            f_mol = calc_molecular_hydrogen_fraction(metallicity, method = method)
            M_total_gas = M_mol_gas
            M_mol_gas = M_total_gas * f_mol
        # 
    elif re.match(r'.*\bGDR\b.*', method, re.IGNORECASE):
        # Assuming a single GDR
        # 
        if GDR is None:
            raise ValueError('Error! Please input \'GDR\' for calc_gas_mass_from_dust_mass() with the method \"GDR\"!')
            #if metallicity is None:
            #    raise ValueError('Error! Please input \'metallicity\' or \'GDR\' for calc_gas_mass_from_dust_mass() with the method \"GDR\"!')
            #else:
            #    if not ('calc_deltaGD_from_metalZ_following_Genzel2015' in dir()):
            #        from calc_delta_GD import calc_deltaGD_from_metalZ_following_Genzel2015
            #    GDR = calc_deltaGD_from_metalZ_following_Genzel2015(metallicity)
        else:
            # 
            method += ', GDR=%.2f'%(GDR)
        # 
        M_mol_gas = M_dust * GDR
        # 
        # if method is like 'Genzel2015a+Krumholz2009'
        if method.find('+') >= 0:
            f_mol = calc_molecular_hydrogen_fraction(metallicity, method = method)
            M_total_gas = M_mol_gas
            M_mol_gas = M_total_gas * f_mol
        # 
    else:
        # 
        raise NotImplementedError('Sorry! The gas mass calibration has not been implemented for the method "%s"! We have following methods: %s'%(method, allowed_methods))
        # 
    # 
    if M_total_gas is None:
        return M_mol_gas, method
    else:
        return M_mol_gas, M_total_gas, method






