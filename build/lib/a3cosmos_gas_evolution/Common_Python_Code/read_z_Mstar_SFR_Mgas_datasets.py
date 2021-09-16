#!/usr/bin/env python
# 
# changes
#     20190208 Cicone+ still use their CO(2-1)/CO(1-0)=1.0 assumption. 
# 

from __future__ import print_function

import os, sys, re, copy, json, time, datetime, shutil, astropy
import numpy as np
from astropy.table import Table, Column, MaskedColumn, hstack
import pandas as pd

if not (os.path.dirname(os.path.abspath(__file__)) in sys.path): sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import apply_cosmology
cosmo = apply_cosmology.cosmo

if sys.version_info.major >= 3:
    long = int
else:
    pass

sys.path.append(os.path.dirname(__file__))

from calc_alpha_CO import ( calc_alphaCO_from_metalZ_following_Wilson1995, 
                            calc_alphaCO_from_metalZ_following_Genzel2015a, 
                            calc_alphaCO_from_metalZ_following_Genzel2015b, 
                            calc_alphaCO_from_metalZ_following_Bolatto2013, 
                            calc_alphaCO_from_metalZ_following_Accurso2017, 
                            calc_alphaCO_from_metalZ_following_Bertemes2018, 
                            calc_alphaCO_from_metalZ_following_Tacconi2018, 
                          )
from calc_delta_GD import ( calc_deltaGD_from_metalZ_following_Leroy2011, 
                            calc_deltaGD_from_metalZ_following_Magdis2012, 
                            calc_deltaGD_from_metalZ_following_RemyRuyer2014a, 
                            calc_deltaGD_from_metalZ_following_RemyRuyer2014b, 
                          )
from calc_fmol import     ( calc_fmol_from_metalZ_following_Krumholz2009, 
                            calc_fmol_from_metalZ_following_Dave2016, 
                            calc_fmol_from_metalZ_following_Popping2014, 
                          )
from calc_metal_Z import  ( calc_metalZ_from_FMR_following_Genzel2015_Eq12a, 
                            calc_metalZ_from_FMR_following_Mannucci2010_Eq4, 
                            convert_metalZ_M08_to_metalZ_PP04_N2_polynomial, 
                            convert_metalZ_KK04_to_metalZ_PP04, 
                            calc_metalZ_from_FMR_with_dzliu_selection, 
                            calc_metalZ_from_FMR_following_Kewley2008_PP04_O3N2,
                          )






# 
# def 
# 
def calc_Sargent2014_sSFR(z, lgMstar=10.5, DeltaMS=0.0):
    return 0.095 * 10**(-0.21*(lgMstar-numpy.log10(5e10))) * numpy.exp(2.05*z/(1.0+0.16*z**1.54)) * 10**(DeltaMS)

def calc_Speagle2014_sSFR(cosmoAge, lgMstar=10.5, DeltaMS=0.0):
    return 10**((0.84 - 0.026*cosmoAge) * lgMstar - (6.51 - 0.11*cosmoAge)) / 10**(lgMstar) * 1e9 * 10**(DeltaMS)

def calc_Scoville2017_sSFR(z, lgMstar=10.5, DeltaMS=0.0):
    lgMstar_ref = 10.5
    SFR_MS_ref = 10**(0.59*lgMstar_ref-5.77)* np.power((1.0+z),(0.22*lgMstar_ref+0.59))
    SFR_MS = SFR_MS_ref * 10**(1.72-np.log10(1+np.power(10**(lgMstar-10.31),-1.07))) / 10**(1.72-np.log10(1+np.power(10**(lgMstar_ref-10.31),-1.07)))
    sSFR_MS = SFR_MS / 10**(lgMstar) * 1e9 # Gyr
    return sSFR_MS

def calc_sSFR_MS(lgMstar, z, cosmoAge=None):
    if cosmoAge is None:
        cosmoAge = cosmo.age(z).value
    sSFR_MS = calc_Speagle2014_sSFR(cosmoAge, lgMstar)
    # 
    #sSFR_MS = calc_Sargent2014_sSFR(z, lgMstar)
    # 
    return sSFR_MS



# 
# mask_data
# 
def mask_dataset(input_data, mask_CPA = True, mask_SED = True, mask_IMG = True, mask_known_zspec = False):
    # the input_data should be a dict with 'ID', ''
    # 
    # mask valid sources
    # 
    mask_valid_sources = (input_data['z']>0)
    
    if mask_SED:
        if os.path.isfile('datatable_discarded_sources_by_SED.txt'):
            list_SED = Table.read('datatable_discarded_sources_by_SED.txt', format='ascii.commented_header') # must use 'ascii.commented_header' otherwise got a bug
            mask_SED = np.isin(input_data['ID'], list_SED.columns[0].data) # SED ALMA band chisq identified spurious sources
            mask_valid_sources = np.logical_and( mask_valid_sources, ~mask_SED )
    
    if mask_CPA:
        if os.path.isfile('datatable_discarded_sources_by_CPA.txt'):
            list_CPA = Table.read('datatable_discarded_sources_by_CPA.txt', format='ascii.commented_header') # must use 'ascii.commented_header' otherwise got a bug
            mask_CPA = np.isin(input_data['ID'], list_CPA.columns[0].data) # counterpart association identified spurious sources
            mask_valid_sources = np.logical_and( mask_valid_sources, ~mask_CPA )
    
    if mask_IMG:
        if os.path.isfile('datatable_discarded_sources_by_IMG.txt'):
            list_IMG = Table.read('datatable_discarded_sources_by_IMG.txt', format='ascii.commented_header') # must use 'ascii.commented_header' otherwise got a bug
            mask_IMG = np.isin(input_data['ID'], list_IMG.columns[0].data) # bad image
            mask_valid_sources = np.logical_and( mask_valid_sources, ~mask_IMG )
    
    if mask_known_zspec:
        if os.path.isfile('datatable_known_zspec.txt'):
            list_known_zspec = Table.read('datatable_known_zspec.txt', format='ascii.commented_header') # must use 'ascii.commented_header' otherwise got a bug
            mask_known_zspec = np.isin(input_data['ID'], list_known_zspec.columns[0]) # known_zspec
            mask_valid_sources = np.logical_and( mask_valid_sources, mask_known_zspec )
    
    print('selecting %d data after masking' % (np.sum(mask_valid_sources)))
    output_data = copy.copy(input_data)
    for keyname in output_data:
        if not np.isscalar(input_data[keyname]):
            output_data[keyname] = np.array(input_data[keyname])[mask_valid_sources]
    
    return output_data



# 
# def
# 
def calc_metal_Z_high_z_method(M_star, SFR, z):
    #return calc_metalZ_from_FMR_with_dzliu_selection(M_star, SFR, z)
    return calc_metalZ_from_FMR_following_Genzel2015_Eq12a(M_star, z)
    
def calc_metal_Z_local_galaxy_method(M_star, SFR, z):
    #return calc_metalZ_from_FMR_with_dzliu_selection(M_star, SFR, z)
    return convert_metalZ_M08_to_metalZ_PP04_N2_polynomial(calc_metalZ_from_FMR_following_Mannucci2010_Eq4(M_star, SFR))
    #return calc_metalZ_from_FMR_following_Kewley2008_PP04_O3N2(M_star)




# 
# read data
# 
def read_datasets():
    # 
    # read tables
    # 
    datasets = []
    
    #which_gas_mass_calibration = 'dzliu 850'
    #which_gas_mass_calibration = 'H17' # Hughes+2017
    
    if True:
        ds = {}
        ds['label'] = 'This work (A3COSMOS)'
        ds['color'] = 'gold'
        ds['facecolor'] = 'gold'
        ds['edgecolor'] = 'k'
        ds['edgelinewidth'] = 0.5
        ds['alpha'] = 1.0
        ds['marker'] = 'o'
        ds['markersize'] = 15
        #ds['datatable'] = '/Users/dzliu/Work/AlmaCosmos/Samples/20180720/datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas.fits'
        #ds['datatable'] = '/Users/dzliu/Work/AlmaCosmos/Samples/20181203/datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas.fits'
        #ds['datatable'] = '/Users/dzliu/Work/AlmaCosmos/Samples/20181203/datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas_with_RemyRuyer2014_GDR_with_Genzel2015_Eq12a_MZR/datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas.fits'
        #ds['datatable'] = '/Users/dzliu/Work/AlmaCosmos/Samples/20181203/datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas.fits' # with_RemyRuyer2014_GDR_with_Genzel2015_Eq12a_MZR
        #print('Reading \'/Users/dzliu/Work/AlmaCosmos/Samples/20181203/datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas.fits\', updated on 2019-03-01 with RemyRuyer2014 GDR, KMT09 fmol, Genzel2015Eq12a MZR(z). ')
        #ds['datatable'] = '/Users/dzliu/Work/AlmaCosmos/Samples/20181203/datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas.fits' # updated on 2019-03-01 with RemyRuyer2014 GDR, KMT09 fmol, Genzel2015Eq12a MZR(z). 
        #ds['datatable'] = '/Users/dzliu/Work/AlmaCosmos/Samples/20181203/datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas_with_specz.fits' # 20190307 selected spec-z subsample
        tbfile = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_a3cosmos'+os.sep+'dataset_v20180801'+os.sep+'datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas.fits'
        tbinfo = open(os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_a3cosmos'+os.sep+'dataset_v20180801'+os.sep+'datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas.info.txt').readline().rstrip()
        print('Reading \'%s\' %s'%(tbfile, tbinfo))
        ds['datatable'] = tbfile
        tb = Table.read(ds['datatable'])
        mask = (tb['ID']!=850535) #<20190715># SED fitting Qz=99 not get prioritized bug
        tb = tb[mask]
        ds['ID'] = tb['ID']
        ds['z'] = tb['z']
        ds['SFR'] = tb['SFR']
        ds['Mstar'] = tb['Mstar']
        #ds['Mmol'] = tb['M_mol_gas_Method4_500_this_work'] # involves f_mol
        #ds['Mmol'] = tb['M_mol_gas_Method4_850_this_work'] # involves f_mol
        #ds['Mmol'] = tb['M_mol_gas_Method4_850_this_work'] #<20190122_with_Leroy_GDR_with_Genzel_Eq12a_MZR_with_KMT09_fmol>#
        #ds['Mmol'] = tb['M_mol_gas_Method4_850_this_work'] #<20190301_with_RemyRuyer_GDR_with_Genzel_Eq12a_with_dzliu_limit_MZR_with_KMT09_fmol>#
        #ds['Mmol'] = tb['M_mol_gas_Method2_850_Scoville2017']
        ds['Mmol'] = tb['M_mol_gas_Method2_850_Hughes2017'] # FINAL GOOD TO USE
        #ds['Mmol'] = tb['M_total_gas_Method1_GDR']
        #ds['Mmol'] = tb['M_mol_gas_Method3_500_Groves2015']
        #ds['Mmol'] = tb['deltaGas4'][mask] * tb['Mstar']
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        ds['deltaGasErr'] = np.sqrt( ((1.0/tb['SNRObs']))**2 + (0.2)**2 ) * ds['deltaGas'] # dust continuum SNRObs, gas conversion 0 dex scatter (Hughes2017), stellar mass 0.2 dex uncertainty
        ds['tauDeplErr'] = np.sqrt( ((1.0/tb['SNRObs']))**2 + (0.1)**2 ) * ds['tauDepl'] # dust continuum SNRObs, gas conversion 0 dex scatter (Hughes2017), SFR 0.1 dex uncertainty
        
        #if False:
        #    mask = (tb['SNRObs']>=2.5)
        #    for tcolname in ds:
        #        ds[tcolname] = ds[tcolname][mask]
        #    ds = mask_dataset(ds) # for A3COSMOS dataset, we need to mask more sources
        
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Tacconi+2018' # only extracted PHBISS1 & 2
        ds['color'] = 'seagreen'
        ds['facecolor'] = 'seagreen'
        ds['edgecolor'] = 'none'
        ds['alpha'] = 0.7
        ds['marker'] = 'D'
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_Tacconi2018_PHIBSS2/datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas_with_Survey_Number_GE_1.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = (tb['Mstar']>0)
        tb = tb[mask]
        ds['ID'] = tb['ID']
        ds['z'] = tb['z']
        ds['SFR'] = tb['SFR']
        ds['Mstar'] = tb['Mstar'] # 
        #ds['Mmol'] = tb['deltaGas'] * tb['Mstar'] # they used metallicity-dependent alphaCO. For metallicity they used Genzel2015_Eq12a. 
        ds['MetalZ'] = calc_metal_Z_high_z_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['Mmol'] = (tb['deltaGas'] * tb['Mstar']) / calc_alphaCO_from_metalZ_following_Tacconi2018(calc_metalZ_from_FMR_following_Genzel2015_Eq12a(ds['Mstar'], ds['z'])) * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        # using calc_alphaCO_from_metalZ_following_Genzel2015b() will have much higher alphaCO than using calc_alphaCO_from_metalZ_following_Genzel2015b(). 
        # comparing to their computed gas mass, my gas mass is slightly higher
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        ds['deltaGasErr'] = np.sqrt( (0.3)**2 + (0.2)**2 ) * ds['deltaGas'] # CO 0.3 dex uncertainty, gas conversion 0 dex uncertainty (alphaCO), stellar mass 0.2 dex uncertainty
        ds['tauDeplErr'] = np.sqrt( (0.3)**2 + (0.1)**2 ) * ds['tauDepl'] # CO 0.3 dex uncertainty, gas conversion 0 dex uncertainty (alphaCO), SFR 0.1 dex uncertainty
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Saintonge+2017'
        ds['color'] = 'blue'
        ds['facecolor'] = 'blue'
        ds['edgecolor'] = 'none'
        ds['alpha'] = 0.7
        ds['marker'] = '+'
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_Saintonge2017_xCOLDGASS/xCOLDGASS_PubCat.fits'
        tb = Table.read(ds['datatable'], format='fits')
        #mask = np.logical_and(np.logical_and(np.logical_and(np.logical_and(tb['SN_CO']>=3, tb['LOGMSTAR']>0), tb['LCO_COR']>0), tb['LOGSFR_BEST']>-99), tb['Z_PP04_O3N2']>0)
        mask = np.logical_and(np.logical_and(np.logical_and(tb['SN_CO']>=3, tb['LOGMSTAR']>0), tb['LCO_COR']>0), tb['LOGSFR_BEST']>-99)
        tb = tb[mask]
        ds['ID'] = tb['ID']
        ds['z'] = tb['Z_SDSS']
        ds['SFR'] = np.power(10.0, tb['LOGSFR_BEST'])
        ds['Mstar'] = np.power(10.0, tb['LOGMSTAR'])
        #ds['Mmol'] = tb['LCO_COR'] * tb['XCO_A17']
        ds['Mmol_Saintonge2017'] = tb['LCO_COR'] * tb['XCO_A17']
        ds['alphaCO_Saintonge2017'] = tb['XCO_A17']
        #ds['LPrmCO10'] = tb['LCO_COR']
        ds['LPrmCO10'] = 23.5 * tb['ICO_COR'] * np.pi/(4*np.log(2)) * (22.0)**2 * (tb['LUMDIST'])**2 * np.power((1.+tb['Z_SDSS']),-3)
        ds['LPrmCO10_err'] = 23.5 * tb['ICO_COR_ERR'] * np.pi/(4*np.log(2)) * (22.0)**2 * (tb['LUMDIST'])**2 * np.power((1.+tb['Z_SDSS']),-3)
        mask2 = (tb['Z_PP04_O3N2']>0) # we use their metalZ if valid, otherwise compute from MZR or FMR with calc_metal_Z_local_galaxy_method()
        ds['MetalZ'] = calc_metal_Z_local_galaxy_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['MetalZ'][mask2] = tb['Z_PP04_O3N2'][mask2] # their metalZ are derived from optical nebular lines (SDSS) with PP04 O3N2 calibration 
        ds['MetalZ_Mannucci2010_Eq4_Method'] = calc_metalZ_from_FMR_following_Mannucci2010_Eq4(ds['Mstar'], ds['SFR'])
        ds['MetalZ_Kewley2008_Method'] = calc_metalZ_from_FMR_following_Kewley2008_PP04_O3N2(ds['Mstar'])
        ds['alphaCO'] = calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['Mmol'] = (ds['LPrmCO10']) * ds['alphaCO']
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        ds['deltaGasErr'] = np.sqrt( (tb['ICO_COR_ERR']/tb['ICO_COR'])**2 + (0.2)**2 ) * ds['deltaGas'] # CO obs uncertainty, gas conversion 0 dex scatter (alphaCO), stellar mass 0.2 dex uncertainty
        ds['tauDeplErr'] = np.sqrt( (tb['ICO_COR_ERR']/tb['ICO_COR'])**2 + (tb['LOGSFR_ERR'])**2 ) * ds['tauDepl'] # CO obs uncertainty, gas conversion 0 dex scatter (alphaCO), SFR obs uncertainty
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
        # 
        #dsout = {}
        #for key in ['z','SFR','Mstar','Mmol','alphaCO','Mmol_Saintonge2017','alphaCO_Saintonge2017','MetalZ','MetalZ_Kewley2008_Method']:
        #    dsout[key] = copy.copy(ds[key])
        #for key in ['Z_PP04_N2','Z_PP04_O3N2','Z_MZR']:
        #    dsout[key] = copy.copy(tb[key][mask].data)
        #tbout = Table(dsout)
        #tbout.write('check_Saintonge2017_sample.fits', format='fits', overwrite=True)
    
    if True:
        ds = {}
        ds['label'] = 'Saintonge+2017 uplims'
        ds['color'] = 'blue'
        ds['facecolor'] = 'blue'
        ds['edgecolor'] = 'none'
        ds['alpha'] = 0.7
        ds['marker'] = 'uplims'
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_Saintonge2017_xCOLDGASS/xCOLDGASS_PubCat.fits'
        tb = Table.read(ds['datatable'], format='fits')
        #mask = np.logical_and(np.logical_and(np.logical_and(np.logical_and(tb['SN_CO']<3, tb['LOGMSTAR']>0), tb['RMS_CO']>0), tb['LOGSFR_BEST']>-99), tb['Z_PP04_O3N2']>0)
        #mask = np.logical_and(np.logical_and(np.logical_and(tb['SN_CO']>=3, tb['LOGMSTAR']>0), tb['LCO_COR']>0), tb['LOGSFR_BEST']>-99)
        mask = np.logical_and(np.logical_and(np.logical_and(tb['SN_CO']<3, tb['LOGMSTAR']>0), tb['LCO_COR']>0), tb['LOGSFR_BEST']>-99)
        tb = tb[mask]
        ds['ID'] = tb['ID']
        ds['z'] = tb['Z_SDSS']
        ds['SFR'] = np.power(10.0, tb['LOGSFR_BEST'])
        ds['Mstar'] = np.power(10.0, tb['LOGMSTAR'])
        #ds['Mmol'] = tb['LCO_COR'] * tb['XCO_A17']
        #ds['LPrmCO10'] = tb['LCO_COR']
        ds['LPrmCO10'] = 23.5 * 5.0*tb['RMS_CO'] * 300/(2*np.sqrt(2*np.log(2))) * np.sqrt(np.pi) * np.pi/(4*np.log(2)) * (22.0)**2 * (tb['LUMDIST'])**2 * np.power((1.+tb['Z_SDSS']),-3)
        #    ------       TODO: 5.0 sigma upper limits, assumed 300 km/s, 1D Gaussian integral = sqrt(pi) * a * c, c = FWHM / (2sqrt(2ln(2)))
        #    ------             their RMS_CO are measured in a 20km/s channel-width data cube
        mask2 = (tb['Z_PP04_O3N2']>0) # we use their metalZ if valid, otherwise compute from MZR or FMR with calc_metal_Z_local_galaxy_method()
        ds['MetalZ'] = calc_metal_Z_local_galaxy_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['MetalZ'][mask2] = tb['Z_PP04_O3N2'][mask2] # their metalZ are derived from optical nebular lines (SDSS) with PP04 O3N2 calibration 
        ds['alphaCO'] = calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['Mmol'] = (ds['LPrmCO10']) * ds['alphaCO']
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Villanueva+2017'
        ds['color'] = 'orangered'
        ds['facecolor'] = 'orangered'
        ds['edgecolor'] = 'none'
        ds['alpha'] = 0.7
        ds['marker'] = 'x'
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_Villanueva2017_VALES/datatable_Villanueva2017_VALES_Survey.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = np.logical_and(tb['lgMH2']>0, tb['e_lgMH2']>0)
        tb = tb[mask]
        ds['ID'] = tb['ID_GAMA']
        ds['z'] = tb['z_spec']
        ds['SFR'] = np.power(10, tb['lgLIR'])/1e10
        ds['Mstar'] = np.power(10, tb['lgMstar']) # Chabrier 2003 IMF
        #ds['Mmol'] = np.power(10, tb['lgMH2']) # they assumed constant αCO = 4.6
        ds['LPrmCO10'] = 3.25e7 * tb['SCO10'] / 115.271**2 * (cosmo.luminosity_distance(tb['z_spec']).value)**2 / (1.+tb['z_spec'])
        ds['MetalZ'] = calc_metal_Z_local_galaxy_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['alphaCO'] = calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['Mmol'] = tb['LPrmCO10'] * ds['alphaCO']
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        ds['deltaGasErr'] = np.sqrt((tb['e_SCO10']/tb['SCO10'])**2 + (tb['e_lgMstar'])**2) * ds['deltaGas'] # CO obs uncertainty, gas conversion 0 dex scatter (alphaCO), stellar mass obs uncertainty
        ds['tauDeplErr'] = np.sqrt((tb['e_SCO10']/tb['SCO10'])**2 + (tb['e_SFR']/tb['SFR'])**2) * ds['tauDepl'] # CO obs uncertainty, gas conversion 0 dex scatter (alphaCO), SFR obs uncertainty
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Villanueva+2017 uplims'
        ds['color'] = 'orangered'
        ds['facecolor'] = 'orangered'
        ds['edgecolor'] = 'orangered'
        ds['alpha'] = 0.7
        ds['marker'] = 'uplims'
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_Villanueva2017_VALES/datatable_Villanueva2017_VALES_Survey.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = np.logical_and(np.logical_or(tb['lgMH2']<=0, tb['e_lgMH2']<=0), tb['e_SCO10']>0)
        tb = tb[mask]
        ds['ID'] = tb['ID_GAMA']
        ds['z'] = tb['z_spec']
        ds['SFR'] = np.power(10, tb['lgLIR'])/1e10
        ds['Mstar'] = np.power(10, tb['lgMstar']) # Chabrier 2003 IMF
        #ds['Mmol'] = np.power(10, tb['lgMH2'])
        ds['LPrmCO10'] = 3.25e7 * 5.0*tb['e_SCO10'] / 115.271**2 * (cosmo.luminosity_distance(tb['z_spec']).value)**2 / (1.+tb['z_spec'])
        #    ------         TODO: 5.0 sigma upper limits, 
        ds['MetalZ'] = calc_metal_Z_local_galaxy_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['alphaCO'] = calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['Mmol'] = tb['LPrmCO10'] * ds['alphaCO']
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Bertemes+2018'
        ds['color'] = 'lime'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'lime'
        ds['alpha'] = 0.7
        ds['marker'] = '^'
        ds['markersize'] = 15
        #ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_Bertemes2018_Stripe82/datatable_Bertemes2018_Stripe82_reformatted.csv'
        #tb = Table.read(ds['datatable'], format='csv')
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatable_z0_Bertemes2018_Stripe82/datatable_Bertemes2018_Stripe82_big_reformatted_by_dzliu_xmatched_to_SDSS_DR7_MPA_JHU_catalog.fits' # 'datatables_z_deltaGas'+os.sep+ <TODO>
        tb = Table.read(ds['datatable'], format='fits')
        ds['ID'] = tb['ID']
        ds['z'] = tb['z']
        ds['SFR'] = np.power(10.0, tb['logSFR'])
        #ds['Mstar'] = np.power(10.0, tb['logMstar']) # Chabrier 2003 IMF
        ds['Mstar'] = np.power(10.0, tb['AVG_LOGMSTAR']) # here I replaced with MPA-JHU stellar mass which is about 0.1 dex larger, and assumed Kroupa IMF. Saintonge+2017 Mstar are exactly from MPA-JHU stellar mass AVG_LOGMSTAR. 
        #ds['Mmol'] = np.power(10, tb['logMgas_via_CO']) # they already have metallicity-dependent alphaCO, but ...
        ds['MetalZ'] = tb['MetalZ'] # their metalZ are derived from optical nebular lines (SDSS) with PP04 O3N2 calibration # calc_metal_Z_local_galaxy_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['Mmol'] = np.power(10.0, tb['logMgas_via_CO']) / calc_alphaCO_from_metalZ_following_Bertemes2018(tb['MetalZ']) * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        # my calc_metalZ_from_FMR_with_dzliu_selection() metallicity are in general lower than their tb['metalZ'], so the inferred alphaCO is even higher...
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        ds['deltaGasErr'] = np.sqrt((tb['E_logMgas_via_CO']/tb['logMgas_via_CO'])**2 + (0.2)**2) * ds['deltaGas'] # CO obs uncertainty, gas conversion 0 dex scatter (alphaCO), stellar mass 0.2 dex uncertainty
        ds['tauDeplErr'] = np.sqrt((tb['E_logMgas_via_CO']/tb['logMgas_via_CO'])**2 + (0.1)**2) * ds['tauDepl'] # CO obs uncertainty, gas conversion 0 dex scatter (alphaCO), SFR 0.1 dex uncertainty
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Lee+2017'
        ds['color'] = 'green'
        ds['facecolor'] = 'green'
        ds['edgecolor'] = 'none'
        ds['alpha'] = 0.7
        ds['marker'] = '2' # tri_up
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_Lee2017/Lee2017_CO32_z0p5_sample.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        ds['ID'] = tb['ID']
        ds['z'] = tb['z']
        ds['SFR'] = tb['SFR']
        ds['Mstar'] = np.power(10, tb['lgMstar']) # Chabrier 2003 IMF
        ds['MetalZ'] = ds['z'] * 0.0 - 99 #<TODO># no metallicity info
        ds['Mmol'] = tb['deltaGas'] * ds['Mstar'] # they used Genzel+2015 alphaCO
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        ds['deltaGasErr'] = np.sqrt((tb['e_deltaGas']/tb['deltaGas'])**2 + (tb['e_lgMstar'])**2) * ds['deltaGas'] # CO obs uncertainty, gas conversion 0 dex scatter (alphaCO), stellar mass obs uncertainty
        ds['tauDeplErr'] = np.sqrt((tb['e_deltaGas']/tb['deltaGas'])**2 + (tb['e_lgLIR'])**2) * ds['tauDepl'] # CO obs uncertainty, gas conversion 0 dex scatter (alphaCO), SFR obs uncertainty
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Bauermeister+2013'
        ds['color'] = 'royalblue'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'royalblue'
        ds['alpha'] = 0.7
        ds['marker'] = '<'
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z0.5_Bauermeister2013/datatable_Bauermeister2013_EGNOG_Survey.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = (tb['Mgas']>3.0*tb['e_Mgas'])
        tb = tb[mask]
        ds['ID'] = tb['Source_Name']
        ds['z'] = tb['z']
        ds['SFR'] = tb['SFR']
        ds['Mstar'] = np.power(10.0, tb['lgMstar']) # SDSS DR7, MPA-JHU stellar mass catalog, Kauffmann+2003 - Kroupa 2001 IMF
        #ds['Mmol'] = tb['Mgas'] # bimodal alphaCO: 1.36 * 3.2 for normal galaxies, 1.36 * 0.8(?) for starburst galaxies. 
        ds['MetalZ'] = calc_metal_Z_local_galaxy_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['alphaCO'] = calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['alphaCO_old'] = tb['Mgas'] * 0.0 + 1.36*3.2 # see their Eq.(4) and text after Eq.(4)
        ds['Mmol'] = (tb['Mgas'] / ds['alphaCO_old']) * ds['alphaCO']
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        ds['deltaGasErr'] = np.sqrt((tb['e_Mgas']/tb['Mgas'])**2 + (np.max([tb['e_lgMstar_hi'],-tb['e_lgMstar_lo']],axis=0))**2) * ds['deltaGas'] # CO obs uncertainty, gas conversion 0 dex scatter (alphaCO), stellar mass obs uncertainty
        ds['tauDeplErr'] = np.sqrt((tb['e_Mgas']/tb['Mgas'])**2 + (np.max([tb['e_SFR_hi'],-tb['e_SFR_lo']],axis=0)/tb['SFR'])**2) * ds['tauDepl'] # CO obs uncertainty, gas conversion 0 dex scatter (alphaCO), SFR obs uncertainty
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
        # --
        # -- here to check the differences -- max 0.39 dex
        #print('')
        #dsprint = {}
        #dsprint['lgMmol'] = np.log10(ds['Mmol'])
        #dsprint['lgMgas'] = np.log10(tb['Mgas'])
        #dsprint['difference'] = (dsprint['lgMmol'] - dsprint['lgMgas'])
        #dsprint['alphaCO'] = ds['alphaCO']
        #dsprint['alphaCO_old'] = ds['alphaCO_old']
        #tbprint = Table(dsprint)
        #print(tbprint)
        #print('')
    
    if True:
        ds = {}
        ds['label'] = 'Cicone+2017' # 2017A&A...604A..53C
        ds['color'] = 'cyan'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'cyan'
        ds['alpha'] = 0.7
        ds['marker'] = 's'
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z0_Cicone2017/Cicone2017_combined_table_by_dzliu.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask2 = np.logical_or(tb['Sample']=='A',tb['Sample']=='I2')
        tb['LCO_CORR'][mask2] = tb['LCO_CORR'][mask2] # / 0.8 # convert CO(2-1) to CO(1-0) #<TODO>#
        tb['e_LCO_CORR'][mask2] = tb['e_LCO_CORR'][mask2] # / 0.8 # convert CO(2-1) to CO(1-0) #<TODO>#
        mask = (tb['flag_LCO_CORR']!='<')
        tb = tb[mask]
        ds['ID'] = tb['ID']
        ds['z'] = tb['z_SDSS']
        ds['SFR'] = 10**(tb['lgSFR']) # MPA-JHU catalog
        ds['Mstar'] = 10**(tb['lgMstar']) # MPA-JHU catalog
        #ds['MetalZ'] = calc_metal_Z_local_galaxy_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['MetalZ'] = tb['Z_PP04_O3N2']
        ds['Mmol'] = (tb['LCO_CORR']) * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        ds['deltaGasErr'] = np.sqrt((tb['e_LCO_CORR']/tb['LCO_CORR'])**2 + (tb['e_lgMstar'])**2) * ds['deltaGas'] # CO obs uncertainty, gas conversion 0 dex scatter (alphaCO), stellar mass obs uncertainty
        ds['tauDeplErr'] = np.sqrt((tb['e_LCO_CORR']/tb['LCO_CORR'])**2 + (tb['e_lgSFR'])**2) * ds['tauDepl'] # CO obs uncertainty, gas conversion 0 dex scatter (alphaCO), SFR obs uncertainty
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Cicone+2017 uplims' # 2017A&A...604A..53C
        ds['color'] = 'cyan'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'cyan'
        ds['alpha'] = 0.7
        ds['marker'] = 'uplims'
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z0_Cicone2017/Cicone2017_combined_table_by_dzliu.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask2 = np.logical_or(tb['Sample']=='A',tb['Sample']=='I2')
        tb['LCO_CORR'][mask2] = tb['LCO_CORR'][mask2] # / 0.8 # convert CO(2-1) to CO(1-0)
        tb['e_LCO_CORR'][mask2] = tb['e_LCO_CORR'][mask2] # / 0.8 # convert CO(2-1) to CO(1-0)
        mask = (tb['flag_LCO_CORR']=='<')
        tb = tb[mask]
        ds['ID'] = tb['ID']
        ds['z'] = tb['z_SDSS']
        ds['SFR'] = 10**(tb['lgSFR']) # MPA-JHU catalog
        ds['Mstar'] = 10**(tb['lgMstar']) # MPA-JHU catalog
        #ds['MetalZ'] = calc_metal_Z_local_galaxy_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['MetalZ'] = tb['Z_PP04_O3N2']
        ds['Mmol'] = (tb['e_LCO_CORR']/3.0*5.0) * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ']) # 3-sigma uplims -> 5-sigma uplims
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Combes+2013' # Combes 2013 - 2013A&A...550A..41C
        ds['color'] = 'none'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'magenta'
        ds['alpha'] = 0.7
        ds['marker'] = 'd' # thin_diamond
        ds['markersize'] = 12
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z0.8_Combes2013/Combes2013_combined_table_by_dzliu.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = np.logical_and(np.logical_and(tb['flag_SCO']!='<', tb['lgMstar']>0), tb['LineName']=='CO(2-1)')
        tb = tb[mask]
        ds['ID'] = tb['ID']
        ds['z'] = tb['z']
        ds['SFR'] = np.power(10, tb['lgLIR'])/1e10
        ds['Mstar'] = np.power(10, tb['lgMstar'] - 0.238) # Salpeter 1955 IMF
        ds['LPrmCO10'] = 3.25e7 * tb['SCO'] / 230.538**2 * (cosmo.luminosity_distance(ds['z']).value)**2 / (1.+ds['z']) / 0.8 # dividing 0.8 convert LPrmCO21 to LPrmCO10
        ds['MetalZ'] = calc_metal_Z_local_galaxy_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['alphaCO'] = calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['Mmol'] = ds['LPrmCO10'] * ds['alphaCO']
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Combes+2013 uplims' # Combes 2013 - 2013A&A...550A..41C
        ds['color'] = 'magenta'
        ds['facecolor'] = 'magenta'
        ds['edgecolor'] = 'magenta'
        ds['alpha'] = 0.7
        ds['marker'] = 'uplims'
        ds['markersize'] = 12
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z0.8_Combes2013/Combes2013_combined_table_by_dzliu.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = np.logical_and(np.logical_and(tb['flag_SCO']=='<', tb['lgMstar']>0), tb['LineName']=='CO(2-1)')
        tb = tb[mask]
        ds['ID'] = tb['ID']
        ds['z'] = tb['z']
        ds['SFR'] = np.power(10, tb['lgLIR'])/1e10
        ds['Mstar'] = np.power(10, tb['lgMstar'] - 0.238) # Salpeter 1955 IMF
        ds['LPrmCO10'] = 3.25e7 * (tb['e_SCO']/3.0*5.0) / 230.538**2 * (cosmo.luminosity_distance(ds['z']).value)**2 / (1.+ds['z']) / 0.8 # dividing 0.8 convert LPrmCO21 to LPrmCO10
        #    ------       TODO: 3.0 sigma -> 5.0 sigma upper limits, assumed 300 km/s except for ID 16 which has known DV according to the paper
        ds['MetalZ'] = calc_metal_Z_local_galaxy_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['alphaCO'] = calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['Mmol'] = ds['LPrmCO10'] * ds['alphaCO']
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Kirkpatrick+2014' # 
        ds['color'] = 'lightgreen'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'lightgreen'
        ds['alpha'] = 0.7
        ds['marker'] = '>'
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z0.1_Kirkpatrick2014/datatable_Kirkpatrick2014.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = (tb['LPrmCO10']>0)
        tb = tb[mask]
        ds['ID'] = tb['Source_Name']
        ds['z'] = tb['z']
        ds['SFR'] = np.power(10.0, tb['lgLIR_SF']) / 1e10
        ds['Mstar'] = np.power(10, tb['lgMstar']) # 
        ds['MetalZ'] = calc_metal_Z_local_galaxy_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['Mmol'] = tb['LPrmCO10'] * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Cortzen+2019' # 
        ds['color'] = 'cyan'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'cyan'
        ds['alpha'] = 0.7
        ds['marker'] = '$c$'
        ds['markersize'] = 12
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z0.1_Cortzen2019/datatable_Cortzen2019_with_deltaMS.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = np.logical_and(tb['LPrmCO10']>0, tb['eLPrmCO10']>0)
        tb = tb[mask]
        ds['ID'] = tb['Source_Name']
        ds['z'] = tb['z']
        ds['SFR'] = np.power(10.0, tb['lgLIR']) / 1e10
        ds['Mstar'] = np.power(10, tb['lgMstar']) # 
        ds['MetalZ'] = calc_metal_Z_local_galaxy_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['Mmol'] = tb['LPrmCO10'] * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Spilker+2018'
        ds['color'] = 'red'
        ds['facecolor'] = 'red'
        ds['edgecolor'] = 'none'
        ds['alpha'] = 0.7
        ds['marker'] = '1' # tri_down
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_Spilker2018/Spilker2018_CO21_z0p7_sample.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = (tb['Mmol']>3.0*tb['eMmol'])
        tb = tb[mask]
        ds['ID'] = tb['ID']
        ds['z'] = tb['z']
        ds['SFR'] = tb['SFR']
        ds['Mstar'] = tb['Mstar'] # Chabrier 2003 IMF
        ds['MetalZ'] = ds['z'] * 0.0 - 99 #<TODO># no metallicity info
        ds['Mmol'] = tb['Mmol']
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Lisenfeld+2017' # 2017A&A...607A.110L
        ds['color'] = 'pink'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'pink'
        ds['alpha'] = 0.7
        ds['marker'] = 'D'
        ds['markersize'] = 12
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z0_Lisenfeld2018/tab2_lourdes_reformatted_by_dzliu.dat'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = (tb['det_co_code']>0)
        tb = tb[mask]
        ds['ID'] = tb['Name']
        ds['z'] = tb['z']
        ds['SFR'] = tb['sfr']
        ds['Mstar'] = np.power(10, tb['log_Mstar']) # Kroupa IMF
        #ds['Mmol'] = np.power(10, tb['log_Mh2']) # assumed X_CO_MW = 2.0e20, i.e., alphaCO=4.3
        ds['MetalZ'] = calc_metal_Z_local_galaxy_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['Mmol'] = np.power(10, tb['log_Mh2']) / 4.3 * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if False:
        ds = {}
        ds['label'] = 'Scoville+2016' # 
        ds['color'] = 'red'
        ds['facecolor'] = 'red'
        ds['edgecolor'] = 'none'
        ds['alpha'] = 0.7
        ds['marker'] = 'x'
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_Scoville2016/datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = (tb['Mstar']>0)
        tb = tb[mask]
        ds['ID'] = tb['ID']
        ds['z'] = tb['z']
        ds['SFR'] = tb['SFR']
        ds['Mstar'] = tb['Mstar'] # 
        ds['MetalZ'] = ds['z'] * 0.0 - 99 #<TODO># no metallicity info
        ds['Mmol'] = tb['deltaGas'] * tb['Mstar']
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Magdis+2012a' # 
        ds['color'] = 'lightgreen'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'lightgreen'
        ds['alpha'] = 0.7
        ds['marker'] = '^'
        ds['markersize'] = 35
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z1.5_Magdis2012/datatable_Magdis2012.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = (tb['Mstar']>0)
        tb = tb[mask]
        ds['ID'] = tb['Source']
        ds['z'] = tb['z']
        ds['SFR'] = tb['SFR'] # Chabrier 2003 IMF
        ds['Mstar'] = tb['Mstar'] # Chabrier 2003 IMF
        ds['MetalZ'] = ds['z'] * 0.0 - 99 #<TODO># no metallicity info
        ds['Mmol'] = tb['Mgas'] # already had their optimized metallicity-dependent alphaCO (~2.4-13.4)
        #ds['MetalZ'] = calc_metal_Z_high_z_method(ds['Mstar'], ds['SFR'], ds['z'])
        #ds['Mmol'] = tb['LPrmCO10'] * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Magdis+2012b' # 
        ds['color'] = 'green'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'green'
        ds['alpha'] = 0.7
        ds['marker'] = 'v'
        ds['markersize'] = 35
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z3.0_Magdis2012b/datatable_Magdis2012b.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = (tb['Mstar']>0)
        tb = tb[mask]
        ds['ID'] = tb['Source']
        ds['z'] = tb['z']
        ds['SFR'] = tb['SFR'] # Chabrier 2003 IMF
        ds['Mstar'] = tb['Mstar'] # Chabrier 2003 IMF
        ds['MetalZ'] = ds['z'] * 0.0 - 99 #<TODO># no metallicity info
        ds['Mmol'] = tb['Mgas'] # already had their optimized metallicity-dependent alphaCO (~2.4-13.4)
        #ds['MetalZ'] = calc_metal_Z_high_z_method(ds['Mstar'], ds['SFR'], ds['z'])
        #ds['Mmol'] = tb['LPrmCO10'] * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Magdis+2017' # 
        ds['color'] = 'darkgreen'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'darkgreen'
        ds['alpha'] = 0.7
        ds['marker'] = '>'
        ds['markersize'] = 35
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z3.0_Magdis2017/datatable_Magdis2017.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = (tb['Mstar']>0)
        tb = tb[mask]
        ds['ID'] = tb['Source']
        ds['z'] = tb['z']
        ds['SFR'] = tb['SFR'] # Chabrier 2003 IMF
        ds['Mstar'] = tb['Mstar'] # Chabrier 2003 IMF
        ds['MetalZ'] = ds['z'] * 0.0 - 99 #<TODO># no metallicity info
        ds['Mmol'] = tb['Mgas'] # already had their optimized metallicity-dependent alphaCO (~2.4-13.4)
        #ds['MetalZ'] = calc_metal_Z_high_z_method(ds['Mstar'], ds['SFR'], ds['z'])
        #ds['Mmol'] = tb['LPrmCO10'] * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Kaasinen+2019' # 
        ds['color'] = 'magenta'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'magenta'
        ds['alpha'] = 0.7
        ds['marker'] = 'H' # hexagon2
        ds['markersize'] = 18
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z2_Kaasinen2019/Kaasinen2019_CO10_z2_sample.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = np.logical_and(np.logical_and(np.logical_and(tb['Mstar']>0, tb['SFR']>0), tb['LPrmCO10']>0), tb['z_CO10']>0)
        tb = tb[mask]
        ds['ID'] = tb['ID']
        ds['z'] = tb['z_CO10']
        ds['SFR'] = tb['SFR']
        ds['Mstar'] = tb['Mstar']
        ds['MetalZ'] = calc_metal_Z_high_z_method(ds['Mstar'], ds['SFR'], ds['z'])
        ds['Mmol'] = tb['LPrmCO10'] * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ']) # They assumed alphaCO=6.5. They observed CO(1-0).
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9 # Gyr
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'Tan+2014' # 
        ds['color'] = 'cyan'
        ds['facecolor'] = 'cyan'
        ds['edgecolor'] = 'cyan'
        ds['alpha'] = 0.7
        ds['marker'] = 'd' # thin_diamond
        ds['markersize'] = 25
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z4.0_Tan2014/datatable_Tan2014.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = (tb['Mstar']>0)
        tb = tb[mask]
        ds['ID'] = tb['Source']
        ds['z'] = tb['z']
        ds['SFR'] = 10**(tb['lgLIR']) / 1e10 # Chabrier 2003 IMF
        ds['Mstar'] = tb['Mstar'] # Chabrier 2003 IMF
        ds['MetalZ'] = ds['z'] * 0.0 - 99 #<TODO># no metallicity info
        ds['Mmol'] = tb['MH2']
        #ds['MetalZ'] = calc_metal_Z_high_z_method(ds['Mstar'], ds['SFR'], ds['z'])
        #ds['Mmol'] = tb['LPrmCO10'] * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if False:
        ds = {}
        ds['label'] = 'Schinnerer+2016' # 
        ds['color'] = 'magenta'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'magenta'
        ds['alpha'] = 0.7
        ds['marker'] = '*'
        ds['markersize'] = 55
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z3.0_Schinnerer2016/datatable_Schinnerer2016_averaged.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header') # this table contains: z sSFR deltaGas deltaGas_err tauDepl
        ds['ID'] = np.array(['0'])
        ds['z'] = tb['z']
        ds['SFR'] = 10**(tb['lgSFR']) # according to the paper abstract, the sample has an averaged SFR of 2.0 in log10. Chabrier 2003 IMF
        ds['Mstar'] = 10**(tb['lgMstar']) # according to the paper abstract, the sample has an averaged stellar mass of 10.7 in log10. Chabrier 2003 IMF
        ds['MetalZ'] = ds['z'] * 0.0 - 99 #<TODO># no metallicity info
        ds['Mmol'] = tb['deltaGas'] * ds['Mstar'] # revert deltaGas (==\mu_{gas} in the paper) to Mmol
        #ds['MetalZ'] = calc_metal_Z_high_z_method(ds['Mstar'], ds['SFR'], ds['z'])
        #ds['Mmol'] = tb['LPrmCO10'] * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['sSFR'] = tb['sSFR']
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = tb['deltaGas']
        ds['tauDepl'] = tb['tauDepl']
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if False:
        ds = {}
        ds['label'] = 'Liu+ (GOODSN 1mm sample)' # dzliu NOEMA 1mm
        ds['color'] = 'magenta'
        ds['facecolor'] = 'magenta'
        ds['edgecolor'] = 'magenta'
        ds['alpha'] = 0.8
        ds['marker'] = 'o'
        ds['markersize'] = 55
        ds['datatable'] = '/Users/dzliu/Work/DeepFields/Works/GOODSN_z35_1mm/Run_SED_Fitting/Calc_Gas_Mass/datatable_ID_RA_Dec_z_Mstar_SFR_sSFR_with_deltaGas.fits'
        tb = Table.read(ds['datatable'], format='fits') # this table contains: z sSFR deltaGas deltaGas_err tauDepl
        ds['ID'] = tb['ID']
        ds['z'] = tb['z']
        ds['SFR'] = tb['SFR']
        ds['Mstar'] = tb['Mstar']
        ds['Mmol'] = tb['deltaGas'] * ds['Mstar']
        ds['sSFR'] = tb['sSFR']
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = tb['deltaGas']
        ds['tauDepl'] = 1.0/tb['sSFR'] # Gyr
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'DGS' # DGS, Remy-Ruyer+2014, 2015
        ds['color'] = 'magenta'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'magenta'
        ds['alpha'] = 0.7
        ds['marker'] = '*'
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatable_z0_DGS/RemyRuyer,2015A&A...582A.121R,Mdust,20190206.csv'
        #tb = Table.read(ds['datatable'], format='csv')
        tb = pd.read_csv(ds['datatable'])
        mask = np.logical_and(np.logical_and(np.logical_and(tb['Mstar']>0, tb['L_IR']>0), tb['MH2']>0), tb['Metallicity']>0)
        tb = tb[mask]
        ds['ID'] = tb['Source'].values
        ds['z'] = tb['L_IR'].values*0.0 + 0.00001
        ds['SFR'] = tb['L_IR'].values / 1e10
        ds['Mstar'] = tb['Mstar'].values
        ds['MetalZ'] = tb['Metallicity'].values
        #ds['Mmol'] = tb['MH2'].values # directly take their metallicity-dependent MH2 -- this leads to some too high tauDepl
        ds['Mmol'] = tb['MH2_L68'].values / 3.16 * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9 # Gyr
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'HRS' # HRS
        ds['color'] = 'magenta'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'magenta'
        ds['alpha'] = 0.7
        ds['marker'] = 'v' # triangle_down
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatable_z0_HRS/datatable_Andreani2018_Table1_with_Hughes2013_metallicity.fits'
        tb = Table.read(ds['datatable'], format='fits')
        mask = np.logical_and(np.logical_and(np.logical_and(tb['lgM_star']>0, tb['lgL_IR']>0), tb['lgM_H2']>0), tb['Metallicity']>0)
        tb = tb[mask]
        ds['ID'] = tb['ID']
        ds['z'] = (10**tb['lgL_IR'])*0.0 + 0.003600 # Virgo Cluster redshift from NED
        ds['SFR'] = (10**tb['lgL_IR']) / 1e10
        ds['Mstar'] = (10**tb['lgM_star'])
        ds['MetalZ'] = tb['Metallicity']
        ds['Mmol'] = (10**tb['lgM_H2']) / 3.6 * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9 # Gyr
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    if True:
        ds = {}
        ds['label'] = 'KINGFISH' # Groves+2015 KINGFISH
        ds['color'] = 'magenta'
        ds['facecolor'] = 'none'
        ds['edgecolor'] = 'magenta'
        ds['alpha'] = 0.7
        ds['marker'] = 'p' # pentagon
        ds['markersize'] = 15
        ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatable_z0_KINGFISH/datatable_Groves2015_FLUX_500.txt'
        tb = Table.read(ds['datatable'], format='ascii.commented_header')
        mask = np.logical_and(np.logical_and(np.logical_and(tb['logMstar']>0, tb['SFR']>0), tb['MH2']>0), tb['Metallicity_KK04']>0)
        tb = tb[mask]
        ds['ID'] = tb['Name_KINGFISH']
        ds['z'] = tb['z']
        ds['SFR'] = tb['SFR']
        ds['Mstar'] = (10**tb['logMstar'])
        ds['MetalZ'] = convert_metalZ_KK04_to_metalZ_PP04(tb['Metallicity_KK04'])
        ds['Mmol'] = tb['MH2'] / 4.4 * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ']) # they assumed alphaCO=4.4, R21=0.8
        ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
        ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
        ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
        ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
        ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9 # Gyr
        datasets.append(ds)
        print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    #if True:
    #    ds = {}
    #    ds['label'] = 'Jiang+2015' # Jiang 2015 - CO in local galaxies - 2015ApJ...799...92J.pdf
    #    ds['color'] = 'cyan'
    #    ds['facecolor'] = 'none'
    #    ds['edgecolor'] = 'cyan'
    #    ds['alpha'] = 0.7
    #    ds['marker'] = 'x'
    #    ds['markersize'] = 15
    #    ds['datatable'] = os.path.dirname(os.path.dirname(__file__))+os.sep+'Data_Tables'+os.sep+'datatables_z_deltaGas'+os.sep+'datatable_z0_Jiang2015/Jiang2015_table2.txt'
    #    tb = Table.read(ds['datatable'], format='fits')
    #    mask = np.logical_and(np.logical_and(np.logical_and(tb['lgMstar']>0, tb['lgLIR']>0), tb['lgMH2']>0), tb['Metallicity']>0)
    #    tb = tb[mask]
    #    ds['ID'] = tb['ID']
    #    ds['z'] = (10**tb['lgL_IR'])*0.0 + 0.003600 # Virgo Cluster redshift from NED
    #    ds['SFR'] = (10**tb['lgL_IR']) / 1e10
    #    ds['Mstar'] = (10**tb['lgM_star'])
    #    ds['MetalZ'] = tb['Metallicity']
    #    ds['Mmol'] = (10**tb['lgM_H2']) / 3.6 * calc_alphaCO_from_metalZ_following_Tacconi2018(ds['MetalZ'])
    #    ds['sSFR'] = ds['SFR'] / ds['Mstar'] * 1e9
    #    ds['sSFR_MS'] = calc_Speagle2014_sSFR(cosmo.age(ds['z']).value, np.log10(ds['Mstar']))
    #    ds['DeltaMS'] = np.log10(ds['sSFR'] / ds['sSFR_MS'])
    #    ds['deltaGas'] = ds['Mmol'] / ds['Mstar']
    #    ds['tauDepl'] = ds['Mmol'] / ds['SFR'] / 1e9 # Gyr
    #    datasets.append(ds)
    #    print('Read "%s" (%d data)'%(ds['datatable'], len(ds['deltaGas'])))
    
    # 
    # print some key info
    print('Dataset & z_range & lgMstar_range & N_data')
    for ds in datasets:
        # fix Masked Array problem
        for key in ds:
            if type(ds[key]) in [np.ma.core.MaskedArray, MaskedColumn]:
                ds[key] = np.array(ds[key].tolist()) # fix Masked Array problem
            if type(ds[key]) in [Column]:
                ds[key] = ds[key].data
        # compute lgMstar
        ds['lgMstar'] = np.log10(ds['Mstar'])
        print('%s & %0.2f -- %0.2f & %0.2f -- %0.2f & %d \\\\'%(ds['label'], np.min(ds['z']), np.max(ds['z']), np.min(ds['lgMstar']), np.max(ds['lgMstar']), len(ds['deltaGas'])))
    
    
    
    # 
    return datasets









def common_meta_info_columns():
    return ['label', 'color', 'facecolor', 'edgecolor', 'alpha', 'marker', 'markersize', 'datatable']







def common_data_array_columns():
    return ['ID', 'z', 'SFR', 'Mstar', 'MetalZ', 'Mmol', 'sSFR', 'sSFR_MS', 'DeltaMS', 'deltaGas', 'tauDepl']







def merge_datasets(datasets, nouplims=False, savetofile=''):
    # 
    if type(datasets) is dict:
        dataset = datasets
        return dataset
    # 
    new_dataset = {}
    new_dataset['labels'] = []
    meta_info_columns = []
    data_array_columns = []
    # 
    isfirst = True
    for i,dataset in enumerate(datasets):
        if not (type(dataset) is dict):
            raise Exception('Error! Each input dataset should be a dict!')
        # skip uplims dataset
        if nouplims == True:
            if dataset['label'].find('uplims') >=0 :
                continue
        # array_length
        array_length = 0
        # 
        for t in dataset:
            # skip empty-value keys
            if dataset[t] is None:
                #print('dataset[t] is None', dataset['label'], t) # edgelinewidth
                continue
            # check scalar meta info or array-like data array
            if np.isscalar(dataset[t]):
                if isfirst:
                    meta_info_columns.append(t) #<TODO># copy meta only from the first dataset
                    pass # do not copy scalars
                    #new_dataset[t] = dataset[t] 
            else:
                # skip empty-value keys
                # single-item array counted as meta as well
                #if len(dataset[t]) == 1
                #else
                # 
                # copy or concatenate
                if isfirst:
                    data_array_columns.append(t)
                    new_dataset[t] = copy.copy(dataset[t])
                    array_length = len(dataset[t])
                else:
                    if t in data_array_columns:
                        new_dataset[t] = np.concatenate( (new_dataset[t], dataset[t]) )
                        array_length = len(dataset[t])
                    else:
                        print('merge_datasets(): Warning! Skipping column "%s" in dataset "%s"'%(t, dataset['label']))
        # 
        if array_length > 0:
            #print('merge_datasets(): Adding %d data with label %s'%(array_length, dataset['label']))
            new_dataset['labels'].extend(np.repeat(dataset['label'], array_length))
        # 
        isfirst = False
    # 
    # save to file if the user has provided 'savetofile'
    if savetofile != '':
        if os.path.isfile(savetofile):
            shutil.move(savetofile, savetofile+'.backup')
        #for t in new_dataset.keys():
        #    print(t, type(new_dataset[t][0]), len(new_dataset[t]))
        dump_data_table = Table(new_dataset)
        for t in dump_data_table.colnames:
            if dump_data_table[t].dtype.kind == 'O':
                dump_data_table[t] = [str(tval) for tval in dump_data_table[t].data.tolist()]
            #print(t, type(dump_data_table[t][0]), dump_data_table[t].dtype, len(dump_data_table[t]))
        dump_data_table.write(savetofile, format='fits')
        print('merge_datasets(): Dumped to "%s"!'%(savetofile))
    # 
    return new_dataset







def write_datasets_to_file(datasets, output_filename):
    if os.path.isfile(output_filename):
        print('Found existing "%s"! Backing it up as "%s"!'%(output_filename, output_filename+'.backup'))
        shutil.move(output_filename, output_filename+'.backup')
        if os.path.isfile(output_filename+'.meta.txt'):
            shutil.move(output_filename+'.meta.txt', output_filename+'.meta.txt'+'.backup')
    # 
    if type(datasets) is dict:
        datasets = [datasets]
    # 
    for dataset in datasets:
        if not (type(dataset) is dict):
            raise Exception('Error! Each input dataset should be a dict!')
        # 
        meta_info = {}
        data_arrays = {}
        for t in dataset:
            if not np.isscalar(dataset[t]):
                if len(dataset[t]) > 1:
                    data_arrays[t] = dataset[t]
                else:
                    meta_info[t] = dataset[t]
            else:
                meta_info[t] = dataset[t]
    # 
    output_data_table = Table(data_arrays)
    if output_filename.endswith('fits') or output_filename.endswith('FITS'):
        output_data_table.meta = meta_info
        output_data_table.write(output_filename, format='fits')
    else:
        output_data_table.write(output_filename, format='ascii.fixed_width', delimiter='  ', bookend=True)
        with open(output_data_table, 'r') as fp:
            fp.seek(0)
            fp.write('#')
    print('Output to "%s"!'%(output_filename))
    # 
    meta_info['datetime'] = datetime.datetime.now().strftime('%Y-%m-%d %Hh%Mm%Ss')+' '+time.localtime().tm_zone
    with open(output_data_table+'.meta.txt', 'w') as fp:
        json.dump(meta_info, fp, sort_keys=True, indent=4)
    print('Output to "%s"!'%(output_filename+'.meta.txt'))






def count_z_number():
    datasets = read_datasets()
    total_complementary_sample = 0
    low_z_complementary_sample = 0
    for dataset in datasets:
        if dataset['label'].find('A3COSMOS')>=0:
            continue
        low_z_complementary_sample += np.count_nonzero(dataset['z']<0.1)
        total_complementary_sample += len(dataset['z'])
    print('total complementary sample = %d'%(total_complementary_sample))
    print('z < 0.1 = %d (%0.2f%%)'%(low_z_complementary_sample, 100.0*low_z_complementary_sample/total_complementary_sample))














