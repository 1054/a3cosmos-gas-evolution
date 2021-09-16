

from __future__ import print_function # for Python<=2.X
#from __future__ import with_statement # for Python<=2.5
#from __future__ import nested_scopes # for Python<=2.1
#from pkg_resources import parse_version
from distutils.version import LooseVersion, StrictVersion
import os, sys, re
if sys.version_info.major>=3:
    import inspect
    inspect_func_signature = inspect.signature
else:
    import funcsigs
    inspect_func_signature = funcsigs.signature



# 
# import functions (keep these the same as in __init__.py)
# 
from .Common_Python_Code import get_cosmic_mol_gas_density
from .Common_Python_Code import calc_gas_fraction
from .Common_Python_Code import calc_galaxy_main_sequence
from .Common_Python_Code import calc_gas_depletion_time
from .Common_Python_Code import calc_alpha_CO
from .Common_Python_Code import calc_delta_GD
from .Common_Python_Code import calc_metal_Z
from .Common_Python_Code import calc_fmol
from .Common_Python_Code import calc_gas_mass_from_dust
from .Common_Python_Code import apply_cosmology
cosmo = apply_cosmology.cosmo




# 
# load all modules
# 
def load_all_modules():
    load_all_functions_in_one_module(get_cosmic_mol_gas_density,  r'^get_cosmic_mol_gas_density_.*',                        print_to_screen = False)
    load_all_functions_in_one_module(get_cosmic_mol_gas_density,  r'^plot_cosmic_mol_gas_density.*',                        print_to_screen = False)
    load_all_functions_in_one_module(calc_gas_fraction,           r'^calc_gas_fraction_.*',                                 print_to_screen = False)
    load_all_functions_in_one_module(calc_galaxy_main_sequence,   r'^calc_SFR_MS_.*',                                       print_to_screen = False)
    load_all_functions_in_one_module(calc_gas_depletion_time,     r'^calc_gas_depletion_time_.*',                           print_to_screen = False)
    load_all_functions_in_one_module(calc_alpha_CO,               r'^calc_alphaCO_from_metalZ_following_.*',                print_to_screen = False)
    load_all_functions_in_one_module(calc_delta_GD,               r'^calc_deltaGD_from_metalZ_following_.*',                print_to_screen = False)
    load_all_functions_in_one_module(calc_metal_Z,                r'^(calc_metalZ_from_FMR_following_|convert_metalZ_).*',  print_to_screen = False)
    load_all_functions_in_one_module(calc_fmol,                   r'^(calc_fmol_from_metalZ_following_).*',                 print_to_screen = False)
    load_all_functions_in_one_module(calc_gas_mass_from_dust,     r'^(calc_gas_mass_from_dust_).*',                         print_to_screen = False)




# 
# load and print all functions in one module
# 
def load_all_functions_in_one_module(t_module, func_name_pattern = '', print_to_screen = False):
    # get function list in the input t_module
    t_func_list = [ t_func_name \
                    for t_func_name in dir(t_module) \
                        if re.match(func_name_pattern, t_func_name, re.IGNORECASE) \
                            and hasattr(getattr(t_module, t_func_name), '__call__') ]
    # sort function names
    t_func_list.sort(key=LooseVersion)
    # prepare printing prefix
    t_print_prefix = '    '
    # get function name maximum string length
    t_func_nchar = 0
    for t_func_name in t_func_list:
        t_func_nchar = max([len(t_func_name), t_func_nchar])
    # loop function names and print function arguments
    for t_func_name in t_func_list:
        if print_to_screen == True:
            print(t_print_prefix, end='')
            print(('%%-%ds'%(t_func_nchar))%(t_func_name), end='')
        t_func_self = getattr(t_module,t_func_name)
        # 
        # The commented code below are for getting t_func_args without import inspect. The problem is len(t_func_defaults) != t_func_narg. 
        #if sys.version_info.major < 2 or (sys.version_info.major==2 and sys.version_info.minor<=5):
        #    t_func_narg = t_func_self.func_code.co_argcount
        #    t_func_args = t_func_self.func_code.co_varnames
        #    t_func_defaults = t_func_self.func_defaults
        #else:
        #    t_func_narg = t_func_self.__code__.co_argcount
        #    t_func_args = t_func_self.__code__.co_varnames
        #    t_func_defaults = t_func_self.__defaults__
        # 
        # Get t_func_args with inspect. Still has the problem that len(t_func_defaults) != t_func_narg. 
        #t_func_spec = inspect.getargspec(t_func_self) # inspect.getargspec() is deprecated since Python 3.0, use inspect.signature() or inspect.getfullargspec()
        #t_func_spec = inspect.getfullargspec(t_func_self)
        #t_func_narg = len(t_func_spec.args)
        #t_func_args = t_func_spec.args
        #t_func_defaults = t_func_spec.defaults
        # 
        # Loop function arguments and print
        #for t_func_iarg in range(t_func_narg):
        #    if t_func_iarg == 0: print('(', end='')
        #    else: print(', ', end='')
        #    print(t_func_args[t_func_iarg], end='')
        #    if t_func_defaults[t_func_iarg] != None:
        #        print(' = %s'%(t_func_defaults[t_func_iarg]), end='')
        #    print(')', end='')
        # 
        # Directly print str(inspect.signature()) is the new way
        if print_to_screen == True:
            print(' ', end='')
            print(str(inspect_func_signature(t_func_self)), end='')
            print('')
        # 
        # import
        #print('globals()['+t_func_name+']')
        #globals()[t_func_name] = t_func_self
        #locals()[t_func_name] = t_func_self
        module = sys.modules[__name__]
        #print(module)
        if not hasattr(module, t_func_name):
            #print('setattr(module, '+t_func_name+', t_func_self)')
            setattr(module, t_func_name, t_func_self)





# 
# list all functions
# 
def help():
    print('Definitions:')
    print('    z        = galaxy\'s redshift')
    print('    SFR      = galaxy\'s star formation rate in units of solar mass per year')
    print('    Mstar    = galaxy\'s stellar mass in units of solar mass')
    print('    lgMstar  = log10(Mstar)')
    print('    sSFR     = SFR/Mstar*1e9 is a galaxy\' specific-SFR in units of Gyr^{-1}')
    print('    SFR_MS   = calc_SFR_MS(z, lgMstar) is a galaxy\'s main sequence SFR depending on its z and Mstar')
    print('    DeltaMS  = log10(SFR/SFR_MS) is a galaxy\'s offset to the galaxy main sequence')
    print('    deltaGas = Mmolgas/Mstar is a galaxy\'s molecular gas to stellar mass ratio')
    print('    muGas    = deltaGas')
    print('    fGas     = Mmolgas/(Mstar+Mmolgas) is a galaxy\'s molecular gas to total mass fraction')
    print('    tauGas   = Mmolgas/SFR*1e9 is a galaxy\'s molecular gas depletion time in units of Gyr')
    print('    tauDepl  = tauGas')
    print('    MetalZ   = 12+log10(O/H) is a galaxy\'s gas phase metallicity')
    print('    fmol     = Mmolgas/(Matomicgas+Mmolgas) is a galaxy\'s molecular fraction compared to total atomic+molecular gas')
    print('    GDR      = Mtotalgas/Mdust is a galaxy\'s total atomic+molecular gas to dust mass ratio')
    print('    deltaGDR = GDR')
    print('')
    print('Functions:')
    print('    # To get cosmic cold molecular gas evolution curve (D. Liu et al. 2019b Fig. 15): ')
    load_all_functions_in_one_module(get_cosmic_mol_gas_density, r'^get_cosmic_mol_gas_density_.*', print_to_screen = True) # print all function within this module
    print('')
    print('    # To calculate gas fraction (deltaGas or muGas) with z, Mstar and SFR using gas scaling relations (note that here is muGas not fGas): ')
    load_all_functions_in_one_module(calc_gas_fraction, r'^calc_gas_fraction_.*', print_to_screen = True) # print all function within this module
    print('')
    print('    # To calculate SFR with z and Mstar using galaxy main sequence correlations: ')
    load_all_functions_in_one_module(calc_galaxy_main_sequence, r'^calc_SFR_MS_.*', print_to_screen = True) # print all function within this module
    print('')
    print('    # To calculate gas depletion time (tauDepl or tauGas) with z, Mstar and SFR using gas scaling relations: ')
    load_all_functions_in_one_module(calc_gas_depletion_time, r'^calc_gas_depletion_time_.*', print_to_screen = True) # print all function within this module
    print('')
    
    print('Examples:')
    print('    import a3cosmos_gas_evolution as a3g')
    print('    a3g.calc_gas_fraction_A3COSMOS(z=3.0, lgMstar=10.5, DeltaMS=0.5)')
    print('    = %s'%(calc_gas_fraction_A3COSMOS(z=3.0, lgMstar=10.5, DeltaMS=0.5)))
    print('    a3g.calc_gas_depletion_time_A3COSMOS(z=3.0, lgMstar=10.5, DeltaMS=0.5)')
    print('    = %s'%(calc_gas_depletion_time_A3COSMOS(z=3.0, lgMstar=10.5, DeltaMS=0.5)))
    print('    # which means that a galaxy at redshift 3.0 with stellar mass 10^{10.5} solMass and 0.5 dex above the main sequence')
    print('    # has a molecular gas mass 3.84 times its stellar mass, and the molecular gas will be consumed possibly via ')
    print('    # star formation activities within only 0.32 billion years.')
    print('    ')
    
    
    #'calc_galaxy_main_sequence.py'
    #'calc_gas_depletion_time.py'
    #'calc_gas_fraction.py'
    #'calc_star_formation_law.py'
    #'calc_stellar_mass_function.py'
    #'calc_stellar_mass_function.py'
    #'calc_alpha_CO.py'
    #'calc_cosmic_star_formation_rate_density.py'
    #'calc_delta_GD.py'
    #'calc_fmol.py'
    #'calc_gas_mass_from_dust.py'




# 
# main function
# 
if __name__ == '__main__':
    help()


