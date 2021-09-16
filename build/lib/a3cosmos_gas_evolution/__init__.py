name = 'a3cosmos_gas_evolution'





# 
# import functions
# 
from .Common_Python_Code import calc_galaxy_main_sequence
from .Common_Python_Code import calc_gas_fraction
from .Common_Python_Code import calc_gas_depletion_time
from .Common_Python_Code import calc_alpha_CO
from .Common_Python_Code import calc_delta_GD
from .Common_Python_Code import calc_metal_Z
from .Common_Python_Code import calc_fmol
from .Common_Python_Code import calc_gas_mass_from_dust
from .Common_Python_Code import apply_cosmology
cosmo = apply_cosmology.cosmo






from .__main__ import help

from .__main__ import load_all_modules

load_all_modules()

from .__main__ import *

