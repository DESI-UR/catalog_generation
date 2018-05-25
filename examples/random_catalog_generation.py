import sys
sys.path.append("../py/")
sys.path.append("..")

from generation_tools import GeneralTools

import matplotlib.pyplot as plt

import numpy as np

from astropy.io import fits


# Generating Mock Catalogs
# This example shows the use of individual components of the mock generation tools. Each step is explained followed
# by the use of the specific function. Here we use a configuration file which is similar to the original "../data/catget.cfg".
# The only difference is the pathing for the input files in the configuration file. In the future, we will try to remove the hardcoded
# paths from the configuration file.

gt = GeneralTools("../data/catgen_nb.cfg")

# Generating flat galaxies
# The function (generate_galaxies) is used to generate the galaxies in the mock catalog. User can choose either to use a
# uniform distribution in r in [r_min, r_max] or use the distribution of r from the provided template. The example below show
# the latter use. For the former, one needs to provide the argument 'uniform=True' in the function.
# 

r_flat, theta_flat, phi_flat = gt.generate_flat_galaxies()
print("Number of center galaxies: {}".format(len(r_flat)))

# Before finalizing the mock, distances are converted to redshift.
# Then z-acceptance is applied to the generated galaxies
r2z_function = gt.generate_LUT_r2z()
# first, the center galaxies:
z_flat     = r2z_function(r_flat)
z_acceptance = z_flat<=gt.z_max 
z_flat     = z_flat[z_acceptance]
theta_flat = theta_flat[z_acceptance]
phi_flat   = phi_flat[z_acceptance]

all_zs     = z_flat
all_thetas = theta_flat
all_phis   = phi_flat
all_types  = np.full(len(z_flat), 4)

gt.write_to_fits(col1=all_phis, col2=all_thetas, col3=all_zs, col4=np.ones(len(all_zs)), col5=all_types, filename="random.fits", coordinates=1)
                 
