
# coding: utf-8

# In[1]:


import sys
sys.path.append("../py/")
sys.path.append("..")

from generation_tools import GeneralTools

import matplotlib.pyplot as plt

import numpy as np

from astropy.io import fits

# Generating Mock Catalogs
# This example shows the use of individual components of the mock generation tools. Each step is explained
# followed by the use of the specific function. Here we use a configuration file which is similar to the
# original "../data/catget.cfg". The only difference is the pathing for the input files in the configuration file.
# In the future, we will try to remove the hardcoded paths from the configuration file.

gt = GeneralTools("../data/catgen_nb.cfg")

for i in range(10):

    # Generating central and flat galaxies
    # The function (generate_galaxies) is used to generate the galaxies in the mock catalog. User can choose either
    # to use a uniform distribution in r in [r_min, r_max] or use the distribution of r from the provided template.
    # The example below show the latter use. For the former, one needs to provide the argument 'uniform=True' in the function.
    # The number of galaxies that will be generated is defined by $num\_obs$
    r_center, theta_center, phi_center = gt.generate_galaxies(num_obs=gt.n_center)
        
    # Generating the clumps
    # The clumps that mimick the clumping due to dark matter is introduced via $f(r)=A\left(r_0/r\right)^\gamma$
    # where A is the normalization factor, $r_0$ and $\gamma$ are the input parameters.
    # For generating the clumps, the user has to provide the center galaxy coordinates. These coordinates are
    # used in the code. Also in the code, flat galaxies are generated and used for injecting clumps. All the
    # input parameters are defined in the configuration file.
    center_clumps, flat_clumps, flats = gt.generate_clumps(r_center, theta_center, phi_center)
    
    # Final Stage
    # Before finalizing the mock, distances are converted to redshift.
    # Then z-acceptance is applied to the generated galaxies
    r2z_function = gt.generate_LUT_r2z()
    # first, the center galaxies:
    z_center     = r2z_function(r_center)
    z_acceptance = gt.check_z_acceptance(z_center)
    z_center     = z_center[z_acceptance]
    theta_center = theta_center[z_acceptance]
    phi_center   = phi_center[z_acceptance]
    # second, the flat galaxies
    z_flat       = r2z_function(flats[0])
    z_acceptance = gt.check_z_acceptance(z_flat)
    z_flat       = z_flat[z_acceptance]
    theta_flat   = flats[1][z_acceptance]
    phi_flat     = flats[2][z_acceptance]
    # third, the center clumps
    z_center_clumps     = r2z_function(center_clumps[0])
    z_acceptance        = gt.check_z_acceptance(z_center_clumps)
    z_center_clumps     = z_center_clumps[z_acceptance]
    theta_center_clumps = center_clumps[1][z_acceptance]
    phi_center_clumps   = center_clumps[2][z_acceptance]
    # fourth, the flat clumps
    z_flat_clumps     = r2z_function(flat_clumps[0])
    z_acceptance      = gt.check_z_acceptance(z_flat_clumps)
    z_flat_clumps     = z_flat_clumps[z_acceptance]
    theta_flat_clumps = flat_clumps[1][z_acceptance]
    phi_flat_clumps   = flat_clumps[2][z_acceptance]

    #
    all_zs     = np.append(z_center, np.append(z_center_clumps, np.append(z_flat_clumps, z_flat)))
    all_thetas = np.append(theta_center, np.append(theta_center_clumps, np.append(theta_flat_clumps, theta_flat)))
    all_phis   = np.append(phi_center, np.append(phi_center_clumps, np.append(phi_flat_clumps, phi_flat)))
    all_types  = np.append(np.full(len(z_center), 0),
                           np.append(np.full(len(z_center_clumps), 2),
                                     np.append(np.full(len(z_flat_clumps), 3), np.full(len(z_flats), 4))))
    gt.write_to_fits(col1=all_phis, col2=all_thetas, col3=all_zs, col4=np.ones(len(all_rs)), col5=all_types,
                     filename="random{:02d}.fits".format(i), coordinates=1)
