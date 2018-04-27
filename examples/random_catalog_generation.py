
# coding: utf-8

# In[1]:


import sys
sys.path.append("../py/")
sys.path.append("..")

from generation_tools import GeneralTools

import matplotlib.pyplot as plt

import numpy as np

from astropy.io import fits


# #### Generating Mock Catalogs
# This example shows the use of individual components of the mock generation tools. Each step is explained followed by the use of the specific function.
# Here we use a configuration file which is similar to the original "../data/catget.cfg". The only difference is the pathing for the input files in the configuration file. In the future, we will try to remove the hardcoded paths from the configuration file.

gt = GeneralTools("../data/catgen_nb.cfg")

# ##### Generating central and flat galaxies
# The function (generate_galaxies) is used to generate the galaxies in the mock catalog. User can choose either to use a uniform distribution in r in [r_min, r_max] or use the distribution of r from the provided template. The example below show the latter use. For the former, one needs to provide the argument 'uniform=True' in the function.
# 
# The number of galaxies that will be generated is defined by $num\_obs$

for i in range(10):

    r, theta, phi = gt.generate_galaxies(num_obs=gt.n_center)

    # #### Generating the clumps
    # The clumps that mimick the clumping due to dark matter is introduced via $f(r)=A\left(r_0/r\right)^\gamma$ where A is the normalization factor, $r_0$ and $\gamma$ are the input parameters.
    # 
    # For generating the clumps, the user has to provide the center galaxy coordinates. These coordinates are used in the code. Also in the code, flat galaxies are generated and used for injecting clumps. All the input parameters are defined in the configuration file.

    center_clumps, flat_clumps, flats = gt.generate_clumps(r, theta, phi)
    
    all_rs     = np.append(r, np.append(center_clumps[0], np.append(flat_clumps[0], flats[0])))
    all_zs     = gt.r2z(r=all_rs)
    all_thetas = np.append(theta, np.append(center_clumps[1], np.append(flat_clumps[1], flats[1])))
    all_phis   = np.append(phi, np.append(center_clumps[2], np.append(flat_clumps[2], flats[2])))
    all_types  = np.append(np.full(len(r), 0),
                           np.append(np.full(len(center_clumps[0]), 2),
                                     np.append(np.full(len(flat_clumps[0]), 3), np.full(len(flats[0]), 4))))
    gt.write_to_fits(col1=all_phis, col2=all_thetas, col3=all_zs, col4=np.ones(len(all_rs)), col5=all_types,
                     filename="random{:02d}.fits".format(i), coordinates=1)

