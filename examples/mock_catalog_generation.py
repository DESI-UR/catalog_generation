import sys
sys.path.append("../py")

import pickle

from generation_tools import GeneralTools

# Generating Mock Catalogs
# This example shows the use of individual components of the mock generation tools. Each step is explained followed
# by the use of the specific function. Here we use a configuration file which is similar to the original "../data/catget.cfg".
# The only difference is the pathing for the input files in the configuration file. In the future, we will try to remove the hardcoded
# paths from the configuration file.
gt = GeneralTools(sys.argv[1])

# Generating central galaxies
# The function (generate_galaxies) is used to generate the galaxies in the mock catalog. User can choose either to use a
# uniform distribution in r in [r_min, r_max] or use the distribution of r from the provided template. The example below show
# the latter use. For the former, one needs to provide the argument 'uniform=True' in the function.
gt.generate_center_galaxies()

# Generating the rim galaxies
# The rim galaxies are added to the mock using a gaussian distribution centered at the galaxies. The acceptance is also
# applied to the generated rim galaxies using the completeness provided.
# The number of rim galaxies is read from the configuration file. Later, we will add the feature to define the number of
# rim galaxies through the function keeping the default to use the configuration file.
gt.generate_rim()

# Generating the clumps
# The clumps that mimick the clumping due to dark matter is introduced via $f(r)=A\left(r_0/r\right)^\gamma$ where A
# is the normalization factor, $r_0$ and $\gamma$ are the input parameters.
# For generating the clumps, the user has to provide the center galaxy coordinates. These coordinates are
# used in the code. Also in the code, flat galaxies are generated and used for injecting clumps. All the input parameters
# are defined in the configuration file.
gt.generate_clumps(add_clumps_to_rims=True)

# This is an example of how one can plot a center galaxy with its rims and clumps
# Another example about this plotting option is on a python notebook
try:
    gt.catalog.plot("cen_0")
except:
    skip_plot = True

# Save the output
# there are two options for user. One output file is in fits format to be used with other TPCF codes
# the other output is in pickle format and it stores detailed information about the mock catalog generated
gt.write_to_fits()

# If the user want to store everything is the catalog structure (galaxies with their corresponding rims, clumps and all)
# the user can uncomment the following line. Be cautious that this requires considerable memory and it may fail if used on
# simple laptop or desktop computer
#gt.write_to_pickle()
