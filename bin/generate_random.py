#!/usr/bin/env python

import argparse

from paramock.paramock import GeneralTools as generalTools

parser = argparse.ArgumentParser(description="")
parser.add_argument("-c", "--config",
                    type=str,
                    help="Configuration filename",
                    required=True)
parser.add_argument("-o", "--output",
                    type=str,
                    help="Output filename (this overwrites the output filename in the configuration file)",
                    required=False,
                    default=None)
args   = parser.parse_args()

# Generating Mock Catalogs

# Instantiate the module with the user defined configuration file
# Detailed information about individual methids are given below as well
gt = generalTools(args.config)

# Generating flat galaxies
# The function (generate_galaxies) is used to generate the galaxies in the mock catalog. User can choose either to use a
# uniform distribution in r in [r_min, r_max] or use the distribution of r from the provided template. The example below show
# the latter use. For the former, one needs to provide the argument 'uniform=True' in the function.
r_flat, theta_flat, phi_flat = gt.generate_flat_galaxies(is_random=True)

# Save the output
# there are two options for user. One output file is in fits format to be used with other TPCF codes
# the other output is in pickle format and it stores detailed information about the mock catalog generated
gt.write_to_fits(filename=args.output, is_random=True)
