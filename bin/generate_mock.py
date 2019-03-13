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
                    help="Output filename (this overwrites the output filename in the configuration file",
                    required=False,
                    default=None)
parser.add_argument("-d", "--diagnostics",
                    type=bool,
                    help="Flag to get diagnostic messages",
                    required=False,
                    default=False)
parser.add_argument("-a", "--add_to_rims",
                    type=bool,
                    help="Flag to add clumping galaxies around rim galaxies",
                    required=False,
                    default=False)
parser.add_argument("-e", "--extended",
                    type=bool,
                    help="Flag to write the extended output for the fits file",
                    required=False,
                    default=False)
parser.add_argument("-p", "--pickle",
                    type=bool, help="Flag to write the extended pickle file",
                    required=False,
                    default=False)
parser.add_argument("-u", "--use_style",
                    type=bool, help="Flag to use the provided plotting style",
                    required=False,
                    default=True)
args   = parser.parse_args()

# Generating Mock Catalogs

# Instantiate the module with the user defined configuration file
# Detailed information about individual methids are given below as well
gt = generalTools(args.config, args.diagnostics, use_style=args.use_style)

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
gt.generate_clumps(args.add_to_rims)

# Save the output
# there are two options for user. One output file is in fits format to be used with other TPCF codes
# the other output is in pickle format and it stores detailed information about the mock catalog generated
gt.write_to_fits(filename=args.output, save_extended=args.extended)

# If the user want to store everything is the catalog structure (galaxies with their corresponding rims, clumps and all)
# the user can uncomment the following line. Be cautious that this requires considerable memory and it may fail if used on
# simple laptop or desktop computer
if args.pickle:
    if args.output is not None:
        gt.write_to_pickle((args.output).replace(".fits", ".pkl"))
    else:
        gt.write_to_pickle()
