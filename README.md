# Fast Mock Catalog Generation

[![Build Status](https://travis-ci.com/DESI-UR/catalog_generation.svg?branch=master)](https://travis-ci.com/DESI-UR/catalog_generation)
[![CodeFactor](https://www.codefactor.io/repository/github/desi-ur/catalog_generation/badge/master)](https://www.codefactor.io/repository/github/desi-ur/catalog_generation/overview/master)

<!-- toc -->
* [Requirements](#requirements)
* [Installtion](#installation)
* [Introduction](#introduction)
* [Usage](#usage)
  * [Using Executables](#using-executables)
  * [Using Modules](#using-modules)
<!-- tocstop -->

## Requirements
The version of the python packages listed below are the tested version. Any recent versions of the packages should technically work.
* [Python](https://www.python.org/) (version >=3.5)
* NumPy  (version 1.13.1)
* HealPy (version 1.11.0)
* configparser (version 3.5.0)
* AstroPy (version 1.2.1)
* SciPy (version 0.19.1)
* Matplotlib (version 2.0.0)

Generating the mock requires a completeness map and a redshift distribution. An example completeness map is provided at `examples/data/example_coverage.npz`. As for the redshift distribution, user can use a data challenge file which has the information for individual galaxies or use a fits file similar to the one provided at `examples/data/example_zcat.fits`. 

## Introduction

To fully simulate a signal due Baryon Acoustic Oscillations (BAO) in galaxy distribution one needs to run a full simulation of the Universe evolving from the era of decoupling to present day. It is very CPU consuming and thus limits the range of parameters that can be tested. At the same time future high statistics experiments will demand a control over the systematic uncertainties down to sub-percent level. Here, we suggest a model that would significantly reduce this computational time, while having main features present in the observation data. 

The mock and random catalog generation requires two distributions: z distribution, P(r), and a completeness map of an experiment, R<sub>ang</sub>(&theta;, &phi;).

We start by randomly generating N<sub>center</sub> new galaxies according to the product of R<sub>ang</sub>(&theta;, &phi;) and P(r). These galaxies form centers of BAO ”bubbles”. Then, for each BAO center we add N<sub>rim</sub> galaxies located at a distance r<sub>BAO</sub> from the center and uniformly distributed over solid angle around the BAO center. Added galaxies can be discarded to correct for the acceptance in r and angles (&theta;, &phi;).

Then, N<sub>flat</sub> random galaxies are generated according to the product of R<sub>ang</sub>(&theta;, &phi;) and P(r). After this, some ”clumping” is introduced around NR<sub>clump</sub> galaxies drawn randomly from the newly-created mock catalog to mimic galaxy clusters, with BAO centers having a different number of clumping partners because of the potential access of dark matter around the BAO centers. Clumping follows the following distribution:
f(r) = Ax(r<sub>0</sub>/r)&gamma;
where r is the distance of the center of the cluster, A is a normalization parameter and r<sub>0</sub> and &gamma; are input parameters.

## Usage

### Using executables

There are two executables provided within the package. The insight of the executables is explain in the following section. The arguments for the executable are shown below:
```
usage: generate_mock.py [-h] -c CONFIG [-o OUTPUT] [-d DIAGNOSTICS]
                        [-a ADD_TO_RIMS] [-e EXTENDED] [-p PICKLE] [-u STYLE]

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
                        Configuration filename
  -o OUTPUT, --output OUTPUT
                        Output filename (this overwrites the output filename
                        in the configuration file
  -d DIAGNOSTICS, --diagnostics DIAGNOSTICS
                        Flag to get diagnostic messages
  -a ADD_TO_RIMS, --add_to_rims ADD_TO_RIMS
                        Flag to add clumping galaxies around rim galaxies
  -e EXTENDED, --extended EXTENDED
                        Flag to write the extended output for the fits file
  -p PICKLE, --pickle PICKLE
                        Flag to write the extended pickle file
  -u STYLE, --style STYLE
                        Flag to use the provided plotting style
```
`--config` argument is required for the configuration to be used.
Other arguments are optional and the explanation is provided with `--help`.

### Using modules

The usage is shown in `bin/generate_mock.py` and `bin/generate_random.py`. However, the user can write own executables or use them as template for python notebook use.

In order to import the module, following lines are required in the python script:
```
from paramock.paramock import GeneralTools as generalTools
```

The examples shown here needs the object created with a configuration file. Example configuration file is at `data/catgen.cfg`
```
gen_tool = generationTools(CONFIGURATION_FILE, diagnostics=DIAGNOSTICS)
```

In the configuration file, user needs to provide either a sample catalog with redshifts or a redshift distribution as `datafile`. The mock and random catalog to be generated within a window given in `angmask`. For the users who do not have such a file or want to generate a catalog in a custom region in the sky, a python script, `examples/generate_completeness_map.py`, is included in the project. This script generates a completeness map for a region of interest. For more complicated completeness maps, users need to write their own scripts considering the data structure shown in `examples/generate_completeness_map.py`

### Generating center galaxies
One needs to generate the center galaxies to be used for introducing the BAO galaxies. The number of center galaxies are defined in the configuration file. The example line below shows how to generate center galaxies.
```
gen_tool.generate_center_galaxies()
```

### Generating the rim galaxies
The BAO signal is introduced around the center galaxies with a gaussian distribution (&mu;, &sigma;)=(r<sub>BAO</sub>, sigma<sub>r_BAO</sub>). To generate them, one needs to use the `generate_rim` function. The position of the center galaxies are pulled from the catalog object within `gen_tool`.
```
gen_tool.generate_rim()
```

### Generating the clump and flat galaxies
In addition to the center and rim galaxies, randomy distributed flat galaxies needed to be introduced. It is automatically done when the clumps are being generated. Then, a number of clump galaxies are injected centered around some selected galaxies. The clumps are categorized into two groups: center seeded clumps and flat seeded clumps. The defition of center seeded clumps is left open for the user. User can use rims for the center seed clumps. The second group, flat seeded clumps, are selected among the flat galaxies generated. In the example below, rim galaxies are added to the center seeded clumps by setting `add_clumps_to_rims=True` (which is False by default).
```
gen_tool.generate_clumps(add_clumps_to_rims=True)
```
There is no difference in generating both group of clumps, they both use the same power law distribution for their location with respect to the seed galaxy of a clump. The only difference is the number of clumps around a certain type of object. This is controlled by the two parameters in the configuration file: `frac_f2c` and `frac_c2r`. `frac_f2c` defines the fraction of flat seeded clumps to center seeded clumps and `frac_c2r` defines the fraction of center galaxies used with respect to the rim galaxies for the center seed galaxies. If the parameters are not defined in the configuration, `gen_tool` defaults to `None` and use random distributions for the selections

### Saving the output
The output can be saved in two different formats, fits and pickle. `.fits` output is structured to run analysis easily and the fields stored can be extended to have the name and the parent name of the galaxies with the flag `save_extended`. The following line shows how to save the output in fits format.
```
gen_tool.write_to_fits(FILENAME, save_extended=True)
```
If the filename is not provided when calling the function, the module will use the filename defined in the configuration file. The second output type, `.pkl`, stores the whole catalog as an object so that the user can use the functions in the catalog and galaxy modules easily for postprocessing.
```
gen_tool.write_to_pickle(FILENAME)
```
Similarly, when the filename is not provided, the filename in the configuration file will be used after changing the extension to `.pkl`.



