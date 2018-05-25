# Fast Mock Catalog Generation

<!-- toc -->
- [Requirements](#requirements)
- [Introduction](#introduction)
- [Usage](#usage)
<!-- tocstop -->

## Requirements
The version of the python packages listed below are the tested version. Any recent versions of the packages should technically work.
* [Python](https://www.python.org/) (version >2.7)
* NumPy  (version 1.13.1)
* HealPy (version 1.11.0)
* configparser (version 3.5.0)
* AstroPy (version 1.2.1)
* SciPy (version 0.19.1)
* Matplotlib (version 2.0.0)

Generating the mock requires a completeness map and a redshift distribution. An example completeness map is provided at `data/completeness.tar.gz` that user needs to extract. As for the redshift distribution, user can use a data challenge file which has the information for individual galaxies or use a fits file similar to the one provided at `data/example.fits`. 

## Introduction

To fully simulate a signal due Baryon Acoustic Oscillations (BAO) in galaxy distribution one needs to run a full simulation of the Universe evolving from the era of decoupling to present day. It is very CPU consuming and thus limits the range of parameters that can be tested. At the same time future high statistics experiments will demand a control over the systematic uncertainties down to sub-percent level. Here, we suggest a model that would significantly reduce this computational time, while having main features present in the observation data. 

The mock and random catalog generation requires two distributions: z distribution, P(r), and a completeness map of an experiment, R<sub>ang</sub>(&theta;, &phi;).

We start by randomly generating N<sub>center</sub> new galaxies according to the product of R<sub>ang</sub>(&theta;, &phi;) and P(r). These galaxies form centers of BAO ”bubbles”. Then, for each BAO center we add N$_{rim}$ galaxies located at a distance r<sub>BAO</sub> from the center and uniformly distributed over solid angle around the BAO center. Added galaxies can be discarded to correct for the acceptance in r and angles (&theta;, &phi;).

Then, N$_{flat}$ random galaxies are generated according to the product of R<sub>ang</sub>(&theta;, &phi;) and P(r). After this, some ”clumping” is introduced around NR$_{clump}$ galaxies drawn randomly from the newly-created mock catalog to mimic galaxy clusters, with BAO centers having a different number of clumping partners because of the potential access of dark matter around the BAO centers. Clumping follows the following distribution:
f(r) = Ax(r<sub>0</sub>/r)&gamma;
where r is the distance of the center of the cluster, A is a normalization parameter and r<sub>0</sub> and &gamma; are input parameters.

## Usage
The usage is shown in `mock_catalog_generation.py` and `random_catalog_generation.py`.

In order to import packages from the `py` directory, following lines are required in the python script:
```
import sys
sys.path.append("../py")
from generation_tools import GenerationTools
```

The examples shown here needs the object created with a configuration file. Example configuration file is at `data/catgen.cfg`
```
gen_tool = GenerationTools(CONFIGURATION_FILE)
```

In the configuration file, user needs to provide either a sample catalog with redshifts or a redshift distribution as `datafile`. The mock and random catalog to be generated within a window given in `angmask`. For the users who do not have such a file or want to generate a catalog in a custom region in the sky, a python script, `examples/generate_completeness_map.py`, is included in the project. This script generates a completeness map for a region of interest. For more complicated completeness maps, users need to write their own scripts considering the data structure shown in `examples/generate_completeness_map.py`

### Generating center galaxies
One needs to generate the center galaxies to be used for introducing the BAO galaxies. The number of center galaxies can be either entered a number or the number defined in the configuration file. The example line below shows how to generate center galaxies using the configuration file. 
```
r_center, theta_center, phi_center = gen_tool(num_obs=gen_tool.n_center)
```

### Generating the rim galaxies
The BAO signal is introduced around the center galaxies with a gaussian distribution ($\mu$, $\sigma$)=(r_BAO, sigma_r_BAO). To generate them, one needs to use the `generate_rim` function with the position of the center galaxies:
```
r_rim, theta_rim, phi_rim = gen_tool.generate_rim(r_center, theta_center, phi_center)
```

### Generating the clump and flat galaxies
In addition to the center and rim galaxies, randomy distributed flat galaxies are introduced as shown in the example below. Then, a number of clump galaxies are injected centered around the all the galaxies. The clumps are categorized into two groups: center clumps and flat clumps. The defition of center clumps is left open for the user. User can use rims for the center clumps or user can leave them for the second group, flat clumps. In the example below, rim galaxies are considered as centers for center clumps. There is no difference in generating both group of clumps, they both use the same power law distribution for their location with respect to the center of a clump. The only difference is the number of clumps around a certain type of object.
```
center_clumps, flat_clumps, flats = gt.generate_clumps(np.append(r_center, r_rim), np.append(theta_center, theta_rim), np.append(phi_center, phi_rim))
```

### Final steps
The galaxies are all generated using a completeness map so that user will only get galaxies in the field of view of an experiment. However, no redshift acceptance is applied until now. The example here shows how a user can apply the redshift acceptance test:
```
r2z_function = gt.generate_LUT_r2z()
z            = r2z_function(r)
z_acceptance = gt.check_z_acceptance(z)
z            = z[z_acceptance]
theta        = theta[z_acceptance]
phi          = phi[z_acceptance]
```
The first line is an interpolator for z-to-r conversion which is faster than using the astropy package for all the values.