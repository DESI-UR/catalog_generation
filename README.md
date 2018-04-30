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

## Introduction

To fully simulate a signal due Baryon Acoustic Oscillations (BAO) in galaxy distribution one needs to run a full simulation of the Universe evolving from the era of decoupling to present day. It is very CPU consuming and thus limits the range of parameters that can be tested. At the same time future high statistics experiments will demand a control over the systematic uncertainties down to sub-percent level. Here, we suggest a model that would significantly reduce this computational time, while having main features present in the observation data. 

The mock and random catalog generation requires two distributions: z distribution, P(r), and a completeness map of an experiment, R$_{ang}$($\theta$, $\phi$).

We start by randomly generating N$_{center}$ new galaxies according to the product of R$_{ang}$($\theta$, $\phi$) and P(r). These galaxies form centers of BAO ”bubbles”. Then, for each BAO center we add N$_{rim}$ galaxies located at a distance r$_{BAO}$ from the center and uniformly distributed over solid angle around the BAO center. Added galaxies can be discarded to correct for the acceptance in r and angles ($\theta$, $\phi$).

Then, N$_{flat}$ random galaxies are generated according to the product of R$_{ang}$($\theta$, $\phi$) and P(r). After this, some ”clumping” is introduced around NR$_{clump}$ galaxies drawn randomly from the newly-created mock catalog to mimic galaxy clusters, with BAO centers having a different number of clumping partners because of the potential access of dark matter around the BAO centers. Clumping follows the following distribution:
f(r) = Ax(r$_0$/r)$^\gamma$
where r is the distance of the center of the cluster, A is a normalization parameter and r$_0$ and $\gamma$ are input parameters.

## Usage
The usage is shown in `mock_catalog_generation.py` and `random_catalog_generation.py`.

In order to import packages from the `py` directory, following lines are required in the python script:
```
import sys
sys.path.append("../py")
from generation_tools import GenerationTools
```

The examples shown here needs the object created with a configuration file. Example configuration file is at `../data/catgen.cfg`
```
gen_tool = GenerationTools(CONFIGURATION_FILE)
```

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
```
center_clumps, flat_clumps, flats = gt.generate_clumps(np.append(r_center, r_rim), np.append(theta_center, theta_rim), np.append(phi_center, phi_rim))
```

### Diagnostics
```

```