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

# Generating central and flat galaxies
# The function (generate_galaxies) is used to generate the galaxies in the mock catalog. User can choose either to use a
# uniform distribution in r in [r_min, r_max] or use the distribution of r from the provided template. The example below show
# the latter use. For the former, one needs to provide the argument 'uniform=True' in the function.
# 
# The number of galaxies that will be generated is defined by $num\_obs$

r_center, theta_center, phi_center = gt.generate_center_galaxies()
print("Number of center galaxies: {}".format(len(r_center)))

# Generating the rim galaxies
# The rim galaxies are added to the mock using a gaussian distribution centered at the galaxies. The acceptance is also
# applied to the generated rim galaxies using the completeness provided.
# The number of rim galaxies is read from the configuration file. Later, we will add the feature to define the number of
# rim galaxies through the function keeping the default to use the configuration file.

r_rim, theta_rim, phi_rim = gt.generate_rim(r_center, theta_center, phi_center)
print("Number of rim galaxies: {}".format(len(r_rim)))

# Generating the clumps
# The clumps that mimick the clumping due to dark matter is introduced via $f(r)=A\left(r_0/r\right)^\gamma$ where A
# is the normalization factor, $r_0$ and $\gamma$ are the input parameters.
# For generating the clumps, the user has to provide the center galaxy coordinates. These coordinates are
# used in the code. Also in the code, flat galaxies are generated and used for injecting clumps. All the input parameters
# are defined in the configuration file.

center_clumps, flat_clumps, flats = gt.generate_clumps(np.append(r_center, r_rim), np.append(theta_center, theta_rim), np.append(phi_center, phi_rim))
print("Number of flat galaxies: {}".format(len(flats[0])))
print("Number of center clump galaxies: {}".format(len(center_clumps[0])))
print("Number of flat clump galaxies: {}".format(len(flat_clumps[0])))

# Final Stage
# Before finalizing the mock, distances are converted to redshift.
# Then z-acceptance is applied to the generated galaxies
r2z_function = gt.generate_LUT_r2z()
# first, the center galaxies:
z_center     = r2z_function(r_center)
z_acceptance = z_center<=gt.z_hi #gt.check_z_acceptance(z_center)
z_center     = z_center[z_acceptance]
theta_center = theta_center[z_acceptance]
phi_center   = phi_center[z_acceptance]
# second, the rim galaxies
z_rim        = r2z_function(r_rim)
z_acceptance = gt.check_z_acceptance(z_rim)
z_rim        = z_rim[z_acceptance]
theta_rim = theta_rim[z_acceptance]
phi_rim   = phi_rim[z_acceptance]
# third, the flat galaxies
z_flat       = r2z_function(flats[0])
z_acceptance = z_flat<=gt.z_hi#gt.check_z_acceptance(z_flat)
z_flat       = z_flat[z_acceptance]
theta_flat   = flats[1][z_acceptance]
phi_flat     = flats[2][z_acceptance]
# fourth, the center clumps
z_center_clumps     = r2z_function(center_clumps[0])
z_acceptance        = gt.check_z_acceptance(z_center_clumps)
z_center_clumps     = z_center_clumps[z_acceptance]
theta_center_clumps = center_clumps[1][z_acceptance]
phi_center_clumps   = center_clumps[2][z_acceptance]
# fifth, the flat clumps
z_flat_clumps     = r2z_function(flat_clumps[0])
z_acceptance      = gt.check_z_acceptance(z_flat_clumps)
z_flat_clumps     = z_flat_clumps[z_acceptance]
theta_flat_clumps = flat_clumps[1][z_acceptance]
phi_flat_clumps   = flat_clumps[2][z_acceptance]
# lastly, we duplicate the z-template and apply acceptance test for diagnostic purposes
template_z   = np.array(gt.template_z)
template_z   = template_z[template_z<gt.z_hi]

"""
# The results
# The distributions and the positions of the galaxies are plotted here. As can be seen, the acceptance is applied during the generation.
plt.figure(1, figsize=(16,8))

plt.subplot(241)
plt.hist(template_z, bins=int(np.sqrt(len(template_z))), alpha=0.3, normed=True, label="Template")
plt.hist(z_center, bins=int(np.sqrt(len(z_center))), alpha=0.3, normed=True, label="BAO Centers")
plt.legend(loc=2)
plt.xlabel("z [redshift]")

plt.subplot(242)
plt.hist(template_z, bins=int(np.sqrt(len(template_z))), alpha=0.3, normed=True, label="Template")
plt.hist(z_rim, bins=int(np.sqrt(len(z_rim))), alpha=0.3, normed=True, label="BAO rims")
plt.legend(loc=2)
plt.xlabel("z [redshift]")

plt.subplot(243)
plt.hist(template_z, bins=int(np.sqrt(len(template_z))), alpha=0.3, normed=True, label="Template")
plt.hist(z_center_clumps, bins=int(np.sqrt(len(z_center_clumps))), alpha=0.3, normed=True, label="Center clumps")
plt.legend(loc=2)
plt.xlabel("z [redshift]")

plt.subplot(244)
plt.hist(template_z, bins=int(np.sqrt(len(template_z))), alpha=0.3, normed=True, label="Template")
plt.hist(z_flat_clumps, bins=int(np.sqrt(len(z_flat_clumps))), alpha=0.3, normed=True, label="Flat clumps")
plt.legend(loc=2)
plt.xlabel("z [redshift]")

plt.subplot(245)
plt.hist(template_z, bins=int(np.sqrt(len(template_z))), alpha=0.3, normed=True, label="Template")
plt.hist(z_flat, bins=int(np.sqrt(len(z_flat))), alpha=0.3, normed=True, label="Flats")
plt.legend(loc=2)
plt.xlabel("z [redshift]")

ax = plt.subplot(246)
ax.hist(theta_center, alpha=0.3, normed=True)
ax.hist(theta_rim, alpha=0.3, normed=True)
ax.hist(theta_center_clumps, alpha=0.3, normed=True)
ax.hist(theta_flat_clumps, alpha=0.3, normed=True)
ax.hist(theta_flat, alpha=0.3, normed=True)
ax.set_xlabel(r'$\theta$ [deg]')

ax = plt.subplot(247)
ax.hist(phi_center, alpha=0.7, label="center galaxies", normed=True)
ax.hist(phi_rim, alpha=0.3, label="rim galaxies", normed=True)
ax.hist(phi_center_clumps, alpha=0.3, label="center clumps", normed=True)
ax.hist(phi_flat_clumps, alpha=0.3, label="flat clumps", normed=True)
ax.hist(phi_flat, alpha=0.3, label="flat galaxies", normed=True)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlabel(r'$\phi$ [deg]')
plt.savefig("distributions.pdf", bbox_inches='tight')
plt.show()
"""
"""
plt.figure(2, figsize=(12,12))
plt.subplot(221)
plt.plot(theta_center, phi_center, 'o', alpha=0.2, label="center galaxies")
plt.plot(theta_rim, phi_rim, 'o', alpha=0.2, label="rim galaxies")
plt.xlim(0, 360)
plt.ylim(-90, 90)
plt.xlabel(r"$\theta$ [deg]")
plt.ylabel(r"$\phi$ [deg]")
plt.legend()
plt.grid(True)
plt.subplot(222)
plt.plot(theta_center, phi_center, 'o', alpha=0.2, label="center galaxies")
plt.plot(center_clumps[1], center_clumps[2], 'o', alpha=0.2, label="central clumps")
plt.xlim(0, 360)
plt.ylim(-90, 90)
plt.xlabel(r"$\theta$ [deg]")
plt.ylabel(r"$\phi$ [deg]")
plt.legend()
plt.grid(True)
plt.subplot(223)
plt.plot(theta_center, phi_center, 'o', alpha=0.2, label="center galaxies")
plt.plot(flat_clumps[1], flat_clumps[2], 'o', alpha=0.2, label="flat clumps")
plt.xlim(0, 360)
plt.ylim(-90, 90)
plt.xlabel(r"$\theta$ [deg]")
plt.ylabel(r"$\phi$ [deg]")
plt.legend()
plt.grid(True)
plt.subplot(224)
plt.plot(theta_center, phi_center, 'o', alpha=0.2, label="center galaxies")
plt.plot(flats[1], flats[2], 'o', alpha=0.2, label="flat galaxies")
plt.xlim(0, 360)
plt.ylim(-90, 90)
plt.xlabel(r"$\theta$ [deg]")
plt.ylabel(r"$\phi$ [deg]")
plt.legend()
plt.grid(True)
plt.savefig("maps.pdf")
plt.show()
"""

all_zs     = np.append(z_center, np.append(z_rim, np.append(z_center_clumps, np.append(z_flat_clumps, z_flat))))
all_thetas = np.append(theta_center, np.append(theta_rim, np.append(theta_center_clumps, np.append(theta_flat_clumps, theta_flat))))
all_phis   = np.append(phi_center, np.append(phi_rim, np.append(phi_center_clumps, np.append(phi_flat_clumps, phi_flat))))
all_types  = np.append(np.full(len(z_center), 0),
                       np.append(np.full(len(z_rim), 1),
                                 np.append(np.full(len(z_center_clumps), 2),
                                           np.append(np.full(len(z_flat_clumps), 3), np.full(len(z_flat), 4)))))
gt.write_to_fits(col1=all_phis, col2=all_thetas, col3=all_zs, col4=np.ones(len(all_zs)), col5=all_types,
                 filename="mock.fits", coordinates=1)
