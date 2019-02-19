"""
Modules to generate random and mock catalogs

Contributors: Tolga Yapici [tyapici@ur.rochester.edu]
"""

import numpy as np
import healpy as hp
import math, random
import matplotlib as mpl
import matplotlib.pyplot as plt
import configparser
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u
import sys
import os
import pickle
from scipy.interpolate import interp1d

from galaxy import galaxy
from catalog import catalog

import time

NSIDE = -1
DEG2RAD = np.pi/180.0
RAD2DEG = 1./DEG2RAD

class GeneralTools():
    """
    Initialize the tools using a configuration file

    Parameters
    -----------
    configFile : string
        The configuration filename with the parameters to use

    Returns
    -------
    None
    """
    def __init__(self, configFile, diagnostics=False, acceptance=True):
        self.diagnostics = diagnostics
        self.acceptance  = acceptance
        self.config_file = configFile
        self.catalog     = catalog()
        self.getConfig()
        self.get_template()
        self.get_mask()
        
    """
    Function to change the diagnostics level. Even though it can be set during initialization, this
    function makes it possible to change the diagnostics during operation

    Parameters
    ----------
    diagnostics : boolean

    Returns
    -------
    None

    """
    def setDiagnostics(self, diagnostics):
        self.diagnostics = diagnostics
        

    """
    Function to read and process a line in configuration. 
    This is needed because the user may input parameters with dimensions.
    For example, the user can either define the BAO radius in Mpc or in h^-1*Mpc
    
    Parameters
    ----------
    config_line: string
         a line read by configparser to be processed for unit and dimension

    Returns
    -------
    quantity: astropy.units.Quantity
         value of the parameter in the configparser line. If the unit is given, it will be taken into account
    """
    def read_config_line(self, config_line):
        try:
            quantity = u.Quantity(config_line)
        except Exception as e:
            raise Exception('the value provided for the following line is invalid'
                            '{}'
                            'Please check for correct syntax, either a real number or a quantity with unit'.format(config_line))
        return quantity
    
    """
    Function to set the parameters of the object using the configuration file.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    def getConfig(self):
        self.config     = configparser.RawConfigParser()
        self.config.read(self.config_file)

        # some of the numbers to use in the generation
        self.num_random = int(self.config.get("N Catalogs", "nrandom"))
        self.num_mock   = int(self.config.get("N Catalogs", "nmock"))
        
        # read in the cosmology parameters to generate the mocks
        H0             = self.read_config_line(self.config.get("Cosmology", "H_0"))
        if H0.unit == u.dimensionless_unscaled:
            # if the unit is not given for H0, then it should be h0.
            # that needs to be changed to H0
            self.h0    = H0
            self.H0    = H0 * 100.
        else:
            # if the unit is given for H0, then it should be H0.
            self.h0    = H0.to("Mpc").value/100.
            self.H0    = H0.to("Mpc").value
        self.omega_m0  = float(self.config.get("Cosmology", "omega_m0"))
        self.omega_b0  = float(self.config.get("Cosmology", "omega_b0"))
        self.omega_de0 = float(self.config.get("Cosmology", "omega_de0"))
        self.temp_cmb  = float(self.config.get("Cosmology", "T_cmb"))
        self.nu_eff    = float(self.config.get("Cosmology", "nu_eff"))
        self.m_nu      = np.array(self.config.get("Cosmology", "m_nu").split(","), dtype=float)
        self.cosmo     = FlatLambdaCDM(name="Custom Cosmology", H0=self.H0, Om0=self.omega_m0, Ob0=self.omega_b0,
                                       Tcmb0=self.temp_cmb, Neff=self.nu_eff, m_nu=self.m_nu*u.eV)
        
        # output filenames
        self.fname_random = self.config.get("File Names", "randout")
        self.fname_mock   = self.config.get("File Names", "mockout")
        self.datafile     = self.config.get('File Names','datafile')
        self.ang_mask     = self.config.get('File Names','angmask')
        
        # the following variable defined the coordinate system, 0: theta, phi, 1: RA, DEC
        self.coordinates  = int(self.config.get("Data Format", "coordinates"))

        # BAO related parameters
        r_BAO               = self.read_config_line(self.config.get('Gen Params','r_bao'))
        if r_BAO.unit == u.dimensionless_unscaled:
            # we need to change r_BAO so that it will have the unit of Mpc
            self.r_BAO      = r_BAO/self.h0
        else:
            self.r_BAO      = r_BAO.to("Mpc").value
        sigma_r_BAO         = self.read_config_line(self.config.get('Gen Params','sigr_bao'))
        if sigma_r_BAO.unit == u.dimensionless_unscaled:
            self.sigma_r_BAO= sigma_r_BAO/self.h0
        else:
            self.sigma_r_BAO= sigma_r_BAO.to("Mpc").value

        # Clumping parameters
        self.gamma          = float(self.config.get('Gen Params','gamma'))
        r_0                 = self.read_config_line(self.config.get('Gen Params','r_0'))
        if r_0.unit == u.dimensionless_unscaled:
            self.r_0        = r_0/self.h0
        else:
            self.r_0        = r_0.to("Mpc").value
        self.n_rand         = int(self.config.get('Gen Params','n_rand'))
        self.n_center       = int(self.config.get('Gen Params','n_center'))
        self.n_rim          = int(self.config.get('Gen Params','n_rim'))
        self.n_flat         = int(self.config.get('Gen Params','n_flat'))
        self.nr_clump       = int(self.config.get('Gen Params','nr_clump'))
        try:
            self.clump_dist = int(self.config.get('Gen Params','clump_dist'))
        except Exception as e:
            self.clump_dist = 0
        try:
            self.frac_f2c   = int(self.config.get('Gen Params','frac_f2c'))
        except Exception as e:
            self.frac_f2c   = None
        try:
            self.frac_c2r   = int(self.config.get('Gen Params','frac_c2r'))
        except Exception as e:
            self.frac_c2r   = None
        try:
            self.a          = float(self.config.get('Gen Params','a'))
            self.b          = float(self.config.get('Gen Params','b'))
        except Exception as e:
            self.a          = 1.
            self.b          = 1.
        self.n_clump        = int(self.config.get('Gen Params','n_clump'))
        self.n_clump_center = int(self.config.get('Gen Params','n_centerclump'))
        try:
            self.z_min      = float(self.config.get("Gen Params", "z_min"))
            self.z_max      = float(self.config.get("Gen Params", "z_max"))
        except Exception as e:
            print("using predefined z bounds [0.4, 0.7]")
            self.z_min      = 0.4
            self.z_max      = 0.7
        self.r_max      = self.cosmo.comoving_distance(self.z_max).value
        self.r_min      = self.cosmo.comoving_distance(self.z_min).value
        
        self.anisotropy_scale = np.array([self.a, self.b, self.a])

        print("anisotropy scale is : {}".format(self.anisotropy_scale))

    """
    Function to generate a lookup table for z to r conversion
    
    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    def generate_LUT_z2r(self):
        z_min = 0.1
        z_max = 3.0
        zs = np.linspace(z_min, z_max, 1000)
        rs = [r.value for r in self.cosmo.comoving_distance(zs)]
        interpolator = interp1d(zs, rs, bounds_error=False, fill_value=-1.)
        return interpolator

    """
    Function to generate a lookup table for r to z conversion
    
    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    def generate_LUT_r2z(self):
        z_min = 0.1
        z_max = 3.0
        zs = np.linspace(z_min, z_max, 1000)
        rs = [r.value for r in self.cosmo.comoving_distance(zs)]
        interpolator = interp1d(rs, zs, bounds_error=False, fill_value=-9999.)
        return interpolator            

    """
    Function to read the template r distribution for mock and random generation

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    def get_template(self):
        self.hdulist      = fits.open(self.datafile)
        self.data_z       = self.hdulist[1].data['Z']
        #template_cut = np.array([(self.data_z[i] < self.z_max+0.125) and (self.data_z[i] > self.z_min) for i in range(len(self.data_z))])
        template_cut = np.logical_and(self.data_z<(self.z_max+0.125), self.data_z>self.z_min)
        self.template_z   = self.data_z[template_cut]
        try:
            self.template_w = self.hdulist[1].data['HistoW'][template_cut]
        except Exception as e:
            self.template_w = np.full(len(self.template_z), fill_value=1.)
        self.template_r   = [r.value for r in self.cosmo.comoving_distance(self.template_z)]
        self.template_r_len = len(self.template_r)
        self.template_r_min = np.amin(self.template_r)
        self.template_r_max = np.amax(self.template_r)
                    
    """
    Function to read the completeness map for mock and random generation

    Parameters
    ----------
    diagnostics : boolean

    Returns
    -------
    None
    """
    def get_mask(self, diagnostics=False):
        self.mask             = np.load(self.ang_mask)
        self.completeness     = self.mask.f.CMPLTNSS
        self.completeness_len = len(self.completeness)
        self.nside            = hp.npix2nside(len(self.completeness))
        if diagnostics or self.diagnostics:
            self.check_diagnostics_directory()
            plt.clf()
            hp.mollview(self.completeness, title="Completeness")
            plt.savefig("diagnostics/completeness.pdf")

    """
    Function to  generate r values uniformly distributed with (r_min, r_max)

    Parameters
    ----------
    num_obs : int
        number of r values to generate
    
    diagnostics : boolean

    Returns
    -------
    array(dtype=flat)
        
    """
    def generate_uniform_r(self, num_obs=None, diagnostics=False):
        if num_obs is None:
            num_obs = self.n_center
        flat_r = np.random.uniform(self.template_r_min, self.template_r_max, num_obs)
        # if diagnostics enabled, plot the distribution
        if diagnostics or self.diagnostics:
            self.check_diagnostics_directory()
            plt.hist(flat_r)
            plt.xlabel("r [Mpc]")
            plt.title("Flat r distribution")
            plt.savefig("diagnostics/r_distribution.pdf")        
        return flat_r

    """
    Function to generate r values based on histogrammed data r distribution
    This function is used to generate r distribution depending on the
    template provided
    
    Parameters
    ----------
    num_obs : int
        number of r values to generate
    
    diagnostics : boolean

    Returns
    -------
    array(dtype=flat)
    
    """
    def generate_r(self, nobs, diagnostics=False):
        num_obs = 0
        rlist   = []
        if self.template_w is not None:
            n, r_edges  = np.histogram(self.template_r, bins=self.template_r_len, weights=self.template_w)
        else:
            n, r_edges  = np.histogram(self.template_r, bins=int(np.sqrt(self.template_r_len)))
        n           = n.astype(float)
        weights     = n/np.sum(n)
        r_med       = (r_edges[:-1] + r_edges[1:])/2.
        halfbinsize = (r_edges[1] - r_edges[0])/2.
        while num_obs < nobs:
            curr_r_center = np.random.choice(r_med, p=weights)
            curr_r        = np.random.uniform(curr_r_center-halfbinsize, curr_r_center+halfbinsize)
            rlist.append(curr_r)
            num_obs += 1
        rlist = np.array(rlist)
        # plot the distribution for diagnostics
        if diagnostics or self.diagnostics:
            self.check_diagnostics_directory()
            plt.hist(rlist)
            plt.xlabel("r [Mpc]")
            plt.savefig("diagnostics/r_distribution.pdf")
        return rlist

    def check_diagnostics_directory(self):
        if not os.path.isdir("diagnostics"):
                os.makedirs("diagnostics")

    """
    Function that return the indices of the r values that passes
    the acceptance test
    
    Parameters
    ----------
    r_test : double
        r values to be test for acceptance

    Returns
    -------
    boolean array
        a boolean array which has True values for the r values that passed the acceptance test
        and False values for the r values that did not pass the acceptance test
    
    """
    def check_radial_acceptance(self, r_test):
        template_r   = np.array(self.template_r)
        template_r   = template_r[template_r<=self.r_max]
        num_bins     = int(np.sqrt(len(template_r)))
        n, r_edges   = np.histogram(template_r, bins=num_bins, normed=True)
        n            = n.astype(float)
        r_med        = (r_edges[:-1] + r_edges[1:])/2.
        halfbinsize  = (r_edges[1] - r_edges[0])/2.
        interpolator = interp1d(r_med, n, bounds_error=False, fill_value=0.)

        num_r_test           = len(r_test)
        n_test, r_test_edges = np.histogram(r_test, bins=int(np.sqrt(num_r_test)), normed=True)
        n_test               = n_test.astype(float)
        r_test_med           = (r_test_edges[:-1] + r_test_edges[1:])/2.
        interpolator_test    = interp1d(r_test_med, n_test, bounds_error=False, fill_value=0.)

        p1  = interpolator(r_test)
        p2  = interpolator_test(r_test)
        acc = [np.random.uniform(0, p) for p in p2]

        if self.acceptance:
            return p1>acc
        return np.full(len(p1), True)

    """
    Function that return the indices of the z values that passes
    the acceptance test
    
    Parameters
    ----------
    z_test : double
        z values to be test for acceptance

    Returns
    -------
    boolean array
        a boolean array which has True values for the z values that passed the acceptance test
        and False values for the r values that did not pass the acceptance test
    
    """
    def check_z_acceptance(self, z_test):
        template_z   = np.array(self.template_z)
        if self.template_w is not None:
            template_w   = self.template_w[np.logical_and(template_z<=self.z_max, template_z>=self.z_min)]
            template_z   = template_z[np.logical_and(template_z<=self.z_max, template_z>=self.z_min)]
            num_bins     = len(template_w)
            n, z_edges   = np.histogram(template_z, bins=num_bins, normed=True, weights=template_w,
                                        range=(self.z_min, self.z_max))
        else:
            template_z   = template_z[np.logical_and(template_z<=self.z_max, template_z>=self.z_min)]
            num_bins     = int(np.sqrt(len(template_z)))
            n, z_edges   = np.histogram(template_z, bins=num_bins, normed=True, range=(self.z_min, self.z_max))
        n            = n.astype(float)
        z_med        = (z_edges[:-1] + z_edges[1:])/2.
        halfbinsize  = (z_edges[1] - z_edges[0])/2.
        interpolator = interp1d(z_med, n, bounds_error=False, fill_value=-1.)
        
        num_z_test           = len(z_test)
        n_test, z_test_edges = np.histogram(z_test, bins=int(np.sqrt(num_z_test)),
                                            range=(self.z_min, self.z_max), normed=True)
        n_test               = n_test.astype(float)
        z_test_med           = (z_test_edges[:-1] + z_test_edges[1:])/2.
        interpolator_test    = interp1d(z_test_med, n_test, bounds_error=False, fill_value=0.)
        
        p1  = interpolator(z_test)
        p2  = interpolator_test(z_test)
        acc = [np.random.uniform(0, p) for p in p2]

        passed_acceptance = np.logical_and(np.logical_and(p1>acc, z_test<=self.z_max), z_test>=self.z_min)
        print(np.min(z_test[passed_acceptance]), np.max(z_test[passed_acceptance]))
        
        if self.acceptance:
            print("applying the z acceptance...")
            return passed_acceptance
        print("not applying z acceptance...")
        return np.full(len(p1), True)

    """
    Function that returns randomly generated (theta, phi)'s using the 
    completeness map defined in the configuration

    Parameters
    ----------
    nobs: int
        number of (theta, phi) pair to generate
    diagnostics: bool
        
    Returns
    -------
    array of arrays
        The first column of the returned array is for theta and the second column is for  the phi
        coordinates. These coordinates are later transformed into RA and Dec while writing to the
        fits file.
 
    """
    def generate_uniform_angular_position(self, nobs, diagnostics=False):
        num_obs  = 0
        thetas   = []
        phis     = []
        pixels   = []
        weights  = self.completeness/np.sum(self.completeness)
        pix_idx  = np.arange(self.completeness_len)
        while num_obs < nobs:
            curr_phi   = np.random.uniform(0., 360., 1)
            curr_theta = np.arccos(np.random.uniform(-1., 1., 1))*RAD2DEG-90.
            curr_pix   = hp.ang2pix(self.nside, curr_phi, curr_theta, lonlat=True)
            if self.completeness[curr_pix] == 1:
                thetas.append(curr_theta)
                phis.append(curr_phi)
                pixels.append(curr_pix)
                num_obs += 1
        # plot the distribution for diagnostics
        if diagnostics or self.diagnostics:
            self.check_diagnostics_directory()
            # completeness plot
            plt.clf()
            pixels = np.zeros([12*self.nside**2])
            pixels[pixlist] = 1
            hp.mollview(pixels.astype(int), title="Completeness")
            self.check_diagnostics_directory()
            plt.savefig("diagnostics/generated_completeness.pdf")
            # theta distribution
            plt.clf()
            plt.hist(thetas, bins=int(np.sqrt(self.completeness_len)))
            plt.title(r"$\theta$ distribution")
            plt.xlabel(r'$\theta$ [deg]')
            plt.savefig("diagnostics/generated_theta_dist.pdf")
            # theta distribution
            plt.clf()
            plt.hist(phis, bins=int(np.sqrt(self.completeness_len)))
            plt.title(r"$\phi$ distribution")
            plt.xlabel(r'$\phi$ [deg]')
            plt.savefig("diagnostics/generated_phi_dist.pdf")
        return [np.array(thetas), np.array(phis)]
            
    # generate gaussian distributed r values
    def generate_gaussian(self, mu, sig, nobs):
        rlist = []
        for nob in range(nobs):
            val = random.gauss(mu,sig)
            rlist.append(val)
        rlist = np.array(rlist)
        return rlist

    # generate random positions in the sky
    def generate_galaxies(self, num_obs, uniform=False):
        if uniform:
            r_list = self.generate_uniform_r(num_obs)
        else:
            r_list = self.generate_r(num_obs)
        angles_list = self.generate_uniform_angular_position(num_obs)
        return r_list, angles_list[0], angles_list[1]

    def generate_center_galaxies(self):
        r_center, theta_center, phi_center = self.generate_galaxies(self.n_center)
        center_galaxies = {}
        for i in range(self.n_center):
            center_galaxies["cen_{}".format(i)] = galaxy(theta=theta_center[i], phi=phi_center[i], r=r_center[i], TYPE=0)
        self.catalog.centers = center_galaxies
        return
        
    # function to convert spherical coordinates to cartesian coordinates
    def toCartesianVector(self, r, theta, phi):
        #unit_vector = hp.pix2vec(self.nside, hp.ang2pix(self.nside, theta, phi, lonlat=True))
        unit_vector = hp.ang2vec(phi, theta, lonlat=True)
        curr_vector = r * np.array(unit_vector)
        return curr_vector

    # function to convert spherical coordinates to cartesian coordinates with eccentricity (to mimick anisotropy)
    def toCartesianVector2(self, r, theta, phi, center_theta=None, center_phi=None):
        #unit_vector = hp.pix2vec(self.nside, hp.ang2pix(self.nside, theta, phi, lonlat=True))
        unit_vector = hp.ang2vec(phi, theta, lonlat=True)
        curr_vector = np.asarray(unit_vector) * self.anisotropy_scale
        if (center_theta is not None) and (center_phi is not None):
            angles = hp.vec2ang(curr_vector, lonlat=True)
            # rotate the positions with respect to the corresponding centers
            rotated_angles = angles + np.asarray([center_phi, center_theta])
            # in some cases, angles need to be corrected for the boundaries
            if rotated_angles[0] > 360.:
                rotated_angles[0] -= 360.
            if rotated_angles[0] < 0.:
                rotated_angles[0] += 360.
            if rotated_angles[1] < 0.:
                rotated_angles[1] += 180.
            if rotated_angles[1] > 180.:
                rotated_angles[1] -= 180.
            curr_vector = np.asarray(hp.ang2vec(rotated_angles[1]*DEG2RAD, rotated_angles[0]*DEG2RAD))
        # multiply by the distance to get the proper coordinates
        return r*curr_vector

    # function to calculate the rotation matrix for a given galaxy
    def calculateRotationMatrix(self, theta, phi, lonlat=True):
        curr_position = hp.ang2vec(phi, theta, lonlat=lonlat)[0]
        # for convenience, we will keep the units in radians for now
        rot_x         = np.arccos(np.dot([1., 0., 0.], curr_position))
        rot_y         = np.arccos(np.dot([0., 1., 0.], curr_position))
        rot_z         = np.arccos(np.dot([0., 0., 1.], curr_position))
        R_x           = np.asarray([[1., 0., 0.],
                                    [0., np.cos(rot_x), -np.sin(rot_x)],
                                    [0., np.sin(rot_x), np.cos(rot_x)]])
        R_y           = np.asarray([[np.cos(rot_y),  0., np.sin(rot_y)],
                                    [0.,             1., 0.],
                                    [-np.sin(rot_y), 0., np.cos(rot_y)]])
        R_z           = np.asarray([[np.cos(rot_z), -np.sin(rot_z), 0.],
                                    [np.sin(rot_z), np.cos(rot_z),  0.],
                                    [0.,            0.,             1.]])
        rotation_matrix = np.matmul(R_z, np.matmul(R_y, R_x))
        return rotation_matrix
        
    # function to convert spherical coordinates to cartesian coordinates with eccentricity (to mimick anisotropy)
    def toCartesianVector3(self, r, theta, phi, center_theta=None, center_phi=None):
        unit_vector = hp.ang2vec(phi, theta, lonlat=True)
        curr_vector = np.asarray(r * self.anisotropy_scale * np.asarray(unit_vector))
        if center_theta is not None and center_phi is not None:
            rotationMatrix = self.calculateRotationMatrix(center_theta, center_phi)
            return np.matmul(rotationMatrix, curr_vector.transpose())
        return curr_vector
        
    # function to convert cartesian coordinates to spherical coordinates
    def fromCartesianVector(self, vec):
        r           = np.sqrt( np.sum(np.array(vec)**2) )
        unit_vector = hp.vec2ang(vec, lonlat=True)
        return r, unit_vector[1], unit_vector[0]
        
    # function to add two spherical vectors
    def addVectors(self, vec1, vec2):
        vec1_cartesian  = self.toCartesianVector(*vec1)
        vec2_cartesian  = self.toCartesianVector(*vec2)
        sumvec          = self.fromCartesianVector(vec1_cartesian+vec2_cartesian)
        # this part is necessary to get the added vector to the right quadrant
        # it is not intended for general use but specific to this project
        return sumvec

    def read_generated_file(self, filename):
        hdus = fits.open(filename)
        data = hdus[1].data
        z     = data['z']
        r     = [(self.cosmo.comoving_distance(curr_z)).value for curr_z in z]
        theta = data['theta']
        phi   = data['phi']
        return np.array(r), np.array(theta), np.array(phi)
    
    #
    def generate_rim_from_file(self, filename, diagnostics=False):
        r, theta, phi = self.read_generated_file(filename)
        return self.generate_rim(r, theta, phi, diagnostics)

    #
    def generate_rim(self):
        if self.catalog.rims is None:
            rim_galaxies = {}
        else:
            rim_galaxies = self.catalog.rims
            
        for i in range(self.n_center):
            curr_center_galaxy = self.catalog.centers["cen_{}".format(i)]
            if curr_center_galaxy.childs:
                curr_center_galaxy_childs = curr_center_galaxy.childs
            else:
                curr_center_galaxy_childs = []
            curr_rim_cnt = 0

            while curr_rim_cnt < self.n_rim:
                curr_phi   = np.random.uniform(0., 360., 1)
                curr_theta = np.arccos(np.random.uniform(-1., 1., 1))*RAD2DEG-90.
                curr_r     = self.generate_gaussian(self.r_BAO, self.sigma_r_BAO, 1)
                curr_rim_wrt_center = self.toCartesianVector3(curr_r[0], curr_theta[0], curr_phi[0],
                                                              center_theta=curr_center_galaxy.theta, center_phi=curr_center_galaxy.phi)
                curr_rim   = (self.toCartesianVector(curr_center_galaxy.r, curr_center_galaxy.theta, curr_center_galaxy.phi) + \
                              curr_rim_wrt_center )[0]
                pixel      = hp.vec2pix(self.nside, x=curr_rim[0], y=curr_rim[1], z=curr_rim[2])
                # apply angular acceptance.
                if self.completeness[pixel] == 1:
                    curr_rim = self.fromCartesianVector(curr_rim)
                    rim_galaxies["rim_{}_{}".format(i, curr_rim_cnt)] = galaxy(theta=curr_rim[1], phi=curr_rim[2], r=curr_rim[0],
                                                                               parent="cen_{}".format(i), TYPE=1)
                    curr_center_galaxy_childs.append("rim_{}_{}".format(i, curr_rim_cnt))
                    curr_rim_cnt += 1                
            curr_center_galaxy.childs = curr_center_galaxy_childs
        if self.catalog.rims is None:
            self.catalog.rims = rim_galaxies
        else:
            self.catalog.rims.update(rim_galaxies)
        return

    def generate_flat_galaxies(self, is_random=False):
        if is_random is False:
            r_flat, theta_flat, phi_flat = self.generate_galaxies(self.n_flat)
        else:
            print("generating flat galaxies for random catalog")
            r_flat, theta_flat, phi_flat = self.generate_galaxies(self.n_rand)
        flat_galaxies = {}
        for i, _ in enumerate(r_flat):
            flat_galaxies["flat_{}".format(i)] = galaxy(theta=theta_flat[i], phi=phi_flat[i], r=r_flat[i], TYPE=2)
        self.catalog.flats = flat_galaxies
        return r_flat, theta_flat, phi_flat
    
    def generate_clumps_from_file(self, filename):
        r, theta, phi = self.read_generated_file(filename)
        return self.generate_clumps(r, theta, phi)
    
    def generate_clumps(self, add_clumps_to_rims = False):
        # generate flat galaxies (will be returned and be added to the mocks later)
        if self.catalog.flats is None:
            self.generate_flat_galaxies()
        total_num_galaxies      = self.n_center + self.n_flat
        total_num_centers       = self.n_center
        if add_clumps_to_rims:
            total_num_galaxies += self.n_center * self.n_rim
            total_num_centers  += self.n_center * self.n_rim
        # generate a list to choose among the centers and flats
        if self.frac_f2c is None:
            galaxy_selection  = np.random.choice(np.arange(0, total_num_galaxies, 1), size=self.nr_clump)
            num_center_clumps = len(galaxy_selection[galaxy_selection<total_num_centers])
            num_flat_clumps   = len(galaxy_selection[galaxy_selection>=total_num_centers])
        else:
            num_center_clumps = int(self.nr_clump / ( 1 + self.frac_f2c ))
            num_flat_clumps   = self.nr_clump - self.center_clumps
        num_center_clumps = num_center_clumps if num_center_clumps > 1 else 2
        num_flat_clumps   = num_flat_clumps if num_flat_clumps > 1 else 2        
        # randomly choose indices from centers for the clumps
        # generate the clump positions with respect to their origin (a center galaxu)
        if self.frac_c2r is None:
            galaxy_selection  = np.random.choice(np.arange(0, total_num_centers, 1), size=num_center_clumps)
            center_clumps_idx = galaxy_selection[galaxy_selection<self.n_center]
            rim_clumps_idx    = galaxy_selection[galaxy_selection>=self.n_center]-self.n_center
        else:
            center_clumps_cnt = int(num_center_clumps / ( 1 + self.frac_c2r ))
            rim_clumps_cnt    = num_center_clumps - center_clumps_cnt
            center_clumps_idx = np.random.choice(np.arange(0, self.n_center, 1), size=center_clumps_cnt)
            rim_clumps_idx    = np.random.choice(np.arange(0, self.n_center*self.n_rim, 1), size=rim_clumps_cnt)

        # get the center clump object
        if self.catalog.clumps_center is None:
            clumps_center = {}
        else:
            clumps_center = self.catalog.clumps_center
        
        # generate the clumps around the center galaxies
        for idx in center_clumps_idx:
            # get the seed galaxy for placing clumps
            curr_seed_galaxy = self.catalog.centers["cen_{}".format(idx)]
            # get the list of already existing childs of the seed galaxy
            # if there is none, then use an empty list
            if curr_seed_galaxy.childs:
                curr_seed_galaxy_childs = curr_seed_galaxy.childs
            else:
                curr_seed_galaxy_childs = []
            # generate clump positions with respect to the seed galaxy
            if self.clump_dist == 0:
                n_clump_to_inject = self.n_clump_center
            else:
                n_clump_to_inject = np.random.poisson(self.n_clump_center)
            clump_r          = (self.r_0**self.gamma * (np.random.pareto(self.gamma-1, n_clump_to_inject)))#self.n_clump_center)))
            clump_phi        = np.random.uniform(0., 360., n_clump_to_inject)#self.n_clump_center)
            clump_theta      = np.arccos(np.random.uniform(-1., 1., n_clump_to_inject))*RAD2DEG-90.#self.n_clump_center))*RAD2DEG-90.
            for j in range(n_clump_to_inject):
                # calculate the absolute position of the clump galaxy
                curr_clump   = (self.toCartesianVector(curr_seed_galaxy.r, curr_seed_galaxy.theta, curr_seed_galaxy.phi)[0] + \
                                self.toCartesianVector(clump_r[j], clump_theta[j], clump_phi[j]))
                # get the healpix pixel position to check with completeness map
                pixel      = hp.vec2pix(self.nside, x=curr_clump[0], y=curr_clump[1], z=curr_clump[2])
                if self.completeness[pixel] == 1:
                    # convert the cartesian coordinates to spherical coordinates
                    curr_clump = self.fromCartesianVector(curr_clump)
                    curr_seed_galaxy_childs.append("cenClump_{}_-1_{}".format(idx, j))
                    # append the new generated clump the object
                    clumps_center["cenClump_{}_-1_{}".format(idx, j)] = galaxy(theta=curr_clump[1], phi=curr_clump[2], r=curr_clump[0],
                                                                               parent="cen_{}".format(idx), TYPE=3)
            # replace the child list of the seed galaxy.
            # it is safe, the new childs are appended to the existing list of childs
            curr_seed_galaxy.childs = curr_seed_galaxy_childs

        # generate the clumps around the rim galaxies (if they will be added)
        for idx in rim_clumps_idx:
            # calculate the index of the center galaxy of the corresponding rim galaxy
            seed_idx1 = int((idx-self.n_center) % self.n_center)
            # calculate the index of the rim galaxy with respect to the center galaxy
            seed_idx2 = int((idx-self.n_center) / self.n_center)
            # get the seed galaxy
            curr_seed_galaxy = self.catalog.rims["rim_{}_{}".format(seed_idx1, seed_idx2)]
            # get the list of childs of the seed galaxy
            # if there is none, then use an empty list
            if curr_seed_galaxy.childs:
                curr_seed_galaxy_childs = curr_seed_galaxy.childs
            else:
                curr_seed_galaxy_childs = []
            # generate clump positions with respect to the seed galaxy
            if self.clump_dist == 0:
                n_clump_to_inject = self.n_clump_center
            else:
                n_clump_to_inject = np.random.poisson(self.n_clump_center)
            clump_r          = (self.r_0**self.gamma * (np.random.pareto(self.gamma-1, n_clump_to_inject)))#self.n_clump_center)))
            clump_phi        = np.random.uniform(0., 360., n_clump_to_inject)#self.n_clump_center)
            clump_theta      = np.arccos(np.random.uniform(-1., 1., n_clump_to_inject))*RAD2DEG-90.#self.n_clump_center))*RAD2DEG-90.
            for j in range(n_clump_to_inject):
                # calculate the absolute position of the clump galaxy
                curr_clump   = (self.toCartesianVector(curr_seed_galaxy.r, curr_seed_galaxy.theta, curr_seed_galaxy.phi)[0] + \
                                self.toCartesianVector(clump_r[j], clump_theta[j], clump_phi[j]))
                # get the healpix pixel position to check with completeness map
                pixel      = hp.vec2pix(self.nside, x=curr_clump[0], y=curr_clump[1], z=curr_clump[2])
                if self.completeness[pixel] == 1:
                    # convert the cartesian coordinates to spherical coordinates
                    curr_clump = self.fromCartesianVector(curr_clump)
                    curr_seed_galaxy_childs.append("cenClump_{}_{}_{}".format(seed_idx1, seed_idx2, j))
                    # append the new generated clump the object
                    clumps_center["cenClump_{}_{}_{}".format(seed_idx1, seed_idx2, j)] = galaxy(theta=curr_clump[1], phi=curr_clump[2], r=curr_clump[0],
                                                                                                parent="rim_{}_{}".format(seed_idx1, seed_idx2), TYPE=3)
            # replace the child list of the seed galaxy.
            # it is safe, the new childs are appended to the existing list of childs
            curr_seed_galaxy.childs = curr_seed_galaxy_childs
        # update the catalog object as all the center clumps are generated
        self.catalog.clumps_center = clumps_center
        
        # get the center clump object
        if self.catalog.clumps_flat is None:
            clumps_flat = {}
        else:
            clumps_flat = self.catalog.clumps_flat

        # randomly choose indices from flats for the clumps
        # generate the clump positions with respect to their origin (a flat galaxy)
        rand_flat_idx          = np.random.randint(0, self.n_flat, num_flat_clumps)
        for idx in rand_flat_idx:
            # get the seed galaxy for placing clumps
            curr_seed_galaxy = self.catalog.flats["flat_{}".format(idx)]
            # get the list of already existing childs of the seed galaxy
            # if there is none, then use an empty list
            if curr_seed_galaxy.childs:
                curr_seed_galaxy_childs = curr_seed_galaxy.childs
            else:
                curr_seed_galaxy_childs = []
            # generate clump positions with respect to the seed galaxy
            if self.clump_dist == 0:
                n_clump_to_inject = self.n_clump
            else:
                n_clump_to_inject = np.random.poisson(self.n_clump)            
            clump_r          = (self.r_0**self.gamma * (np.random.pareto(self.gamma-1, n_clump_to_inject)))#self.n_clump)))
            clump_phi        = np.random.uniform(0., 360., n_clump_to_inject)#self.n_clump)
            clump_theta      = np.arccos(np.random.uniform(-1., 1., n_clump_to_inject))*RAD2DEG-90.#self.n_clump))*RAD2DEG-90.
            for j in range(n_clump_to_inject):
                # calculate the absolute position of the clump galaxy
                curr_clump   = (self.toCartesianVector(clump_r[j], clump_theta[j], clump_phi[j]) + \
                                self.toCartesianVector(curr_seed_galaxy.r, curr_seed_galaxy.theta, curr_seed_galaxy.phi))[0]
                # get the healpix pixel position to check with completeness map
                pixel        = hp.vec2pix(self.nside, x=curr_clump[0], y=curr_clump[1], z=curr_clump[2])
                # apply angular acceptance
                if self.completeness[pixel] == 1:
                    # convert the cartesian coordinates to spherical coordinates
                    curr_clump = self.fromCartesianVector(curr_clump)
                    curr_seed_galaxy_childs.append("flatClump_{}_-1_{}".format(idx, j))
                    # append the new generated clump the object
                    clumps_flat["flatClump_{}_-1_{}".format(idx, j)] = galaxy(theta=curr_clump[1], phi=curr_clump[2], r=curr_clump[0],
                                                                              parent="flat_{}".format(idx), TYPE=4)
            # replace the child list of the seed galaxy.
            # it is safe, the new childs are appended to the existing list of childs
            curr_seed_galaxy.childs = curr_seed_galaxy_childs
            # update the catalog object as all the center clumps are generated
        self.catalog.clumps_flat = clumps_flat

    def r2z(self, r):
        return [z_at_value(self.cosmo.comoving_distance, curr_r*u.Mpc) for curr_r in r]

    def write_to_pickle(self, filename=None):
        if filename is None:
            filename = self.fname_mock.replace(".fits", ".pkl")
        config = {}
        config['H_0']         = self.H0
        config['omega_m0']    = self.omega_m0
        config['r_BAO']       = self.r_BAO
        config['sig_r_BAO']   = self.sigma_r_BAO
        config['gamma']       = self.gamma
        config['r_0']         = self.r_0
        config['n_rand']      = self.n_rand
        config['n_center']    = self.n_center
        config['n_rim']       = self.n_rim
        config['n_flat']      = self.n_flat
        config['nr_cl']       = self.nr_clump
        config['n_cl']        = self.n_clump
        config['n_cl_center'] = self.n_clump_center
        config['frac_f2c']    = self.frac_f2c
        config['frac_c2r']    = self.frac_c2r
        save_items = {'catalog': self.catalog, 'config': config}
        pickle.dump(save_items, open(filename, "wb"), protocol=-1)

    def write_to_fits(self, filename=None, is_random=False):
        if filename is None:
            if not is_random:
                filename = self.fname_mock
            else:
                filename = self.fname_random
        rs, ras, decs, types = self.catalog.flatten()
        LUT = self.generate_LUT_r2z()
        zs  = LUT(rs)
        ws  = np.ones(len(zs))
        # acceptance test here. should it be done for individual types?
        if self.acceptance:
            accepted_indices = self.check_z_acceptance(zs)
            zs    = np.asarray(zs)[accepted_indices]
            ras   = np.asarray(ras)[accepted_indices]
            decs  = np.asarray(decs)[accepted_indices]
            types = np.asarray(types)[accepted_indices]
            ws    = np.asarray(ws)[accepted_indices]
        # We also write the output in fits format
        if os.path.isfile(filename):
            print("a file with the designated name already exists... please remove the file first")
            return
        header = fits.Header()
        header['H_0']         = self.H0
        header['omega_m0']    = self.omega_m0
        header['r_BAO']       = self.r_BAO
        header['sig_r_BAO']   = self.sigma_r_BAO
        header['gamma']       = self.gamma
        header['r_0']         = self.r_0
        header['n_rand']      = self.n_rand
        header['n_center']    = self.n_center
        header['n_rim']       = self.n_rim
        header['n_flat']      = self.n_flat
        header['nr_cl']       = self.nr_clump
        header['n_cl']        = self.n_clump
        header['n_cl_center'] = self.n_clump_center
        header['a']           = self.a
        header['b']           = self.b
        col1 = fits.Column(name="z", array=zs, format='E')
        col2 = fits.Column(name="ra", array=ras, format='E')
        col3 = fits.Column(name="dec", array=decs, format='E')
        col4 = fits.Column(name="weight", array=ws, format='E')
        col5 = fits.Column(name="TYPE", array=types, format='J')
        cols = fits.ColDefs([col1, col2, col3, col4, col5])
        hdu  = fits.BinTableHDU.from_columns(cols, header=header)
        hdu.writeto(filename)

def generate_histo_from_cat(inFilename, outFilename):
    # reads the datafile to generate a histogram
    # this works on any DESI datachallange
    # it may need to be modified later for other experiments (maybe)
    hdus            = fits.open(inFilename)
    zs              = np.array(hdus[1].data['Z'])
    entries, edges  = np.histogram(zs, bins=int(np.sqrt(len(zs))))
    # after getting the histogram, generate a new fits file
    # with the bin centers and entries (as weights)
    col1 = fits.Column(name="Z",      array=(edges[:-1]+edges[1:])/2., format="E")
    col2 = fits.Column(name="HistoW", array=entries, format="E")
    cols = fits.ColDefs([col1, col2])
    hdu  = fits.BinTableHDU.from_columns(cols)
    hdu.writeto(outFilename)
