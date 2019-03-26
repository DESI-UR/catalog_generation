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

#from paramock import galaxy
#from paramock import catalog
from galaxy import galaxy
from catalog import catalog

import time

NSIDE = -1
DEG2RAD = np.pi/180.0
RAD2DEG = 1./DEG2RAD

class GeneralTools():

    def __init__(self, configFile, diagnostics=False, acceptance=True, use_style=True):
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
        self.diagnostics = diagnostics
        self.acceptance  = acceptance
        self.config_file = configFile
        self.catalog     = catalog()
        self.getConfig()
        self.get_template()
        self.get_mask()
        self.flat_clump_r_dist, _   = np.histogram([], bins=100, range=(0, 100))
        self.center_clump_r_dist, _ = np.histogram([], bins=100, range=(0, 100))
        self.rim_clump_r_dist, _    = np.histogram([], bins=100, range=(0, 100))
        self.rim_BAO_r_dist, _      = np.histogram([], bins=100, range=(100, 200))
        if use_style:
            import style
        
    def setDiagnostics(self, diagnostics):
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
        self.diagnostics = diagnostics

    def read_config_line(self, config_line):
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
        try:
            quantity = u.Quantity(config_line)
        except Exception as e:
            raise Exception('the value provided for the following line is invalid'
                            '{}'
                            'Please check for correct syntax, either a real number or a quantity with unit'.format(config_line))
        return quantity
    
    def getConfig(self, diagnostics=False):
        """
        Function to set the parameters of the object using the configuration file.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
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
            self.h0    = H0.value
            self.H0    = H0.value * 100.
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
            self.r_BAO      = r_BAO.value/self.h0
        else:
            self.r_BAO      = r_BAO.to("Mpc").value
        sigma_r_BAO         = self.read_config_line(self.config.get('Gen Params','sigr_bao'))
        if sigma_r_BAO.unit == u.dimensionless_unscaled:
            self.sigma_r_BAO= sigma_r_BAO.value/self.h0
        else:
            self.sigma_r_BAO= sigma_r_BAO.to("Mpc").value

        # Clumping parameters
        self.gamma          = float(self.config.get('Gen Params','gamma'))
        r_scale             = self.read_config_line(self.config.get('Gen Params','r_scale'))
        if r_0.unit == u.dimensionless_unscaled:
            self.r_scale    = r_scale.value/self.h0
        else:
            self.r_scale    = r_scale.to("Mpc").value
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
            self.frac_f2c   = float(self.config.get('Gen Params','frac_f2c'))
        except Exception as e:
            self.frac_f2c   = None
        try:
            self.frac_c2r   = float(self.config.get('Gen Params','frac_c2r'))
        except Exception as e:
            self.frac_c2r   = None
        try:
            self.a          = float(self.config.get('Gen Params','a'))
            self.b          = float(self.config.get('Gen Params','b'))
            self.c          = float(self.config.get('Gen Params','c'))
        except Exception as e:
            self.a          = 1.
            self.b          = 1.
            self.c          = 1.
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
        
        self.anisotropy_scale = np.array([self.a, self.b, self.c])

        if self.diagnostics or diagnostics:
            print("Configuration is a follows:\n"\
                  " - Cosmology: H0={}, omegaM0={}\n"\
                  " - Limits: z=({}, {})\n"\
                  " - Anisotropy scale: {}\n"\
                  " - Number of galaxies to be introduced:\n"\
                  "   - N(flat): {}\n"\
                  "   - N(center): {}\n"\
                  "   - N(rims per center): {}\n"\
                  "   - N(clumping): {}\n"\
                  "   - N(clumps per center): {}\n"\
                  "   - N(clumps per flat): {}\n"\
                  "".format(self.H0, self.omega_m0, self.z_min, self.z_max, self.anisotropy_scale,
                            self.n_flat, self.n_center, self.n_rim, self.nr_clump, self.n_clump_center, self.n_clump))

    def generate_LUT_z2r(self):
        """
        Function to generate a lookup table for z to r conversion
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        z_min = 0.1
        z_max = 3.0
        zs = np.linspace(z_min, z_max, 1000)
        rs = [r.value for r in self.cosmo.comoving_distance(zs)]
        interpolator = interp1d(zs, rs, bounds_error=False, fill_value=-1.)
        return interpolator

    def generate_LUT_r2z(self):
        """
        Function to generate a lookup table for r to z conversion
    
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
        z_min = 0.1
        z_max = 3.0
        zs = np.linspace(z_min, z_max, 1000)
        rs = [r.value for r in self.cosmo.comoving_distance(zs)]
        interpolator = interp1d(rs, zs, bounds_error=False, fill_value=-9999.)
        return interpolator            

    def get_template(self):
        """
        Function to read the template r distribution for mock and random generation
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """
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

    def get_mask(self, diagnostics=False):
        """
        Function to read the completeness map for mock and random generation

        Parameters
        ----------
        diagnostics : boolean
        
        Returns
        -------
        None
        """
        self.mask             = np.load(self.ang_mask)
        self.completeness     = self.mask.f.CMPLTNSS
        self.completeness_len = len(self.completeness)
        self.nside            = hp.npix2nside(len(self.completeness))
        if diagnostics or self.diagnostics:
            self.check_diagnostics_directory()
            fig = plt.figure()
            cvals = [0, 1]
            colors= ['white', 'black']
            norm  = plt.Normalize(min(cvals), max(cvals))
            tuples= list(zip(map(norm,cvals), colors))
            cmap  = mpl.colors.LinearSegmentedColormap.from_list("", tuples)
            hp.mollview(self.completeness, title="Completeness", cbar=False, flip='geo', margins=(0,0,0,0), cmap=cmap)
            hp.graticule()
            plt.savefig("diagnostics/completeness.pdf")
            plt.close()
            
    def generate_uniform_r(self, num_obs=None, diagnostics=False):
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
        if num_obs is None:
            num_obs = self.n_center
        flat_r = np.random.uniform(self.template_r_min, self.template_r_max, num_obs)
        # if diagnostics enabled, plot the distribution
        if diagnostics or self.diagnostics:
            self.check_diagnostics_directory()
            fit = plt.figure()
            plt.hist(flat_r, histtype='step', bins=int(np.sqrt(len(flat_r))))
            plt.xlabel("r [Mpc]")
            plt.title("Flat r distribution")
            plt.savefig("diagnostics/flat_r_distribution.pdf")
            plt.close()
        return flat_r

    def generate_r(self, nobs, diagnostics=False):
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
            fig = plt.figure()
            plt.hist(rlist*self.h0, histtype='step', bins=int(np.sqrt(len(rlist))))
            plt.xlabel("r [$h^{-1}$Mpc]")
            plt.savefig("diagnostics/r_distribution.pdf")
            plt.close()
        return rlist

    def check_diagnostics_directory(self):
        if not os.path.isdir("diagnostics"):
                os.makedirs("diagnostics")

    def check_radial_acceptance(self, r_test):
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

    def check_z_acceptance(self, z_test, diagnostics=False):
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
        template_z   = np.array(self.template_z)
        if self.template_w is not None:
            template_w   = self.template_w[np.logical_and(template_z<=self.z_max, template_z>=self.z_min)]
            template_z   = template_z[np.logical_and(template_z<=self.z_max, template_z>=self.z_min)]
            num_bins     = len(template_w)
            n, z_edges   = np.histogram(template_z, bins=num_bins, normed=True, weights=template_w,
                                        range=(self.z_min, self.z_max))
        else:
            template_z   = template_z[np.logical_and(template_z<=self.z_max, template_z>=self.z_min)]
            num_bins     = self.calculate_n_bins(len(template_z)) 
            n, z_edges   = np.histogram(template_z, bins=num_bins, normed=True, range=(self.z_min, self.z_max))
        n            = n.astype(float)
        z_med        = (z_edges[:-1] + z_edges[1:])/2.
        halfbinsize  = (z_edges[1] - z_edges[0])/2.
        interpolator = interp1d(z_med, n, bounds_error=False, fill_value=-1.)
        
        num_z_test           = len(z_test)
        num_z_test_bins      = self.calculate_n_bins(num_z_test) 
        n_test, z_test_edges = np.histogram(z_test, bins=num_z_test_bins,
                                            range=(self.z_min, self.z_max), normed=True)
        n_test               = n_test.astype(float)
        z_test_med           = (z_test_edges[:-1] + z_test_edges[1:])/2.
        interpolator_test    = interp1d(z_test_med, n_test, bounds_error=False, fill_value=0.)
        
        p1  = interpolator(z_test)
        p2  = interpolator_test(z_test)
        acc = [np.random.uniform(0, p) for p in p2]

        passed_acceptance   = np.logical_and(np.logical_and(p1>acc, z_test<=self.z_max), z_test>=self.z_min)
        num_z_accepted_bins = self.calculate_n_bins(len(z_test[passed_acceptance])/2.)
        
        if self.diagnostics or diagnostics:
            # plot the distributions before and after the acceptance test
            # along with the template distribution
            self.check_diagnostics_directory()
            fig = plt.figure()
            n1, b1, _ = plt.hist(template_z, bins=len(template_w), density=True, weights=template_w, range=(self.z_min, self.z_max),
                                 label="Template distribution", histtype='step', linewidth=2)
            b1_mid    = (b1[1:] + b1[:-1])*.5
            n2, b2, _ = plt.hist(z_test, bins=num_z_test_bins, density=True, range=(self.z_min, self.z_max),
                                 label="Distribution before acceptance", histtype='step', linewidth=2)
            b2_mid    = (b2[1:] + b2[:-1])*.5
            n3, b3, _ = plt.hist(z_test[passed_acceptance], bins=num_z_accepted_bins, density=True, range=(self.z_min, self.z_max),
                                 label="Distribution after acceptance", histtype='step', linewidth=2)
            b3_mid    = (b3[1:] + b3[:-1])*.5
            plt.close()
            fig = plt.figure(figsize=(6,6))
            plt.plot(b1_mid, n1, label="Template distribution")
            plt.plot(b2_mid, n2, label="Distribution before acceptance")
            plt.plot(b3_mid, n3, label="Distribution after acceptance")
            plt.xlabel("Redshift [z]")
            plt.ylim(0, 1.25*max(max(n1), max(n2), max(n3)))
            plt.legend()
            plt.savefig("diagnostics/acceptance_stats.pdf")
            plt.close()
            
        if self.acceptance:
            print("applying the z acceptance...")
            if self.diagnostics or diagnostics:
                print("Fraction of galaxies passing the acceptance test is: {:.2f}".format(float(np.sum(passed_acceptance))/len(passed_acceptance)))
            return passed_acceptance
        print("not applying z acceptance...")
        return np.full(len(p1), True)

    def generate_uniform_angular_position(self, nobs, diagnostics=False):
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
        if diagnostics:
            self.check_diagnostics_directory()
            # completeness plot
            fig = plt.figure()
            cvals  = [0, 1]
            colors = ["white", "blue"]
            norm   = plt.Normalize(min(cvals),max(cvals))
            tuples = list(zip(map(norm,cvals), colors))
            cmap   = mpl.colors.LinearSegmentedColormap.from_list("", tuples)
            allpixels = np.zeros([12*self.nside**2])
            print(np.asarray(pixels).flatten())
            allpixels[np.asarray(pixels).flatten()] = 1
            hp.mollview(allpixels.astype(int), cbar=False, title='Completeness', flip='geo', margins=(0,0,0,0), cmap=cmap)
            hp.graticule()
            self.check_diagnostics_directory()
            plt.savefig("diagnostics/generated_completeness.pdf")
            plt.close()
            # theta distribution
            fig = plt.figure()
            plt.hist(thetas, bins=int(np.sqrt(self.completeness_len)), histtype='step')
            plt.title(r"$\theta$ distribution")
            plt.xlabel(r'$\theta$ [deg]')
            plt.savefig("diagnostics/generated_theta_dist.pdf")
            plt.close()
            # theta distribution
            fig = plt.figure()
            plt.hist(phis, bins=int(np.sqrt(self.completeness_len)), histtype='step')
            plt.title(r"$\phi$ distribution")
            plt.xlabel(r'$\phi$ [deg]')
            plt.savefig("diagnostics/generated_phi_dist.pdf")
            plt.close()
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
            center_galaxies["cen_{}".format(i)] = galaxy(theta=theta_center[i], phi=phi_center[i], r=r_center[i], TYPE=0, name="cen_{}".format(i))
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
    def generate_rim(self, diagnostics=False):
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

            curr_r_list = []
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
                                                                               parent="cen_{}".format(i), TYPE=1, name="rim_{}_{}".format(i, curr_rim_cnt))
                    curr_center_galaxy_childs.append("rim_{}_{}".format(i, curr_rim_cnt))
                    curr_rim_cnt += 1
                    curr_r_list.append(curr_r)
            # if diagnostics is enabled, add the current distribution of distances with respect to the respective center
            # into the main histogram
            if self.diagnostics or diagnostics:
                rim_BAO_r_dist, _    = np.histogram(curr_r_list, bins=100, range=(100, 200))
                self.rim_BAO_r_dist += rim_BAO_r_dist
            curr_center_galaxy.childs = curr_center_galaxy_childs
        if self.catalog.rims is None:
            self.catalog.rims = rim_galaxies
        else:
            self.catalog.rims.update(rim_galaxies)
        # if diagnostics is enabled, plot the distribution of the rim galaxies with respect to their centers
        if self.diagnostics or diagnostics:
            bins = np.linspace(100, 200, 101)
            bin_centers = (bins[1:]+bins[:-1])*.5*self.h0
            normed_rim_BAO_r_dist = self.rim_BAO_r_dist/np.sum(self.rim_BAO_r_dist)
            fig = plt.figure(figsize=(6,6))
            plt.plot(bin_centers, normed_rim_BAO_r_dist, label='Distribution in the mock', color='black')
            plt.plot([self.r_BAO*self.h0, self.r_BAO*self.h0], [0, normed_rim_BAO_r_dist.max()], '--', color='black',
                     label=r'r$_{{BAO}}$ ({} $h^{{-1}}Mpc$)'.format(self.r_BAO*self.h0))
            plt.plot([(self.r_BAO-1.2*self.sigma_r_BAO)*self.h0, (self.r_BAO+1.2*self.sigma_r_BAO)*self.h0],
                     [normed_rim_BAO_r_dist.max()/np.e, normed_rim_BAO_r_dist.max()/np.e], ':',
                      color='black', label=r'$\sigma_{{BAO}}$ ({} $h^{{-1}}Mpc$)'.format(self.sigma_r_BAO*self.h0))
            plt.xlabel("Distance [$h^{-1}$Mpc]")
            plt.title('Distribution of the distance of rim galaxies \n with respect to their centers')
            plt.legend(loc=2)
            plt.grid(True)
            plt.ylim(-0.002, 1.5*max(normed_rim_BAO_r_dist))
            plt.savefig("diagnostics/BAO_distribution.pdf")
            plt.close()
            fig = plt.figure(figsize=(6,6))
            plt.plot(bin_centers, normed_rim_BAO_r_dist, label='Distribution in the mock')
            plt.plot([self.r_BAO*self.h0, self.r_BAO*self.h0], [0, normed_rim_BAO_r_dist.max()], 
                     label=r'r$_{{BAO}}$ ({} $h^{{-1}}Mpc$)'.format(self.r_BAO*self.h0))
            plt.plot([(self.r_BAO-1.2*self.sigma_r_BAO)*self.h0, (self.r_BAO+1.2*self.sigma_r_BAO)*self.h0],
                     [normed_rim_BAO_r_dist.max()/np.e, normed_rim_BAO_r_dist.max()/np.e], 
                      label=r'$\sigma_{{BAO}}$ ({} $h^{{-1}}Mpc$)'.format(self.sigma_r_BAO*self.h0))
            plt.xlabel("Distance [$h^{-1}$Mpc]")
            plt.title('Distribution of the distance of rim galaxies \n with respect to their centers')
            plt.legend(loc=2)
            plt.grid(True)
            plt.ylim(-0.002, 1.5*max(normed_rim_BAO_r_dist))
            plt.savefig("diagnostics/BAO_distribution_colored.pdf")
            plt.close()
        return

    def generate_flat_galaxies(self, is_random=False):
        if is_random is False:
            r_flat, theta_flat, phi_flat = self.generate_galaxies(self.n_flat)
        else:
            print("generating flat galaxies for random catalog")
            r_flat, theta_flat, phi_flat = self.generate_galaxies(self.n_rand)
        flat_galaxies = {}
        for i, _ in enumerate(r_flat):
            flat_galaxies["flat_{}".format(i)] = galaxy(theta=theta_flat[i], phi=phi_flat[i], r=r_flat[i], TYPE=2, name="flat_{}".format(i))
        self.catalog.flats = flat_galaxies
        return r_flat, theta_flat, phi_flat

    def generate_clumps(self, add_clumps_to_rims = False, diagnostics=False):
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
            num_flat_clumps   = self.nr_clump - num_center_clumps

        diagnostics = True
        if diagnostics or self.diagnostics:
            print("Number of seed galaxies for clumping: \n"\
                  "- Seed flat galaxies: {} \n"\
                  "- Seed (center+rim) galaxies: {}".format(num_flat_clumps, num_center_clumps))
            
        # randomly choose indices from centers for the clumps
        # generate the clump positions with respect to their origin (a center galaxu)
        try:
            if self.frac_c2r is None:
                galaxy_selection  = np.random.choice(np.arange(0, total_num_centers, 1), size=num_center_clumps)
                center_clumps_idx = galaxy_selection[galaxy_selection<self.n_center]
                rim_clumps_idx    = galaxy_selection[galaxy_selection>=self.n_center]-self.n_center
            else:
                rim_clumps_cnt    = int(num_center_clumps / ( 1 + self.frac_c2r ))
                center_clumps_cnt = num_center_clumps - rim_clumps_cnt
                center_clumps_idx = np.random.choice(np.arange(0, self.n_center, 1), size=center_clumps_cnt)
                rim_clumps_idx    = np.random.choice(np.arange(0, self.n_center*self.n_rim, 1), size=rim_clumps_cnt)
        except ValueError:
            rim_clumps_idx = []
            center_clumps_idx = []
            
        if self.diagnostics or diagnostics:
            print("- Seed center galaxies: {} \n"\
                  "- Seed rim galaxies   : {} \n".format(len(center_clumps_idx), len(rim_clumps_idx), num_flat_clumps))

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
            # we move the calculated distances by 1Mpc to have a realistic clustering
            clump_r          = ( self.r_scale ** self.gamma ) * ( np.random.pareto( self.gamma, n_clump_to_inject ) )
            if self.diagnostics or diagnostics:
                center_clump_r_dist, _    = np.histogram(clump_r, bins=100, range=(0, 100))
                self.center_clump_r_dist += center_clump_r_dist
            clump_phi        = np.random.uniform(0., 360., n_clump_to_inject)
            clump_theta      = np.arccos(np.random.uniform(-1., 1., n_clump_to_inject))*RAD2DEG-90.
            for j in range(n_clump_to_inject):
                # due to random distribution, sometimes we ggenerate locations very far from the
                # seed galaxy. These distance galaxies cannot be considered as part of the
                # cluster, so they are not added. the upper limit here is defined assuming
                # superclusters are of size roughly 200Mpc
                if clump_r[j] > 200.:
                    # these are the ones not to be added
                    continue
                else:
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
                                                                                   parent="cen_{}".format(idx), TYPE=3,
                                                                                   name="cenClump_{}_-1_{}".format(idx, j))
            # replace the child list of the seed galaxy.
            # it is safe, the new childs are appended to the existing list of childs
            curr_seed_galaxy.childs = curr_seed_galaxy_childs

        # generate the clumps around the rim galaxies (if they will be added)
        for idx in rim_clumps_idx:
            # calculate the index of the center galaxy of the corresponding rim galaxy
            seed_idx1 = int((idx-self.n_center) % self.n_center)
            # calculate the index of the rim galaxy with respect to the center galaxy
            seed_idx2 = int((idx-self.n_center) / self.n_center)
            # in some cases, the following part fails. patching for now but will take a look later.
            try:
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
                # we move the calculated distances by mu to have a realistic clustering
                clump_r          = ( self.r_scake ** self.gamma ) * ( np.random.pareto( self.gamma, n_clump_to_inject ) )
                if self.diagnostics or diagnostics:
                    rim_clump_r_dist, _    = np.histogram(clump_r, bins=100, range=(0, 100))
                    self.rim_clump_r_dist += rim_clump_r_dist
                clump_phi        = np.random.uniform(0., 360., n_clump_to_inject)
                clump_theta      = np.arccos(np.random.uniform(-1., 1., n_clump_to_inject))*RAD2DEG-90.
                for j in range(n_clump_to_inject):
                    # due to random distribution, sometimes we ggenerate locations very far from the
                    # seed galaxy. These distance galaxies cannot be considered as part of the
                    # cluster, so they are not added. the upper limit here is defined assuming
                    # superclusters are of size roughly 200Mpc
                    if clump_r[j] > 200.:
                        continue
                    else:
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
                            clumps_center["rimClump_{}_{}_{}".format(seed_idx1, seed_idx2, j)] = galaxy(theta=curr_clump[1], phi=curr_clump[2],
                                                                                                        r=curr_clump[0],
                                                                                                        parent="rim_{}_{}".format(seed_idx1, seed_idx2),
                                                                                                        TYPE=4, name="rimClump_{}_{}_{}".format(seed_idx1, seed_idx2, j))
                # replace the child list of the seed galaxy.
                # it is safe, the new childs are appended to the existing list of childs
                curr_seed_galaxy.childs = curr_seed_galaxy_childs
            except KeyError:
                if self.diagnostics or diagnostics:
                    print("a problem is detected for retrieving current rim. moving on...")

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
            # we move the calculated distances by 1Mpc to have a realistic clustering
            clump_r          = ( self.r_scale ** self.gamma )* ( np.random.pareto( self.gamma-1, n_clump_to_inject ) )
            if self.diagnostics or diagnostics:
                flat_clump_r_dist, _    = np.histogram(clump_r, bins=100, range=(0, 100))
                self.flat_clump_r_dist += flat_clump_r_dist
            clump_phi        = np.random.uniform(0., 360., n_clump_to_inject)
            clump_theta      = np.arccos(np.random.uniform(-1., 1., n_clump_to_inject))*RAD2DEG-90.
            for j in range(n_clump_to_inject):
                # due to random distribution, sometimes we ggenerate locations very far from the
                # seed galaxy. These distance galaxies cannot be considered as part of the
                # cluster, so they are not added. the upper limit here is defined assuming
                # superclusters are of size roughly 200Mpc
                if clump_r[j] > 200.:
                    continue
                else:
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
                                                                                  parent="flat_{}".format(idx), TYPE=5,
                                                                                  name="flatClump_{}_-1_{}".format(idx, j))
            # replace the child list of the seed galaxy.
            # it is safe, the new childs are appended to the existing list of childs
            curr_seed_galaxy.childs = curr_seed_galaxy_childs
            # update the catalog object as all the center clumps are generated
        self.catalog.clumps_flat = clumps_flat

        # if diagnostics is enables, plot the distribution of the distances of clumps with respect to their seeds
        if self.diagnostics or diagnostics:
            values     = np.linspace(0, 100, 101)
            mid_values = (values[1:]+values[:-1])*0.5*self.h0
            self.check_diagnostics_directory()
            fig = plt.figure(figsize=(6,6))
            plt.plot(mid_values, self.flat_clump_r_dist/np.sum(self.flat_clump_r_dist), label="flat seeded clumps", color='black')
            plt.plot(mid_values, self.center_clump_r_dist/np.sum(self.center_clump_r_dist), '--', label="center seeded clumps", color='black')
            plt.plot(mid_values, self.rim_clump_r_dist/np.sum(self.rim_clump_r_dist), '-.', label="rim seeded clumps", color='black')
            plt.yscale("log")
            plt.title('Distribution of the distances of clumping galaxies \n with respect to their seeds')
            plt.xlabel("Distance [$h^{-1}$Mpc]")
            plt.legend()
            plt.grid(True)
            plt.savefig("diagnostics/clump_r_distributions.pdf")
            plt.close()
            fig = plt.figure(figsize=(6,6))
            plt.plot(mid_values, self.flat_clump_r_dist/np.sum(self.flat_clump_r_dist), label="flat seeded clumps")
            plt.plot(mid_values, self.center_clump_r_dist/np.sum(self.center_clump_r_dist), label="center seeded clumps")
            plt.plot(mid_values, self.rim_clump_r_dist/np.sum(self.rim_clump_r_dist), label="rim seeded clumps")
            plt.yscale("log")
            plt.title('Distribution of the distances of clumping galaxies \n with respect to their seeds')
            plt.xlabel("Distance [$h^{-1}$Mpc]")
            plt.legend()
            plt.grid(True)
            plt.savefig("diagnostics/clump_r_distributions_colored.pdf")
            plt.close()
            average_flat_clump_r   = np.average(mid_values, weights=self.flat_clump_r_dist/np.sum(self.flat_clump_r_dist))
            average_center_clump_r = np.average(mid_values, weights=self.center_clump_r_dist/np.sum(self.center_clump_r_dist))
            average_rim_clump_r    = np.average(mid_values, weights=self.rim_clump_r_dist/np.sum(self.rim_clump_r_dist))
            print("Average distance for flat clumps  : {:.2f} h^-1 Mpc [{:.2f} Mpc]".format(average_flat_clump_r, average_flat_clump_r/self.h0))
            print("Average distance for center clumps: {:.2f} h^-1 Mpc [{:.2f} Mpc]".format(average_center_clump_r, average_center_clump_r/self.h0))
            print("Average distance for rim clumps   : {:.2f} h^-1 Mpc [{:.2f} Mpc]".format(average_rim_clump_r, average_rim_clump_r/self.h0))
            
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
        config['r_scale']     = self.r_scale
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

    def calculate_n_bins(self, nEvents, rule="Rice"):
        if rule == "Sturge":
            return int(np.ceil(1+3.322*np.log10(nEvents)))
        elif rule == "Rice":
            return int(np.ceil(np.power(nEvents, 1./3.)*2))
        else:
            print("Problem with the choice of binning method. Falling back to default")
            return int(np.ceil(np.power(nEvents, 1./3.)*2))
        
    def write_to_fits(self, filename=None, is_random=False, save_extended=False):
        if filename is None:
            if not is_random:
                filename = self.fname_mock
            else:
                filename = self.fname_random
        rs, ras, decs, types, parents, names = self.catalog.flatten()
        LUT = self.generate_LUT_r2z()
        zs  = LUT(rs)
        ws  = np.ones(len(zs))
        # acceptance test here. should it be done for individual types?
        if self.acceptance:
            accepted_indices = self.check_z_acceptance(zs)
            zs      = np.asarray(zs)[accepted_indices]
            ras     = np.asarray(ras)[accepted_indices]
            decs    = np.asarray(decs)[accepted_indices]
            types   = np.asarray(types)[accepted_indices]
            parents = np.asarray(parents)[accepted_indices]
            names   = np.asarray(names)[accepted_indices]
            ws      = np.asarray(ws)[accepted_indices]
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
        header['r_scale']     = self.r_scale
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
        if save_extended:
            col6 = fits.Column(name="parent", array=parents, format='A16')
            col7 = fits.Column(name="name", array=names, format='A16')
            cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7])
        else:
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
