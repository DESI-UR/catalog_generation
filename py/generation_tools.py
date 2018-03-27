"""
Nodules to generate random and mock catalogs

Contributors: Dylan and Tolga
"""


import numpy as np
import healpy as hp
import math, random
import matplotlib as mpl
#mpl.use("Agg")
import matplotlib.pyplot as plt
import configparser
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM, z_at_value
import astropy.units as u
import sys
from ROOT import TFile, TTree, gROOT
import os

NSIDE = -1
fbool = False

DEG2RAD = np.pi/180.0
RAD2DEG = 1./DEG2RAD

class GeneralTools():
    
    def __init__(self, configFile, diagnostics=False):
        """
        COMMENTS HERE
        """
        self.diagnostics = diagnostics
        self.config_file = configFile
        self.getConfig()
        self.get_template()
        self.get_mask()


    def setDiagnostics(self, diagnostics):
        self.diagnostics = diagnostics
        
    def getConfig(self):
        self.config     = configparser.RawConfigParser()
        self.config.read(self.config_file)

        # some of the numbers to use in the generation
        self.num_random = int(self.config.get("N Catalogs", "nrandom"))
        self.num_mock   = int(self.config.get("N Catalogs", "nmock"))
        
        # read in the cosmology parameters to generate the mocks
        self.H0       = float(self.config.get("Cosmology", "H_0"))
        self.omega_m0 = float(self.config.get("Cosmology", "omega_m"))
        self.cosmo        = FlatLambdaCDM(self.H0, self.omega_m0)
        
        # output filenames
        self.fname_random = self.config.get("File Names", "randout")
        self.fname_mock   = self.config.get("File Names", "mockout")
        self.datafile     = self.config.get('File Names','datafile')
        self.ang_mask     = self.config.get('File Names','angmask')
        
        # the following variable defined the coordinate system, 0: theta, phi, 1: RA, DEC
        self.coordinates  = int(self.config.get("Data Format", "coordinates"))

        # BAO related parameters
        self.r_BAO          = float(self.config.get('Gen Params','r_bao'))
        self.sigma_r_BAO    = float(self.config.get('Gen Params','sigr_bao'))
        self.gamma          = float(self.config.get('Gen Params','gamma'))
        self.r_0            = float(self.config.get('Gen Params','r_0'))
        self.n_rnd          = int(self.config.get('Gen Params','n_rand'))
        self.n_center       = int(self.config.get('Gen Params','n_center'))
        self.n_rim          = int(self.config.get('Gen Params','n_rim'))
        self.n_flat         = int(self.config.get('Gen Params','n_flat'))
        self.nr_clump       = int(self.config.get('Gen Params','nr_clump'))
        self.n_clump        = int(self.config.get('Gen Params','n_clump'))
        self.n_clump_center = int(self.config.get('Gen Params','n_centerclump'))
        try:
            self.z_lo       = float(self.config.get("Gen Params", "z_lo"))
            self.z_hi       = float(self.config.get("Gen Params", "z_hi"))
        except:
            print("using predefined z bounds")
            self.z_lo       = 0.4
            self.z_hi       = 0.7

    # read the template to use for the mock and random catalogs
    def get_template(self):
        self.hdulist      = fits.open(self.datafile)
        self.data_z       = self.hdulist[1].data['Z']
        template_cut = np.array([(self.data_z[i] < self.z_hi) and (self.data_z[i] > self.z_lo) for i in range(len(self.data_z))])
        self.template_z   = self.data_z[template_cut]
        self.template_r   = [r.value for r in self.cosmo.comoving_distance(self.template_z)]
        self.template_r_len = len(self.template_r)
        self.template_r_min = np.amin(self.template_r)
        self.template_r_max = np.amax(self.template_r)

    # read the mask/completeness for the output
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

    # generate r values uniformly distributed with (r_min, r_max)
    def generate_uniform_r(self, num_obs, diagnostics=False):
        flat_r = np.random.uniform(self.template_r_min, self.template_r_max, num_obs)
        # if diagnostics enabled, plot the distribution
        if diagnostics or self.diagnostics:
            self.check_diagnostics_directory()
            plt.hist(flat_r)
            plt.xlabel("r [Mpc]")
            plt.title("Flat r distribution")
            plt.savefig("diagnostics/r_distribution.pdf")        
        return flat_r
    
    # generate r values based on histogrammed data r distribution
    # this function is used to generate r distribution depending on the
    # template provided
    def generate_r(self, nobs, diagnostics=False):
        num_obs = 0
        rlist   = []
        n, r_edges = np.histogram(self.template_r, bins=int(np.sqrt(self.template_r_len)))
        n = n.astype(float)
        weights    = n/np.sum(n)
        r_med = (r_edges[:-1] + r_edges[1:])/2.
        while num_obs < nobs:
            curr_r = np.random.choice(r_med, p=weights)
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
    
    def generate_uniform_angular_position(self, nobs, diagnostics=False):
        num_obs = 0
        pixlist = []
        weights = self.completeness/np.sum(self.completeness)
        pix_idx = np.arange(self.completeness_len)
        while num_obs < nobs:
            curr_pix = np.random.choice(pix_idx, p=weights)
            pixlist.append(curr_pix)
            num_obs += 1
        angles = hp.pix2ang(ipix=np.array(pixlist), nside=self.nside, lonlat=True)
        # plot the distribution for diagnostics
        if diagnostics or self.diagnostics:
            # completeness plot
            plt.clf()
            pixels = np.zeros([12*self.nside**2])
            pixels[pixlist] = 1
            hp.mollview(pixels.astype(int), title="Completeness")
            self.check_diagnostics_directory()
            plt.savefig("diagnostics/generated_completeness.pdf")
            # theta distribution
            plt.clf()
            print(angles[0])
            plt.hist(angles[0], bins=int(np.sqrt(self.completeness_len)))
            plt.title(r"$\theta$ distribution")
            plt.xlabel(r'$\theta$ [deg]')
            plt.savefig("diagnostics/generated_theta_dist.pdf")
            # theta distribution
            plt.clf()
            plt.hist(angles[1], bins=int(np.sqrt(self.completeness_len)))
            plt.title(r"$\phi$ distribution")
            plt.xlabel(r'$\phi$ [deg]')
            plt.savefig("diagnostics/generated_phi_dist.pdf")
        return angles
            
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

    # function to convert spherical coordinates to cartesian coordinates
    def toCartesianVector(self, r, theta, phi):
        unit_vector = hp.pix2vec(self.nside, hp.ang2pix(self.nside, theta, phi, lonlat=True))
        curr_vector = r * np.array(unit_vector)
        return curr_vector
        
    # function to convert cartesian coordinates to spherical coordinates
    def fromCartesianVector(self, vec):
        r           = np.sqrt( np.sum(np.array(vec)**2) )
        unit_vector = hp.vec2ang(vec, lonlat=True)
        return r, unit_vector[0], unit_vector[1]
        
    # function to add two spherical vectors
    def addVectors(self, vec1, vec2):
        vec1_cartesian  = self.toCartesianVector(*vec1)
        vec2_cartesian  = self.toCartesianVector(*vec2)
        sumvec          = self.fromCartesianVector(vec1_cartesian+vec2_cartesian)
        # this part is necessary to get the added vector to the right quadrant
        # it is not intended for general use but specific to this project
        return sumvec

    #
    def generate_rim(self, r, theta, phi):
        rim_rs     = []
        rim_thetas = []
        rim_phis   = []
        len_n_center = len(r)
        for i in range(len_n_center):
            curr_rim_cnt = 0
            while curr_rim_cnt < self.n_rim:
                curr_theta = np.random.uniform(0., 360., 1)
                curr_phi   = np.random.uniform(-90., 90., 1)
                curr_r     = self.generate_gaussian(self.r_BAO, self.sigma_r_BAO, 1)
                curr_rim   = self.addVectors(np.array([r[i], theta[i], phi[i]]),
                                             np.array([curr_r[0], curr_theta[0], curr_phi[0]]))
                pixel      = hp.ang2pix(self.nside, curr_rim[1], curr_rim[2], lonlat=True)
                # this is the acceptance.
                if self.completeness[pixel] > 0:
                    rim_rs.append(curr_rim[0])
                    rim_thetas.append(curr_rim[1])
                    rim_phis.append(curr_rim[2])
                    curr_rim_cnt += 1
        rim_rs = np.array(rim_rs).flatten()
        rim_thetas = np.array(rim_thetas).flatten()
        rim_phis = np.array(rim_phis).flatten()
        return rim_rs, rim_thetas, rim_phis

    def generate_clumps(self, r_centers, theta_centers, phi_centers, diagnostics=False):
        # generate a list to choose among the centers and flats
        galaxy_selection  = np.random.randint(0, 2, self.nr_clump)
        num_center_clumps = len(galaxy_selection[galaxy_selection==0])
        num_flat_clumps   = len(galaxy_selection[galaxy_selection==1])
        # generate flat galaxies (will be returned and be added to the mocks later)
        r_flat, theta_flat, phi_flat = self.generate_galaxies(self.n_flat)
        # randomly choose indices from centers and flat for the clumps
        rand_center_idx = np.random.randint(0, len(r_centers), num_center_clumps)
        rand_flat_idx   = np.random.randint(0, len(r_flat),    num_flat_clumps)
        selected_r_centers     = r_centers[rand_center_idx]
        selected_theta_centers = theta_centers[rand_center_idx]
        selected_phi_centers   = phi_centers[rand_center_idx]
        selected_r_flats       = r_flat[rand_flat_idx]
        selected_theta_flats   = theta_flat[rand_flat_idx]
        selected_phi_flats     = phi_flat[rand_flat_idx]
        # generate the clump positions with respect to their origin (a center galaxu)
        center_clumps  = []
        for i in range(num_center_clumps):
            curr_r_center     = selected_r_centers[i]
            curr_theta_center = selected_theta_centers[i]
            curr_phi_center   = selected_phi_centers[i]
            clump_r        = (self.r_0**self.gamma * (np.random.pareto(self.gamma-1, self.n_clump_center)))
            clump_theta    = np.random.uniform(0., 360., self.n_clump_center)
            clump_phi      = np.random.uniform(-90., 90., self.n_clump_center)
            for j in range(self.n_clump_center):
                center_clumps.append(self.addVectors(np.array([curr_r_center, curr_theta_center, curr_phi_center]),
                                                     np.array([clump_r[j], clump_theta[j], clump_phi[j]])))

        # generate the clump positions with respect to their origin (a flat galaxu)
        flat_clumps    = []
        for i in range(num_flat_clumps):
            curr_r_flat     = selected_r_flats[i]
            curr_theta_flat = selected_theta_flats[i]
            curr_phi_flat   = selected_phi_flats[i]
            clump_r        = (self.r_0**self.gamma * (np.random.pareto(self.gamma-1, self.n_clump)))
            clump_theta    = np.random.uniform(0., 360., self.n_clump)
            clump_phi      = np.random.uniform(-90., 90., self.n_clump)
            for j in range(self.n_clump):
                flat_clumps.append(self.addVectors(np.array([curr_r_flat, curr_theta_flat, curr_phi_flat]),
                                                   np.array([clump_r[j], clump_theta[j], clump_phi[j]])))
        
        return (np.array(center_clumps)).transpose(), (np.array(flat_clumps)).transpose(), np.array([r_flat, theta_flat, phi_flat])
        
    def write_to_fits(self, col1, col2, col3, col4, filename):
        col_defs = [['phi' 'theta', 'z', 'weight'], ['ra', 'dec', 'z', 'weight']]
        print(self.coordinates)
        use_col_defs = col_defs[self.coordinates]
        # We also write the output in fits format
        if os.path.isfile(filename):
            print("a file with the designated name already exists... please remove the file first")
            return
        col1 = fits.Column(name=use_col_defs[0], array=col1, format='f8')
        col2 = fits.Column(name=use_col_defs[1], array=col2, format='f8')
        col3 = fits.Column(name=use_col_defs[2], array=col3,  format='f8')
        col4 = fits.Column(name=use_col_defs[3], array=col4,  format='f8')
        cols = fits.ColDefs([col1, col2, col3, col4])
        hdu  = fits.BinTableHDU.from_columns(cols)
        hdu.writeto(filename)
