"""
Nodules to generate random and mock catalogs

Contributors: Dylan and Tolga
"""


import numpy as np
import healpy as hp
import math, random
import matplotlib as mpl
mpl.use("Agg")
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

class GeneralTools():
    
    def __init__(self, configFile):
        """
        COMMENTS HERE
        """
        self.config_file = configFile
        self.getConfig()
        self.get_template()
        self.get_mask()
        
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
        if diagnostics:
            plt.clf()
            hp.mollview(self.completeness, title="Completeness")
            if not os.path.isdir("diagnostics"):
                os.makedirs("diagnostics")
            plt.savefig("diagnostics/completeness.pdf")

    # generate r values uniformly distributed with (r_min, r_max)
    def generate_uniform_r(self, num_obs, diagnostics=False):
        flat_r = np.random.uniform(self.template_r_min, self.template_r_max, num_obs)
        # if diagnostics enabled, plot the distribution
        if diagnostics:
            plt.hist(flat_r)
            plt.xlabel("r [Mpc]")
            if not os.path.isdir("diagnostics"):
                os.makedirs("diagnostics")
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
        if diagnostics:
            plt.hist(rlist)
            plt.xlabel("r [Mpc]")
            if not os.path.isdir("diagnostics"):
                os.makedirs("diagnostics")
            plt.savefig("diagnostics/r_distribution.pdf")
        return rlist

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
        if diagnostics:
            # completeness plot
            plt.clf()
            pixels = np.zeros([12*self.nside**2])
            pixels[pixlist] = 1
            hp.mollview(pixels.astype(int), title="Completeness")
            if not os.path.isdir("diagnostics"):
                os.makedirs("diagnostics")
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
        x = r * math.sin(theta) * math.cos(phi)
        y = r * math.sin(theta) * math.sin(phi)
        z = r * math.cos(theta)
        return np.array([x, y, z])

    # function to convert cartesian coordinates to spherical coordinates
    def fromCartesianVector(self, vec):
        r     = np.sqrt( np.sum(np.array(vec)**2) )
        theta = np.degrees(np.arctan(vec[1]/vec[0]))
        phi   = np.degrees(np.arctan(vec[2]/vec[0]))
        return np.array([r, theta, phi])

    # function to add two spherical vectors
    def addVectors(self, vec1, vec2):
        vec1_cartesian = self.toCartesianVector(*vec1)
        vec2_cartesian = self.toCartesianVector(*vec2)
        return self.fromCartesianVector(vec1_cartesian+vec2_cartesian)

    #
    def generate_rim(self, r, theta, phi):
        rim_rs     = []
        rim_thetas = []
        rim_phis   = []
        len_n_center = len(r)
        for i in range(len_n_center):
            curr_theta = np.random.uniform(0., 360., self.n_rim)
            curr_phi   = np.random.uniform(0., 180., self.n_rim)
            curr_r     = self.generate_gaussian(self.r_BAO, self.sigma_r_BAO, self.n_rim)
            for j in range(self.n_rim):
                curr_rim   = self.addVectors(np.array([r[i], theta[i], phi[i]]),
                                             np.array([curr_r[j], curr_theta[j], curr_phi[j]]))
                rim_rs.append(curr_rs)
                rim_thetas.append(curr_thetas)
                rim_phis.append(curr_phis)    
        rim_rs = np.array(rim_rs).flatten()
        rim_thetas = np.array(rim_thetas).flatten()
        rim_phis = np.array(rim_phis).flatten()
        # need to add acceptance for the rim galaxies
        return rim_rs, rim_thetas, rim_phis

    def generate_clumps(self, r_centers, theta_centers, phi_centers, diagnostics=False):
        # generate flat galaxies (will be returned and be added to the mocks later)
        r_flat, theta_flat, phi_flat = self.generate_galaxies(self.n_flat)
        # randomly choose indices from centers and flat for the clumps
        rand_center_idx = np.random.randint(0, len(r_centers), self.n_clump_center)
        rand_flat_idx   = np.random.randint(0, len(r_flat),    self.n_clump)
        # generate the clump positions with respect to their origin (a center galaxu)
        clump__r       = (self.r_0**self.gamma * (np.random.pareto(self.gamma-1, self.n_clump_center)))
        clump_theta    = np.random.uniform(0., 360., self.n_clump_center)
        clump_phi      = np.random.uniform(0., 180., self.n_clump_center)
        center_clumps  = self.addVectors(np.array([r_centers[rand_center_idx], theta_centers[rand_center_idx], phi_centers[rand_center_idx]]),
                                              np.array([clump_r, clump_theta, clump_phi]))
        # generate the clump positions with respect to their origin (a flat galaxu)
        clump_r        = (self.r_0**self.gamma * (np.random.pareto(self.gamma-1, self.n_clump)))
        clump_theta    = np.random.uniform(0., 360., self.n_clump)
        clump_phi      = np.random.uniform(0., 180., self.n_clump)
        flat_clumps    = self.addVectors(np.array([r_flat[rand_flat_idx], theta_flat[rand_flat_idx], phi_flat[rand_flat_idx]]),
                                         np.array([clump_r, clump_theta, clump_phi]))
        # merge two arrays
        clump_rs     = np.append(center_clumps[0], flat_clumps[0])
        clump_thetas = np.append(center_clumps[1], flat_clumps[1])
        clump_phis   = np.append(center_clumps[2], flat_clumps[2])
        #
        return clump_rs, clump_thetas, clump_phis, f_flat, theta_flat, phi_flat
        
    #generate random catalog with data distribution in r,phi,theta
    def generate_random(self, catalog_no=0):
        try:
            from ROOT import Galaxy
        except:
            self.genROOTObject()
            from ROOT import Galaxy
        tfilename   = self.fname_random+"_"+str(catalog_no)+".root"
        output      = TFile(tfilename, "RECREATE")
        output_tree = TTree("data_pol", "data_pol")
        galaxy      = Galaxy()
        if self.coordinates == 0:
            output_tree.Branch("position_pol", galaxy, "phi/F:theta/F:z/F:weight/F")
        else:
            output_tree.Branch("position_pol", galaxy, "ra/F:dec/F:z/F:weight/F")
        rs, thetas, phis = self.generate_galaxies(self.n_rnd)
        zs               = [z_at_value(self.cosmo.comoving_distance, r*u.Mpc) for r in rs]
        weights          = np.array([1.]*len(rs))  # THIS SEEMS WRONG
        galaxy  = Galaxy()
        ras     = phis*(180./np.pi)
        decs    = ((np.pi/2.)-thetas)*(180./np.pi)
        for i in range(self.n_rnd):
            galaxy.z = zs[i]
            if self.coordinates == 0:
                galaxy.phi = phis[i]
                galaxy.theta = thetas[i]
                galaxy.w = weights[i]
            else:
                galaxy.ra = ras[i]
                galaxy.dec = decs[i]
                galaxy.weight = weights[i]
            output_tree.Fill()
        output.Write()
        output.Close()
        fits_filename = self.fname_random+"_"+str(ncat)+".fits"
        if self.coordinates == 0:
            self.write_to_fits(phis, thetas, zs, weights, fits_filename)
        else:
            self.write_to_fits(ras, decs, zs, weights, fits_filename)

    def generate_mock(self, catalog_no=0):
        try:
            from ROOT import Galaxy
        except:
            self.genROOTObject()
            from ROOT import Galaxy
        tfilename   = self.fname_random+"_"+str(catalog_no)+".root"
        output      = TFile(tfilename, "RECREATE")
        output_tree = TTree("data_pol", "data_pol")
        galaxy      = Galaxy()
        if self.coordinates == 0:
            output_tree.Branch("position_pol", galaxy, "phi/F:theta/F:z/F:weight/F")
        else:
            output_tree.Branch("position_pol", galaxy, "ra/F:dec/F:z/F:weight/F")
        r_centers, theta_centers, phi_centers = self.generate_galaxies(self.n_center)
        # TODO:
        # add flats
        # add rims
        # add clumps
            
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

