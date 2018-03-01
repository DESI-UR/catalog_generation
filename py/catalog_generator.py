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

class CatalogGenerator():
    
    def __init__(self, configFile):
        """
        COMMENTS HERE
        """
        self.z_lo        = 0.4
        self.z_hi        = 0.7
        self.H0          = 75.0
        self.omega_m0    = 0.30
        self.config_file = configFile

    def getConfig(self):
        self.config     = configparser.RawConfigParser()
        self.config.read(self.config_file)
        # some of the numbers to use in the generation
        self.num_random = int(self.config.get("N Catalogs", "nrandom"))
        self.num_mock   = int(self.config.get("N Catalogs", "nmock"))
        # output filenames
        self.fname_random = self.config.get("File Names", "randout")
        self.fname_mock   = self.config.get("File Names", "mockout")
        # the following variable defined the coordinate system, 0: theta, phi, 1: RA, DEC
        self.coordinates  = int(self.config.get("Data Format", "coordinates"))
        #
        self.datafile     = self.config.get('File Names','datafile')
        self.ang_mask     = self.config.get('File Names','angmask')
        # cosmology to use
        self.H0           = float(self.config.get('Cosmology','h_0'))
        self.omega_m0     = float(self.config.get('Cosmology','omega_m'))
        self.cosmo        = FlatLambdaCDM(self.H0, self.omega_m0)
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
        #
        try:
            self.z_lo       = float(self.config.get("Gen Params", "z_lo"))
            self.z_hi       = float(self.config.get("Gen Params", "z_hi"))
        except:
            print("using predefined z bounds")
        
    def getLists(self):
        self.hdulist      = fits.open(self.datafile)
        self.zlist        = self.hdulist[1].data['Z']
        self.zcut         = np.array([(self.zlist[i] < self.z_hi) and (self.zlist[i] > self.z_lo) for i in range(len(self.zlist))])
        self.zlist_cut    = self.zlist[self.zcut]
        self.mask         = np.load(self.ang_mask)
        self.completeness = self.mask.f.CMPLTNSS
        self.nside        = hp.npix2nside(len(self.completeness))
        self.radlist      = [r.value for r in self.cosmo.comoving_distance(self.zlist_cut)]
    
    def genROOTObject(self):
        if self.coordinates == 0:
            gROOT.ProcessLine("struct Galaxy {\
            float phi;\
            float theta;\
            float z;\
            float w;\
            };" );
        else:
            gROOT.ProcessLine("struct Galaxy {\
            float ra;\
            float dec;\
            float z;\
            float weight;\
            };" );
            
    def generate_randoms(self):
        try:
            from ROOT import Galaxy
        except:
            self.genROOTObject()
            from ROOT import Galaxy
        for curr_random in range(self.num_random):
            self.random_catgen(curr_random)

    def generate_mocks(self):
        try:
            from ROOT import Galaxy
        except:
            self.genROOTObject()
            from ROOT import Galaxy
        for curr_mock in range(self.num_mock):
            self.mock_catgen(curr_mock)

    # generate r values based on histogrammed data r distribution
    def rgen(self, histo, nobs):
        num_obs = 0
        dlist = []
        num_bins = len(histo[0])
        r_min    = (histo[1])[0]
        r_max    = (histo[1])[num_bins]
        seln  = np.amax(histo[0])
        while num_obs < nobs:
            rsel = random.uniform(r_min,r_max)
            nbin = int((rsel-r_min)*num_bins/r_max)
            nsel = random.uniform(0.,seln) 
            nlim = (histo[0])[nbin] 
            # what is the logic here? it does not make sense.
            if (nsel < nlim):
                dlist.append(rsel)
                num_obs = num_obs + 1
        return dlist

    #generate r values based on uniform spatial distribution
    def urgen(self, rmin, rmax, nobs):
        umin = float(rmin)/rmax
        ulist = np.random.uniform(umin*umin*umin, 1, nobs)
        dlist = [np.power(u,1./3) for u in ulist]
        return dlist

    #generate Gaussian distributed r values
    def gaussrgen(self, mu, sig, nobs):
        llist = []
        for nob in range(nobs):
            leng = random.gauss(mu,sig)
            llist.append(leng)
        return llist

    #put list of r,phi,theta tuples through r acceptance based on data r histogram
    def racc(self, histo, inlist):
        outrlist = []
        outplist = []
        outtlist = []
        nbins = len(histo[0])
        rmin  = (histo[1])[0]
        rmax  = (histo[1])[nbins]
        inlist2 = [[],[],[]]
        for k in range(len(inlist[0])):
            if ((inlist[0])[k] >= rmin) and ((inlist[0])[k] <= rmax):
                (inlist2[0]).append((inlist[0])[k])
                (inlist2[1]).append((inlist[1])[k])
                (inlist2[2]).append((inlist[2])[k])
        seln  = np.amax(histo[0])
        histp = plt.hist(inlist2[0],bins=nbins)
        histq = [float((histo[0])[k])/(histp[0])[k] for k in range(nbins)]
        iind  = np.argmax(histq)
        factr = histq[iind]
        histr = [binval/factr for binval in histq]
        pmin  = np.amin(histp[0])
        for i in range(len(inlist2[0])):
            r = (inlist2[0])[i]
            nsel = random.uniform(0., 1.)
            nbin = int((r-rmin)*nbins/(rmax-rmin))
            if nbin == nbins:
                nbin = nbins - 1
            nlim = histr[nbin]
            if (nsel < nlim):
                outrlist.append((inlist2[0])[i])
                outplist.append((inlist2[1])[i])
                outtlist.append((inlist2[2])[i])
        return [outrlist,outplist,outtlist]

    #randomly generates unit 3-vector
    def generate_unit_vector(self):
        v = [random.gauss(0,1) for i in range(3)]
        factor = 1.0 / math.sqrt(sum(x*x for x in v))
        return np.array([x * factor for x in v])

    #randomly generates nobs unit 3-vectors
    def vecgen(self, nobs):
        vlist = []
        for nob in range(nobs):
            nuvec = self.generate_unit_vector()
            vlist.append(nuvec)
        return vlist

    #checks whether unit vector is within data angular mask
    def isinmap(self, vect):
        pixid = hp.vec2pix(self.nside, vect[0], vect[1], vect[2])
        return (self.completeness[pixid] > 0.)

    #randomly generates nobs unit 3-vectors within data angular mask
    def vecgen2(self, nobs):
        nob = 0
        vlist = []
        while nob < nobs:
            nuvec = self.generate_unit_vector()
            if self.isinmap(nuvec):
                vlist.append(nuvec)
                nob = nob + 1
        return vlist

    #converts list of cartesian 3-vectors to lists of r's,phi's,theta's
    def cart2pol(self, inlist):
        outrlist = []
        outplist = []
        outtlist = []
        for v in inlist:
            r = math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
            ev = [e/r for e in v]
            ang = hp.rotator.vec2dir(ev)
            if (ang[1] >= 0):
                phi = ang[1]
            else:
                phi = 2*np.pi + ang[1]
            theta = ang[0]
            outrlist.append(r)
            outplist.append(phi)
            outtlist.append(theta)
        return [outrlist,outplist,outtlist]

    #puts list of 3-vectors through angular mask acceptance
    def angacc(self, inlist):
        outlist = []
        for v in inlist:
            if self.isinmap(v):
                outlist.append(v)
        return outlist

    #randomly generates nobs 3-vectors according to Gaussian r distribution
    def vgaussgen(self, cen, mu, sig, nobs):
        rl = self.gaussrgen(mu, sig, nobs)
        vl = self.vecgen(nobs)
        vl2 = [np.add(cen,np.multiply(rl[j],vl[j])) for j in range(nobs)]
        return vl2

    #randomly generate nobs 3-vectors according to inverse power r distribution (with lower cutoff at r=1)
    def vinvgen(self, cen, gam, nobs):
        rl = np.random.pareto(self.gamma-1,nobs)+1
        vl = self.vecgen(nobs)
        vl2 = [np.add(cen,np.multiply(rl[j],vl[j])) for j in range(nobs)]
        return vl2

    def generate_flat_galaxies(self):
        flatlistv = self.vecgen2(self.n_rnd)
        try:
            flatlistr = self.rgen(self.hist0, self.n_rnd)
        except:
            self.hist = plt.hist(self.radlist, bins=400)
            flatlistr = self.rgen(self.hist, self.n_rnd)
        flatlist  = [np.multiply(flatlistr[k], flatlistv[k]) for k in range(self.n_rnd)]
        flatlist2 = self.cart2pol(flatlist)
        warr      = np.array([1.]*len(flatlist2[0]))
        pharr     = np.array(flatlist2[1])
        tharr     = np.array(flatlist2[2])
        rlist3    = flatlist2[0]
        zlist3    = [z_at_value(self.cosmo.comoving_distance, r*u.Mpc, zmin=0.3, zmax=0.8) for r in rlist3]
        rlist3    = []
        zarr      = np.array(zlist3)
        zlist3    = []
        flatlist2 = []
    
    #generate random catalog with data distribution in r,phi,theta
    def random_catgen(self, ncat):
        try:
            from ROOT import Galaxy
        except:
            self.genROOTObject()
            from ROOT import Galaxy
        tfilename   = self.fname_random+"_"+str(ncat)+".root"
        output      = TFile(tfilename, "RECREATE")
        output_tree = TTree("data_pol", "data_pol")
        galaxy      = Galaxy()
        if self.coordinates == 0:
            output_tree.Branch("position_pol", galaxy, "phi/F:theta/F:z/F:weight/F")
        else:
            output_tree.Branch("position_pol", galaxy, "ra/F:dec/F:z/F:weight/F")
        flatlistv = self.vecgen2(self.n_rnd)
        try:
            flatlistr = self.rgen(self.hist0, self.n_rnd)
        except:
            self.hist0= plt.hist(self.radlist, bins=400)
            flatlistr = self.rgen(self.hist0, self.n_rnd)
        flatlist  = [np.multiply(flatlistr[k], flatlistv[k]) for k in range(self.n_rnd)]
        flatlist2 = self.cart2pol(flatlist)
        warr      = np.array([1.]*len(flatlist2[0]))
        pharr     = np.array(flatlist2[1])
        tharr     = np.array(flatlist2[2])
        rlist3    = flatlist2[0]
        zlist3    = [z_at_value(self.cosmo.comoving_distance, r*u.Mpc, zmin=0.3, zmax=0.8) for r in rlist3]
        rlist3    = []
        zarr      = np.array(zlist3)
        zlist3    = []
        flatlist2 = []
        galaxy    = Galaxy()
        raarr     = pharr*(180./np.pi)
        decarr    = ((np.pi/2.)-tharr)*(180./np.pi)
        for i in range(self.n_rnd):
            galaxy.z = zarr[i]
            if self.coordinates == 0:
                galaxy.phi = pharr[i]
                galaxy.theta = tharr[i]
                galaxy.w = warr[i]
            else:
                galaxy.ra = raarr[i]
                galaxy.dec = decarr[i]
                galaxy.weight = warr[i]
            output_tree.Fill()
        output.Write()
        output.Close()
        fits_filename = self.fname_random+"_"+str(ncat)+".fits"
        if self.coordinates == 0:
            self.write_to_fits(pharr, tharr, zarr, warr, fits_filename)
        else:
            self.write_to_fits(raarr, decarr, zarr, warr, fits_filename)
            
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

        
    #generate mock catalog with BAO signal and clumping, fitting data r,phi,theta distribution
    def mock_catgen(self, ncat):
        try:
            from ROOT import Galaxy
        except:
            self.genROOTObject()
            from ROOT import Galaxy
        tfilename   = self.fname_mock+"_"+str(ncat)+".root"
        output      = TFile(tfilename, "RECREATE")
        output_tree = TTree("data_pol", "data_pol")
        galaxy      = Galaxy()
        if self.coordinates == 0:
            output_tree.Branch("position_pol", galaxy, "phi/F:theta/F:z/F:weight/F")
        else:
            output_tree.Branch("position_pol", galaxy, "ra/F:dec/F:z/F:weight/F")
        cenlistv = self.vecgen2(self.n_center)
        try:
            cenlistr = self.rgen(self.hist0, self.n_center)
        except:
            self.hist0= plt.hist(self.radlist, bins=400)
            cenlistr = self.rgen(self.hist0, self.n_center)
        cenlist = [np.multiply(cenlistr[k], cenlistv[k]) for k in range(self.n_center)]
        cenlistv = []
        cenlistr = []
        rimlist = []
        for k in range(self.n_center):
            listv = self.vgaussgen(cenlist[k], self.r_BAO, self.sigma_r_BAO, self.n_rim)
            rimlist.extend(listv)
        rimlist2 = self.angacc(rimlist)
        rimlist = []
        flatlistv = self.vecgen2(self.n_flat)
        flatlistr = self.rgen(self.hist0, self.n_flat)
        flatlist = [np.multiply(flatlistr[k],flatlistv[k]) for k in range(self.n_flat)]
        flatlistv = []
        flatlistr = []
        normlist = flatlist + rimlist2
        flatlist = []
        rimlist2 = []
        indlist = random.sample(range(self.n_center+len(normlist)),self.nr_clump) # or nr_clump???
        clumplist = []
        for k in indlist:
            clumps = []
            if k < self.n_center:
                clumps = self.vinvgen(cenlist[k],self.gamma,self.n_clump)
            else:
                clumps = self.vinvgen(normlist[k-self.n_center],self.gamma,self.n_clump_center)
            clumplist.extend(clumps)
        clumplist2 = self.angacc(clumplist)
        clumplist = []
        totlist = cenlist + normlist + clumplist2
        cenlist = []
        normlist = []
        clumplist2 = []
        totlist2 = self.cart2pol(totlist)
        totlist = []
        totlist3 = self.racc(self.hist0,totlist2)
        totlist2 = []
        warr = np.array([1.]*len(totlist3[0]))
        pharr = np.array(totlist3[1])
        tharr = np.array(totlist3[2])
        rlist3 = totlist3[0]
        zlist3 = [z_at_value(self.cosmo.comoving_distance, r*u.Mpc, zmin=0.3, zmax=0.8) for r in rlist3]
        rlist3 = []
        zarr = np.array(zlist3)
        zlist3 = []
        totlist3 = []
        raarr = pharr*(180./np.pi)
        decarr = ((np.pi/2.)-tharr)*(180./np.pi)
        for i in range(self.n_flat):
            galaxy.z = zarr[i]
            if self.coordinates == 0:
                galaxy.phi = pharr[i]
                galaxy.theta = tharr[i]
                galaxy.w = warr[i]
            else:
                galaxy.ra = raarr[i]
                galaxy.dec = decarr[i]
                galaxy.weight = warr[i]
            output_tree.Fill()
        output.Write()
        output.Close()
        # We also write the output in fits format
        fits_filename = self.fname_mock+"_"+str(ncat)+".fits"
        if self.coordinates == 0:
            self.write_to_fits(pharr, tharr, zarr, warr, fits_filename)
        else:
            self.write_to_fits(raarr, decarr, zarr, warr, fits_filename)

