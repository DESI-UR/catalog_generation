#!/usr/bin/env python

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="This script generates a completeness map based on "+\
                                 "the given RA and Dec limits")
parser.add_argument("-r", "--ra",
                    type=float,
                    help="RA limits",
                    nargs=2,
                    required=True)
parser.add_argument("-d", "--dec",
                    type=float,
                    help="Dec limits",
                    nargs=2,
                    required=True)
parser.add_argument("-n", "--nside",
                    type=int,
                    help="Nside for the healpix map",
                    required=False,
                    default=1024,
                    choices=[16, 32, 64, 128, 256, 512, 1024]) 
parser.add_argument("-o", "--output",
                    type=str,
                    help="output filename",
                    required=True)
parser.add_argument("-p", "--plot",
                    type=lambda x: (str(x).lower() in ['true','1', 'yes']),
                    required=False,
                    default=False)
args   = parser.parse_args()

theta_min = float(min(args.dec))
theta_max = float(max(args.dec))
phi_min   = float(min(args.ra))
phi_max   = float(max(args.ra))
nside     = args.nside

RAD2DEG    = 180./np.pi
DEG2RAD    = 1/RAD2DEG
resolution = hp.nside2resol(nside)*RAD2DEG
    
thetas = np.arange(theta_min, theta_max, resolution)
phis   = np.arange(phi_min, phi_max, resolution)
pixels = np.zeros(hp.nside2npix(nside))

for theta in thetas:
    for phi in phis:
        pixels[hp.ang2pix(nside, theta, phi, lonlat=True)] = 1
pixels = np.array(pixels)

if args.plot:
    hp.mollview(pixels, cbar=False, title='Completeness')
    hp.graticule()
    plt.savefig("completeness.pdf")
    plt.show()

np.savez(args.output, CMPLTNSS=pixels)
