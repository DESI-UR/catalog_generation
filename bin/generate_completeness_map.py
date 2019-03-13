
import healpy as hp
import numpy as np
import sys
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="This script generates a completeness map based on "+\
                                 "the given RA and Dec limits")
parser.add_argument("-r", "--ra",     type=float, help="RA limits",  nargs=2, required=True)
parser.add_argument("-d", "--dec",    type=float, help="Dec limits", nargs=2, required=True)
parser.add_argument("-o", "--output", type=str,   help="output filename", nargs=1, required=True)
args   = parser.parse_args()
"""

theta_min = float(sys.argv[1])
theta_max = float(sys.argv[2])
phi_min = float(sys.argv[3])
phi_max = float(sys.argv[4])

try:
    nside = int(sys.argv[5])
except:
    nside = 1024

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

hp.mollview(pixels, cbar=False, title='Completeness')
plt.savefig("completeness.pdf")
plt.show()

np.savez("example_c.npz", CMPLTNSS=pixels)
"""
