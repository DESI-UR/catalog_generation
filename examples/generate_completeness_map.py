
import healpy as hp
import numpy as np
import sys
import matplotlib.pyplot as plt
import argparse

RAD2DEG    = 180./np.pi
DEG2RAD    = 1/RAD2DEG

parser = argparse.ArgumentParser(description="This script generates a completeness map based on "+\
                                 "the given RA and Dec limits")
parser.add_argument("-r", "--ra",     type=float, nargs=2, required=True, help="RA limits  (range: [0, 360])")
parser.add_argument("-d", "--dec",    type=float, nargs=2, required=True, help="Dec limits (range: [-90, 90])")
parser.add_argument("-n", "--nside",  type=int,   help="nside for the generated cocmpleteness", default=1024)
parser.add_argument("-o", "--output", type=str,   help="output filename", required=True)
parser.add_argument("-p", "--plot",   type=int,   help="plot the resukts in mollview", default=0)
args   = parser.parse_args()

ra_min = np.min(args.ra)
ra_max = np.max(args.ra)

dec_min = np.min(args.dec)
dec_max = np.max(args.dec)

nside      = args.nside
resolution = hp.nside2resol(nside)*RAD2DEG

decs   = np.arange(dec_min, dec_max, resolution)
ras    = np.arange(ra_min,  ra_max, resolution)
pixels = np.zeros(hp.nside2npix(nside))

for ra in ras:
    for dec in decs:
        pixels[hp.ang2pix(nside, 90.-ra, -dec, lonlat=True)] = 1
pixels = np.array(pixels)

if args.plot:
    hp.mollview(pixels, cbar=False, title='Completeness')
    plt.savefig("completeness.pdf")
    plt.show()
    
np.savez(args.output, CMPLTNSS=pixels)
