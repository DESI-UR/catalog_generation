import numpy as np
from astropy.io import fits
import glob
import itertools

filenames = glob.glob("random*fits")

ra  = []
dec = []
z   = []
w   = []
t   = []

for filename in filenames:
    hdus = fits.open(filename)
    hdr  = hdus[1].header
    data = hdus[1].data
    ra.append((data['ra']).tolist())
    dec.append((data['dec']).tolist())
    z.append((data['z']).tolist())
    w.append((data['weight']).tolist())
    t.append((data['TYPE']).tolist())
    
ra  = list(itertools.chain.from_iterable(ra))
dec = list(itertools.chain.from_iterable(dec))
z   = list(itertools.chain.from_iterable(z))
w   = list(itertools.chain.from_iterable(w))
t   = list(itertools.chain.from_iterable(t))

col1 = fits.Column(name="ra",  array=ra, format='f8')
col2 = fits.Column(name="dec", array=dec, format='f8')
col3 = fits.Column(name="z",   array=z, format='f8')
col4 = fits.Column(name="weight", array=w, format='f8')
col5 = fits.Column(name="TYPE", array=t, format='f8')
cols = fits.ColDefs([col1, col2, col3, col4, col5])
hdu  = fits.BinTableHDU.from_columns(cols, header=hdr)
hdu.writeto("random.fits")
