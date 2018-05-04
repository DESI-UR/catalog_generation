
# coding: utf-8

# In[1]:


import sys
sys.path.append("../py/")
sys.path.append("..")

import matplotlib as mpl
import matplotlib.pyplot as plt

import numpy as np

from astropy.io import fits

hdus = fits.open(sys.argv[1])
centers       = hdus[1].data[hdus[1].data['TYPE']==0]
rims          = hdus[1].data[hdus[1].data['TYPE']==1]
center_clumps = hdus[1].data[hdus[1].data['TYPE']==2]
flat_clumps   = hdus[1].data[hdus[1].data['TYPE']==3]
flats         = hdus[1].data[hdus[1].data['TYPE']==4]

print("Number of center       galaxies: {}".format(len(centers)))
print("Number of rim          galaxies: {}".format(len(rims)))
print("Number of flat         galaxies: {}".format(len(flats)))
print("Number of center clump galaxies: {}".format(len(center_clumps)))
print("Number of flat clump   galaxies: {}".format(len(flat_clumps)))
print("Number of              galaxies: {}".format(len(centers)+len(rims)+len(flats)+len(center_clumps)+len(flat_clumps)))

plt.plot(centers['z'], 'o')
plt.show()

plt.figure(1, figsize=(16,8))

try:
    plt.subplot(241)
    plt.hist(centers['z'], bins=int(np.sqrt(len(centers['z']))), alpha=0.3, normed=True, label="BAO Centers")
    plt.hist(rims['z'], bins=int(np.sqrt(len(rims['z']))), alpha=0.3, normed=True, label="BAO Rims")
    plt.legend(loc=2)
    plt.xlabel("z")
except:
    print("skipping rims")

plt.subplot(242)
plt.hist(centers['z'], bins=int(np.sqrt(len(centers['z']))), alpha=0.3, normed=True, label="BAO Centers")
plt.hist(center_clumps['z'], bins=int(np.sqrt(len(center_clumps['z']))), alpha=0.3, normed=True, label="Center clumps")
plt.legend(loc=2)
plt.xlabel("z")

plt.subplot(243)
plt.hist(centers['z'], bins=int(np.sqrt(len(centers['z']))), alpha=0.3, normed=True, label="BAO Centers")
plt.hist(flat_clumps['z'], bins=int(np.sqrt(len(flat_clumps['z']))), alpha=0.3, normed=True, label="Flat clumps")
plt.legend(loc=2)
plt.xlabel("z")

plt.subplot(244)
plt.hist(centers['z'], bins=int(np.sqrt(len(centers['z']))), alpha=0.3, normed=True, label="BAO Centers")
plt.hist(flats['z'], bins=int(np.sqrt(len(flats['z']))), alpha=0.3, normed=True, label="Flats")
plt.legend(loc=2)
plt.xlabel("z")

ax = plt.subplot(245)
ax.hist(centers['ra'], alpha=0.3, normed=True)
try:
    ax.hist(rims['ra'], alpha=0.3, normed=True)
except:
    print("skipping rims")
ax.hist(center_clumps['ra'], alpha=0.3, normed=True)
ax.hist(flat_clumps['ra'], alpha=0.3, normed=True)
ax.hist(flats['ra'], alpha=0.3, normed=True)
ax.set_xlabel(r'RA [deg]')
ax = plt.subplot(246)
ax.hist(centers['dec'], alpha=0.7, label="center galaxies", normed=True)
try:
    ax.hist(rims['dec'], alpha=0.3, label="rim galaxies", normed=True)
except:
    print("skipping rims")
ax.hist(center_clumps['dec'], alpha=0.3, label="center clumps", normed=True)
ax.hist(flat_clumps['dec'], alpha=0.3, label="flat clumps", normed=True)
ax.hist(flats['dec'], alpha=0.3, label="flat galaxies", normed=True)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width , box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax.set_xlabel(r'DEC [deg]')
ax.text(600, 0.0040, r'{}'.format(sys.argv[1].replace(".fits", "")), horizontalalignment='right', fontweight='bold')
ax.text(600, 0.0035, r'N$_{{center}}$: {}'.format(len(centers)), horizontalalignment='right')
ax.text(600, 0.0030, r'N$_{{rim}}$: {}'.format(len(rims)), horizontalalignment='right')
ax.text(600, 0.0025, r'N$_{{flat}}$: {}'.format(len(flats)), horizontalalignment='right')
ax.text(600, 0.0020, r'N$_{{center clump}}$: {}'.format(len(center_clumps)), horizontalalignment='right')
ax.text(600, 0.0015, r'N$_{{flat clump}}$: {}'.format(len(flat_clumps)), horizontalalignment='right')
ax.text(600, 0.0010, r'N$_{{total}}$: {}'.format(len(centers)+len(rims)+len(flats)+len(center_clumps)+len(flat_clumps)), horizontalalignment='right')
        
filename = "distributions_" + sys.argv[1].replace(".fits", ".pdf")
plt.savefig(filename, bbox_inches='tight')
