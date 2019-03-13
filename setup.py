#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function

#
# Standard imports
#
import glob
import os
import sys

#
# setuptools' sdist command ignores MANIFEST.in
#
from distutils.command.sdist import sdist as DistutilsSdist
from setuptools import setup, find_packages
from setuptools.command.install import install as InstallCommand

class Install(InstallCommand):
    """ Customized setuptools install command which uses pip."""
    def run(self, *args, **kwargs):
        from pip._internal import main
        main(['install', '.'])
        InstallCommand.run(self, *args, **kwargs)

#
# Begin setup
#
setup_keywords = dict()
#
# THESE SETTINGS NEED TO BE CHANGED FOR EVERY PRODUCT.
#
setup_keywords['name'] = 'ParaMock'
setup_keywords['description'] = 'Parametric Mock Catalog Generation for BAO analysis'
setup_keywords['author'] = 'Tolga Yapici'
setup_keywords['author_email'] = 'tyapici@ur.rochester.edu'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/DESI-UR/catalog_generation'
setup_keywords['version'] = '0.1.0'
#
# END OF SETTINGS THAT NEED TO BE CHANGED.
#
#setup_keywords['version'] = get_version(setup_keywords['name'])

# Set other keywords for the setup function.  These are automated, & should
# be left alone unless you are an expert.
#
# Treat everything in bin/ except *.rst as a script to be installed.
#
if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
        if not os.path.basename(fname).endswith('.rst')]

setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['requires'] = ['Python (>3.5)', 'healpy (>=1.11.0)', 'numpy (>=1.13.1)', 'configparser (>=3.5)', \
                              'astropy (>=1.2.1)', 'scipy (>=0.19.1)', 'matplotlib (>=2.0.0)']
setup_keywords['install_requires'] = ['healpy>=1.11.0', 'numpy>=1.13.1', 'configparser>=3.5', \
                                      'astropy>=1.2.1', 'scipy>=0.19.1', 'matplotlib>=2.0.0']
setup_keywords['zip_safe'] = False
setup_keywords['use_2to3'] = True
setup_keywords['packages'] = find_packages('py')
setup_keywords['package_dir'] = {'':'py'}
setup_keywords['cmdclass'] = {'install': Install,}

#setup_keywords['cmdclass'] = {'version': DesiVersion, 'test': DesiTest, 'sdist': DistutilsSdist}
#setup_keywords['test_suite']='{name}.test.test_suite'.format(**setup_keywords)

#
# Autogenerate command-line scripts.
#
# setup_keywords['entry_points'] = {'console_scripts':['desiInstall = desiutil.install.main:main']}
#
# Add internal data directories.
#
setup_keywords['package_data'] = {'ParaMock': ['data/*',]}
#
# Run setup command.
#
setup(**setup_keywords)
