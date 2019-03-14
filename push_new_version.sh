#!/bin/bash

git clone https://${GH_TOKEN}@github.com/DESI-UR/catalog_generation.git -b development catgen-dev
cd catgen-dev
python py/paramock/versioning.py
git config --global user.email "travis@travis-ci.org"
git config --global user.name "Travis CI"
git commit -m "after the successful build, updating the version number" py/paramock/_version.py
git push
