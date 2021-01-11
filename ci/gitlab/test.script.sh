#!/bin/sh

set -e

. $(dirname $(dirname $0))/functions.sh



python setup.py build_data --inplace bdist_wheel
python -m pip install --find-links=dist gecco[train]
python -m coverage run -p -m unittest discover -vv
