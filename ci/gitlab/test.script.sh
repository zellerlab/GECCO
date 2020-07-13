#!/bin/sh

set -e

. $(dirname $(dirname $0))/functions.sh


python setup.py bdist_wheel
pip install -U dist/*.whl
python -m coverage run -p -m unittest discover -vv
