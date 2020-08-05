#!/bin/sh

set -e

. $(dirname $(dirname $0))/functions.sh



python setup.py bdist_wheel
WHEEL=$(python setup.py --name)-$(python setup.py --version)-py2.py3-none-any.whl
pip install -U "dist/$WHEEL[train]"
python -m coverage run -p -m unittest discover -vv
