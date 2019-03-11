#!/bin/bash

set -x

cp -v LICENSE thermocepstrum/
cp -v README.md thermocepstrum/ 
python3 -m pip install --user --upgrade setuptools wheel

cd thermocepstrum
python3 setup.py sdist bdist_wheel

python3 -m pip install --user --upgrade twine
python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*


