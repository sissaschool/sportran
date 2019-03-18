#!/bin/bash

set -x

python -m pip install --user --upgrade setuptools wheel

cd thermocepstrum
python setup.py sdist bdist_wheel

# to test on this machin
# python -m pip install dist/dokr-0.1-py3-none-any.whl

python -m pip install --user --upgrade twine
python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*


