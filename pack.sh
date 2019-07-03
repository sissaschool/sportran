#!/bin/bash


python -m pip install --user --upgrade setuptools wheel
python -m pip install --user --upgrade twine

python setup.py sdist bdist_wheel

echo
echo "==========================="
echo "| upload file with command |"
echo "==========================="
echo 'python3 -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*'
