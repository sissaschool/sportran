# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from os import path
import sys
import json

if __name__ == '__main__':
    THIS_FOLDER = path.split(path.abspath(__file__))[0]

    with open(path.join(THIS_FOLDER, 'README.md'), 'r') as fh:
        long_description = fh.read()

    with open(path.join(THIS_FOLDER, 'setup.json'), 'r') as info:
        SETUP_JSON = json.load(info)

    with open(path.join(THIS_FOLDER, 'sportran/metadata.json'), 'w') as fh:
        json.dump(SETUP_JSON, fh)

    setup(include_package_data=True,
          packages=find_packages(),
          long_description=long_description,
          long_description_content_type='text/markdown',
          **SETUP_JSON)
