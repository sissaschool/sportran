# -*- coding: utf-8 -*-

from __future__ import absolute_import
from setuptools import setup, find_packages
from os import path
import sys
import json

GUI_dependencies = [
        'future-fstrings>=0.1.0',
        'markdown2>=2.0.0',
        'pillow>=5.4.0',
        'uncertainties>=2.4'
]

if sys.version_info[0] == 2:
    INSTALL_GUI = False
    sys.stderr.write("************* GUI does not support Python 2: it will not be installed. That's a pity! *************")
elif sys.version_info[0] == 3:
    INSTALL_GUI = True


if __name__ == '__main__':
    THIS_FOLDER = path.split(path.abspath(__file__))[0]

    with open(path.join(THIS_FOLDER, 'README.md'), 'r') as fh:
        long_description = fh.read()

    with open(path.join(THIS_FOLDER, 'setup.json'), 'r') as info:
        SETUP_JSON = json.load(info)

    if INSTALL_GUI:
        exclude_packages = []
        SETUP_JSON['install_requires'] += GUI_dependencies
        SETUP_JSON['entry_points']['console_scripts'] += ['thermocepstrum-gui = thermocepstrum_gui.main:run']
        with open(path.join(THIS_FOLDER, 'thermocepstrum/metadata.json'), 'w') as fh:
           json.dump(SETUP_JSON, fh)
    else:
        exclude_packages = ['thermocepstrum_gui*']

    setup(include_package_data=True,
          packages=find_packages(exclude=exclude_packages),
          long_description=long_description,
          long_description_content_type='text/markdown',
          **SETUP_JSON)
