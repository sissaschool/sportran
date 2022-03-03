# -*- coding: utf-8 -*-

from pkg_resources import resource_stream, resource_string
import json

with resource_stream('sportran', 'metadata.json') as JS:
    METADATA = json.load(JS)
dev_state = ''
if 'classifiers' in METADATA:
    dev_state = [x for x in METADATA['classifiers'] if 'Development Status ::' in x]
    if len(dev_state) > 0:
        dev_state = dev_state[0].split('::')[1]
    else:
        dev_state = ''

with resource_stream(__name__, 'languages.json') as JS:
    LANGUAGES = json.load(JS)

ICON = resource_string(__name__, 'icon.gif')

README_MD = resource_string('sportran', 'README.md')
README_GUI_MD = resource_string('sportran_gui', 'README_GUI.md')
