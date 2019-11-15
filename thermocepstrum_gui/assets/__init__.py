from pkg_resources import resource_stream, resource_string
import json

with resource_stream('thermocepstrum', 'metadata.json') as JS:
    METADATA = json.load(JS)

with resource_stream(__name__, 'languages.json') as JS:
    LANGUAGES = json.load(JS)

ICON = resource_string(__name__, 'icon.gif')

README_MD = resource_string('thermocepstrum', 'README.md')
README_GUI_MD = resource_string('thermocepstrum_gui', 'README_GUI.md')
