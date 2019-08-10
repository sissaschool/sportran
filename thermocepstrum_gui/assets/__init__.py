
from pkg_resources import resource_stream, resource_string
import json

with resource_stream('thermocepstrum','setup.json') as JS:
    METADATA=json.load(JS)

ICON=resource_string(__name__,'icon.gif')