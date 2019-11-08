from setuptools import setup, find_packages
import json


try:
    with open('thermocepstrum_gui/README.md', 'r') as fh:
        long_description = fh.read()
except:
    raise RuntimeError("Maybe you must run 'copy_readmes_in_packages.sh' in the parent directory to prepare the file tree?")

if __name__ == '__main__':
    with open('thermocepstrum_gui/setup.json', 'r') as info:
        kwargs = json.load(info)

    setup(include_package_data=True,
          packages=find_packages(),   # exclude=['docs','tests*'] ...
          setup_requires=['reentry'],
          reentry_register=True,
          long_description=long_description,
          long_description_content_type='text/markdown',
          **kwargs)

# https://packaging.python.org/guides/distributing-packages-using-setuptools/#setup-args
# https://setuptools.readthedocs.io/en/latest/setuptools.html#including-data-files
