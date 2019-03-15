import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="thermocepstrum",
    version="0.1",
    author="Loris Ercole, Riccardo Bertossa",
    author_email="lercole@sissa.it",
    description="Cepstral analysis of heat current time series in single and multicomponent substances",
    license="GPL 3",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lorisercole/thermocepstrum",
    install_requires=['numpy','scipy','matplotlib'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Physics",
        "Environment :: Console"
    ],
    package_data={
	 '' : ["grafici_belli.mplstyle"],
 	},
    entry_points={
        'console_scripts': [
            'thermocepstrum-analysis = thermocepstrum.analysis:main'
        ],
    },
)
