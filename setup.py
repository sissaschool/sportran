import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="thermocepstrum",
    version="0.1",
    author="Loris Ercole, Riccardo Bertossa",
    author_email="loris.ercole@epfl.ch",
    description="Cepstral Data Analysis for Green-Kubo transport coefficients",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lorisercole/thermocepstrum",
    keywords="cepstral data analysis thermal conductivity transport coefficients physics green-kubo",
    install_requires=['numpy','scipy','matplotlib'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "Topic :: Scientific/Engineering :: Physics",
        "Environment :: Console"
    ]
)
# scripts=['thermocepstrum']  --> analysis.sh
# OR BETTER: entry_points
# packages=setuptools.find_packages(exclude=['docs','tests*']),  --> exclude these subpackages from release
# py_modules = ['aaaaa']  --> aaaaaa.py is not part of the package
# python_requires='>=2.6'
# data_files = ....
# https://packaging.python.org/guides/distributing-packages-using-setuptools/#setup-args
# https://setuptools.readthedocs.io/en/latest/setuptools.html#including-data-files

