import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

if __name__ == '__main__':
    setuptools.setup(
        include_package_data=True,
        name="thermocepstrum",
        version="0.1.2",
        author="Loris Ercole, Riccardo Bertossa",
        author_email="loris.ercole@epfl.ch",
        description="Cepstral Data Analysis of current time series for Green-Kubo transport coefficients",
        license="GPL 3",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/lorisercole/thermocepstrum",
        keywords="cepstral data analysis thermal conductivity transport coefficients physics green-kubo",
        install_requires=['numpy>=1.9.0','scipy>=0.17.0','matplotlib>=2.0.0'],
        python_requires='>=2.6',
        packages=setuptools.find_packages(), # exclude=['docs','tests*'] ...
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
        ],
        package_data={
             'thermocepstrum.utils' : ["plot_style.mplstyle"],
            },
        entry_points={
            'console_scripts': [
                'thermocepstrum-analysis = thermocepstrum.analysis:main'
            ],
        },
    )

# https://packaging.python.org/guides/distributing-packages-using-setuptools/#setup-args
# https://setuptools.readthedocs.io/en/latest/setuptools.html#including-data-files

