# Examples

In this folder you can find some practical examples using *thermocepstrum* to compute the thermal conductivity of a system.

* the `data` folder contains the trajectories, obtained from classical NVE molecular dynamic simulations, of a silica and a molten NaCl system, that are analysed in the examples.

#### Jupyter notebook examples
These examples show how to use the thermocepstrum package in a Python script, step by step, to analyse a time series and perform cepstral analysis. Some of the tools and plot functions of the code are presented.
* `example_cepstrum_singlecomp_silica.ipynb`: analysis of solid amorphous silica, a one-component system.
* `example_cepstrum_doublecomp_NaCl.ipynb`: analysis of molten NaCl, a two-component system.

#### Command line examples
These examples show how to use the command line program `analysis.py` to perform cepstral analysis in a straightforward. This is an expedient tool to analyse a time series using predefined parameters. Results are produced in the form of numerous text/binary files and pdf plots.
* `example_commandline_NaCl`: analysis of molten NaCl, a two-component system.
