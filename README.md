# thermocepstrum
Code to compute thermal conductivity through cepstral analysis of heat flux time series, as described in papers:
 - Ercole, Marcolongo, Baroni, Sci. Rep. 7, 15835 (2017), https://doi.org/10.1038/s41598-017-15843-2
 - for the multi-component analysis:  Baroni, Bertossa, Ercole, Grasselli, Marcolongo, https://arxiv.org/abs/1802.08006

### usage

The file analysis.py is the main program, that can run from the command line with python 2.7.
It executes all the parts of cepstral analysis producing the final result and a number of figures in a pdf output.
See help (python analysis.py --help) for more informations.

The code can be used as a library, in a jupyter notebook for example, as shown in the example directory.

