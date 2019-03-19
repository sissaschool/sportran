# thermocepstrum
Code to compute thermal conductivity through cepstral analysis of heat flux time series, as described in papers:
 - Ercole, Marcolongo, Baroni, Sci. Rep. 7, 15835 (2017), https://doi.org/10.1038/s41598-017-15843-2
 - for the multi-component analysis:  Baroni, Bertossa, Ercole, Grasselli, Marcolongo, https://arxiv.org/abs/1802.08006

**Acknowledgement**  The development of this software is part of the scientific program of the EU MaX Centre of Excellence for Supercomputing Applications (Grant No. 676598) and has been partly funded through it.

### Usage
The code can be used as a **library**, for example in a Jupyter notebook. 
In the "example" folder you can find some examples. 

Alternatively, you can run the code "analysis.py" from the **command line**.
It can execute most of the cepstral analysis routines, writing the result in a series of data files and PDF plots.
See help (python analysis.py --help) for more informations.

### Requirements
 - Python 2.7
 - Numpy
 - Scipy
 - Matplotlib

### Installation
  1. Clone this repository ```git clone https://github.com/lorisercole/thermocepstrum.git```
  2. Select the branch that you want, for example ```cd thermocepstrum; git checkout beta```
  3. Run the script ```./pack.sh```. You will find in the ```dist``` directory the python package.
  4. Install with pip: for example ```pip install --user dist/thermocepstrum-0.1-py2-none-any.whl```. Dependencies will be automatically downloaded
  5. You are done! Check that the installation is working by trying the command ```thermocepstrum-analysis```. If it is not working, maybe you have to add the folder where you python package are installed (on my system is ```~/.local/bin```) to your ```PATH```. Remember to read the examples in the ```examples``` directory!