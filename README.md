# thermocepstrum
Code to compute thermal conductivity through cepstral analysis of heat flux time series, as described in papers:
 - [Ercole, Marcolongo, Baroni, Sci. Rep. 7, 15835 (2017)](https://doi.org/10.1038/s41598-017-15843-2),
 - (multicomponent systems) [Bertossa, Grasselli, Ercole, Baroni Phys. Rev. Lett. 122, 255901 (2019)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.255901),
 - [Baroni, Bertossa, Ercole, Grasselli, Marcolongo](https://arxiv.org/abs/1802.08006),
 
**Acknowledgment**  The development of this software is part of the scientific program of the EU MaX Centre of Excellence for Supercomputing Applications (Grant No. 676598) and has been partly funded through it.

### Usage
There is a [**GUI**](README_GUI.md) that you can try after pip-installing everything.

The code can be used as a **library**, for example in a Jupyter notebook.
In the [`examples`](examples/) folder you can find some examples.

Alternatively, you can run the code `analysis.py` from the **command line** without any installation procedure.
It can execute most of the cepstral analysis routines, returning the results in a series of data files and PDF plots.
See the [`examples/example_commandline_NaCl`](examples/example_commandline_NaCl/) folder and the help (`python analysis.py --help`) for more information.

### Requirements
#### only cepstral analysis (command line interface)
 - Python 2.7 or 3.x
 - Numpy
 - Scipy
 - Matplotlib
#### graphical user interface
 - all the requirements of cepstral analysis
 - python 3 only
 - tkinter
 - future-fstrings
 - pillow (you may need to run `pip install pillow --upgrade` )
 - tk-html-widgets
 - markdown2

### Installation
We provide two python packages. The core library and the command line interface are in the subdirectory `thermocepstrum`, while the GUI is in the `thermocepstrum_gui` subdirectory. Remember that the GUI requires python 3, while the library and the command line interface don't. That's why we have two packages.
  You can simply pip-install thermocepstrum with `pip install thermocepstrum` from the pypi repository.
  
  Alternatively:

  1. Clone this repository: `git clone https://github.com/lorisercole/thermocepstrum.git`
  2. cd into the new directory created by git
   * system where `setup.sh` does work:
     1. run the bash script `setup.sh` and you are done (skip next steps).
   * system where `setup.sh` does not work:
     1. copy the files `README.md README_GUI.md LICENSE.txt` in the source folders `thermocepstrum/thermocepstrum` and `thermocepstrum_gui/thermocepstrum_gui/`
     2. `cd thermocepstrum; pip install .`
     3. `cd thermocepstrum_gui; pip install`

The copy is needed  because during python packaging it is not possible to include files from parent directories, and we want to have only a single version of the READMEs in the project tree.

You are done! You can check that the installation is working by trying to run the command `thermocepstrum-analysis`.
The graphical user interface is started with the command `thermocepstrum-gui` (only if you installed it).

### Issues
You are strongly encouraged to report any issue on the [official](https://github.com/lorisercole/thermocepstrum/issues) GitHub issues page.

