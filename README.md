# sportran

Code to compute thermal conductivity through cepstral analysis of heat flux time series, as described in papers:

 - (cepstral analysis) [Ercole, Marcolongo, Baroni, *Sci. Rep.* **7**, 15835 (2017)](https://doi.org/10.1038/s41598-017-15843-2)
 - (multicomponent systems) [Bertossa, Grasselli, Ercole, Baroni, *Phys. Rev. Lett.* **122**, 255901 (2019)](https://doi.org/10.1103/PhysRevLett.122.255901) ([arXiv](https://arxiv.org/abs/1808.03341))
 - (review) [Baroni, Bertossa, Ercole, Grasselli, Marcolongo, *Handbook of Materials Modeling* (2018)](https://doi.org/10.1007/978-3-319-50257-1_12-1) ([arXiv](https://arxiv.org/abs/1802.08006))

Developed by Loris Ercole, Riccardo Bertossa, Sebastiano Bisacchi under the supervision of prof. Stefano Baroni

**Acknowledgment**  The development of this software is part of the scientific program of the EU MaX Centre of Excellence for Supercomputing Applications (Grant No. 676598, 824143) and has been partly funded through it.

### Usage
There is a [**GUI**](README_GUI.md) that you can try after installing everything.

The code can be used as a **library**, for example in a Jupyter notebook.
In the [`examples`](examples/) folder you can find some examples.

Alternatively, you can run the code `analysis.py` from the **command line** without any installation procedure.
It can execute most of the cepstral analysis routines, returning the results in a series of data files and PDF plots.
See the [`examples/example_commandline_NaCl`](examples/example_commandline_NaCl/) folder and the help (`python analysis.py --help`) for more information.

### Requirements
#### only cepstral analysis (command line interface)
 - numpy
 - scipy
 - matplotlib

#### graphical user interface
 - all the requirements of cepstral analysis
 - tkinter
 - pillow (you may need upgrade it)
 - tk-html-widgets
 - markdown2

### Installation
  You can simply pip-install sportran downloading it from PyPI with `pip install sportran`.

  Alternatively:

  1. Clone this repository: `git clone https://github.com/lorisercole/sportran.git`
  2. Install the package with pip (dependencies will be automatically downloaded). For example:
```
cd sportran
pip install .
```
  You are all set! You can check that the installation is working by trying to run the command `sportran-analysis`.

  The Graphical User Interface will be installed and can be started with the command `sportran-gui`.

### Issues
  You are strongly encouraged to report any issue on the [official](https://github.com/lorisercole/sportran/issues) GitHub issues page.
