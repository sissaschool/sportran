.. SporTran documentation master file, created by
   sphinx-quickstart on Sun Mar 22 22:31:28 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SporTran's documentation!
====================================


Introduction
============
.. toctree::
   README.md
   README_GUI.md
   CONTRIBUTING.md

Examples
========

Jupyter notebook examples
-------------------------

These examples show how to use the sportran package in a Python script, step by step, to analyse a time series and perform cepstral analysis. Some of the tools and plot functions of the code are presented.

.. toctree::
   example_cepstrum_singlecomp_silica
   example_cepstrum_doublecomp_NaCl
   example_input_formats

Command line examples
---------------------

These examples show how to use the command line program sportran-analysis (analysis.py) to perform cepstral analysis in a straightforward way. This is an expedient tool to analyse a time series using predefined parameters. Results are produced in the form of several text/binary files and pdf plots.
You can find the examples in the `repository <https://github.com/sissaschool/sportran/tree/develop/examples>`_. 

.. toctree::
   EXAMPLE_SILICA.md
   EXAMPLE_NACL.md

LAMMPS input files for energy current generation
------------------------------------------------

.. toctree::
   LAMMPS_SILICA.md

API
===
.. toctree::
   modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
