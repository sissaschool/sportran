#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module is the CLI of the library, that can be called after installing it from the command line
"""

import os
from sys import path, argv
import argparse
import numpy as np

try:
    import sportran as st
except ImportError:
    abs_path = os.path.abspath(__file__)
    tc_path = abs_path[:abs_path.rfind('/')]
    path.append(tc_path[:tc_path.rfind('/')])
    try:
        import sportran as st
    except ImportError:
        raise ImportError('Cannot locate sportran.')

from sportran.utils import log

log.set_method('bash')
log.append_method('file')
from sportran.plotter.cli import CLIPlotter

st.Current.set_plotter(CLIPlotter)
from sportran.plotter import plt   # this imports matplotlib.pyplot
from sportran.plotter import addPlotToPdf, PdfPages

np.set_printoptions(precision=8)


def main():
    """
    --------------------------------------------------------------------------------
      *** SPORTRAN ***  command line interface
    --------------------------------------------------------------------------------
    This script performs the cepstral analysis of a (heat) current.
    Results are written to stdout and a log file, and plots are saved in PDF format.

    INPUT FORMAT:
     - table  : a column-formatted text file, with a header in the same format of LAMMPS.
                The name of the LAMMPS compute can start with c_ and end with [#some_number], the code will recognize
                vectors, and will read automatically all the components.
     - dict   : a Numpy binary file containing a dictionary (e.g. obtained from the script i_o/read_lammps_log.py)
     - LAMMPS : a LAMMPS log file.
                In this case a --run-keyword  must be provided, that identifies the desired 'run' command. This keyword must equal to the comment line placed just before the desired 'run' command (see documentation of i_o/read_lammps_log.py for an example).

    Physical parameters (time step, temperature, volume, units) must be provided.
    The average temperature is computed if a column with the header (or a dictionary key) 'Temp' is found; otherwise you have to specify it.

    You must provide the key that identifies the main current ('-k KEY')
    You can also provide additional currents if your system is a multi-component fluid ('-j CURRENT2 -j CURRENT3'), or you want to decorrelate the main current with respect to them (see PRL).
    (Notice that the output is the same with any number of components. If you have a lots of components, note that you may want to use more than 3 independent processes -- see theory.)

    OUTPUT FILES:
      [output].logfile
          A log of the available information.
      [output].plots.pdf
          A PDF with all the plots generated.
    OUTPUT DATA files (can be text ".dat" or binary ".npy"):
      [output].psd
          freqs [THz], original periodogram, original log(periodogram)
      [output].cospectrum (if any)
          freqs [THz], full matrix cospectrum
      [output].resampled_psd
          freqs [THz], resampled periodogram, resampled log(periodogram)
      [output].cepstral
          cepstral coefficients ck, error(ck), L0(P*), err(L0(P*)), kappa(P*) [W/mK], err(kappa(P*)) [W/mK]
          the line number minus one is the number of cepstral coefficients used (P*).
      [output].cepstrumfiltered_psd
          freqs [THz], cepstrum-filtered periodogram, cepstrum-filtered log(periodogram)

    -------------------------
    Example:
      read and analyze "examples/data/Silica.dat" file. The energy-flux columns are called c_flux[1], c_flux[2], c_flux[3]

        ./analysis "examples/data/Silica.dat" --VOLUME 3130.431110818 --TEMPERATURE 1065.705630 -t 1.0 -k flux1 -u metal -r --FSTAR 28.0 -w 0.5 -o silica_test
    -------------------------
    """
    _epilog = """---
    Enjoy it!
    ---
    Developed by Loris Ercole, Riccardo Bertossa, Sebastiano Bisacchi, under the supervision of prof. Stefano Baroni at SISSA, Via Bonomea, 265 - 34136 Trieste ITALY.

    Please cite these references:
     - Ercole, Marcolongo, Baroni, Sci. Rep. 7, 15835 (2017), https://doi.org/10.1038/s41598-017-15843-2
     - Bertossa, Grasselli, Ercole, Baroni, Phys. Rev. Lett. 122, 255901 (2019), https://doi.org/10.1103/PhysRevLett.122.255901
     - Baroni, Bertossa, Ercole, Grasselli, Marcolongo, Handbook of Materials Modeling (2018), https://doi.org/10.1007/978-3-319-50257-1_12-1

    GitHub:    https://github.com/sissaschool/sportran
    Contact:   loris.ercole@epfl.ch, rbertoss@sissa.it

    Acknowledgment
    The development of this software is part of the scientific program of the EU MaX Centre of Excellence for Supercomputing Applications (Grant No. 676598, 824143) and has been partly funded through it.
    """

    # yapf: disable
    if '--list-currents' in argv:
        print('units and docstrings list for each current:')
        print('=================')
        print(st.current._list_of_currents_and_units(verbose=True))
        print('=================')
        print('currents and units implementation table')
        print(st.current.build_currents_units_table())
        return 0

    parser = argparse.ArgumentParser(description=main.__doc__, epilog=_epilog, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('inputfile', type=str,
            help='input file to read (default format: Table)')

    input_file_group = parser.add_argument_group('Input file format')
    input_file_group.add_argument('--input-format', default='table', type=str, choices=['table', 'dict', 'lammps'],
            help='Format of the input file. (default: table)')
    input_file_group.add_argument('-k', '--mainfluxkey', type=str, required=True,
            help='Name of the column keyword that identifies the first flux in the onsager matrix')
    input_file_group.add_argument('-j', '--add-currents', type=str, default=[], action='append',
            help='Additional current for multi-component fluids. (optional, repeat -j to add more currents)')
    input_file_group.add_argument('-N', '--nsteps', type=int, default=0,
            help='Number of steps to read. (optional, default: 0=all)')
    input_file_group.add_argument('-S', '--start-step', type=int, default=0,
            help='The first step to read. (optional, default: 0=first)')
    input_file_group.add_argument('--cindex', nargs='*', type=int,
            help='Column indexes of the main current to read (0,1,2,...). (optional, default: all)')
    input_file_group.add_argument('--sindex', nargs='*', type=int,
            help='Column indexes of a current to be substracted from the main current. (optional)')
    input_file_group.add_argument('--split', type=int, default=1,
            help='Build a time series with n*m independent processes (n is the number of processes of the original timeseries, m is the number provided with --split). The length of the new time series will be [original length]/m. (optional)')

    lammps_group = input_file_group.add_argument_group('LAMMPS format settings')
    lammps_group.add_argument('--run-keyword', type=str,
            help='Keyword that identifies the run to be read: a specific comment line placed just before the run command (only for "lammps" format)')
    lammps_group.add_argument('--structure', type=str,
            help='LAMMPS data file containing the structure. Read to get Volume. (optional)')

    output_file_group = parser.add_argument_group('Output file format')
    output_file_group.add_argument('-o', '--output', type=str, default='output',
            help='Prefix of the output files')
    output_file_group.add_argument('-O', '--bin-output', action='store_true',
            help='Save also binary files. (optional)')
    output_file_group.add_argument('--no-text-output', action='store_true',
            help='Do not save text files. (optional)')
    output_file_group.add_argument('--no-plot', action='store_true',
            help='Do not save plot files. (optional)')
    output_file_group.add_argument('--bin-output-old', action='store_true',
            help='Use old format for binary files (compatibility). (optional) - *TO BE DEPRECATED*')

    input_params_group = parser.add_argument_group('Physical parameters')
    input_params_group.add_argument('-t', '--timestep', type=float, required=True,
            help='Time step of the data (fs)')
    for parameter in st.current.all_parameters:
        input_params_group.add_argument(f'--{parameter}', type=float,
            help='Usually Angstrom or Kelvins, see description of units and currents implemented available with --list-currents')
    input_params_group.add_argument('-u', '--units', type=str, default='metal',
            choices=st.current.all_units,
            help='Units. (optional, default: metal)')
    input_params_group.add_argument('-C', '--current', type=str, default='heat',
            choices=list(st.current.all_currents.keys()),
            help='Type of currents that is provided to the code. Usually this just changes the conversion factor')
    input_params_group.add_argument('--param-from-input-file-column', type=str,
            action='append', dest='parameters_from_input_file', nargs=2,
            help='in order: header of the column and name of the parameter that will be setted to the average of that column of the input file')
    input_params_group.add_argument('--list-currents', action='store_true',
            help='show the list of currents implemented, the docstrings of the units, then exit')

    analysis_group = parser.add_argument_group('Analysis options')
    analysis_group.add_argument('-r', '--resample', action='store_true',
            help='Resample the time series (using TSKIP or FSTAR). (optional)')
    resamplearg = analysis_group.add_mutually_exclusive_group()
    resamplearg.add_argument('--TSKIP', type=int,
            help='Resampling time period (steps)')
    resamplearg.add_argument('--FSTAR', type=float,
            help='Resampling target Nyquist frequency (THz)')
    cutoff_group = analysis_group.add_mutually_exclusive_group()
    cutoff_group.add_argument('-c', '--corr-factor', type=float, default=1.0,
            help='Correction factor to the AIC. (optional, default: 1.0 = no correction)')
    cutoff_group.add_argument('-P', '--manual-Pstar', type=int,
            help='Manual P* value. (optional, default: use automatic P*)')

    plot_group = parser.add_argument_group('Plot options (optional)')
    plot_group.add_argument('-w', '--psd-filterw', type=float,
            help='Periodogram filter window width (THz)')
    plot_group.add_argument('--plot-conv-max-pstar', type=int,
            help='Max number of P* in the kappa(P*) plot (x)')
    plot_group.add_argument('--plot-conv-max-kappa', type=float,
            help='Max kappa in the kappa(P*) plot (y)')
    plot_group.add_argument('--plot-conv-pstar-tick-interval', type=int,
            help='Tick interval on the x-axis for the kappa(P*) plot')
    plot_group.add_argument('--plot-conv-kappa-tick-interval', type=float,
            help='Tick interval on the y-axis for the kappa(P*) plot')
    plot_group.add_argument('--plot-psd-max-THz', type=float,
            help='Max frequency (THz) in the psd plot (x)')
    plot_group.add_argument('--plot-psd-max-kappa', type=float,
            help='Max kappa (W/m/K) in the psd plot (y)')
    plot_group.add_argument('--plot-psd-THz-tick-interval', type=float,
            help='Tick interval on the x-axis for the psd plot')
    plot_group.add_argument('--plot-psd-kappa-tick-interval', type=float,
            help='Tick interval on the y-axis for the psd plot')

    beta_group = parser.add_argument_group('Testing options')
    beta_group.add_argument('--test-suite-run', action='store_true')
    beta_group.add_argument('--savetxt-format', type=str, default='%.18e',
            help='Format string used by `numpy.savetxt` in the output files')
    # yapf: enable
    args = parser.parse_args()

    run_analysis(args)
    return 0


def run_analysis(args):

    inputfile = args.inputfile
    input_format = args.input_format
    j1_key = args.mainfluxkey
    NSTEPS = args.nsteps
    START_STEP = args.start_step
    jindex = args.cindex
    sindex = args.sindex
    NSPLIT = args.split
    run_keyword = args.run_keyword
    structurefile = args.structure

    output = args.output
    binout = args.bin_output
    binout_old = args.bin_output_old
    do_plot = not args.no_plot
    no_text_out = args.no_text_output

    DT_FS = args.timestep
    parameters = {}
    for parameter in st.current.all_parameters:
        p = getattr(args, parameter)
        if p is not None:
            if p <= 0.:
                raise ValueError(f'{parameter} must be positive')
            parameters[parameter] = p
    parameters_from_input_file = args.parameters_from_input_file if args.parameters_from_input_file else []
    parameters_from_input_file_key = [x[0] for x in parameters_from_input_file]
    parameters_from_input_file_name = [x[1] for x in parameters_from_input_file]
    units = args.units
    current_type = args.current

    resample = args.resample
    TSKIP = args.TSKIP
    FSTAR = args.FSTAR
    corr_factor = args.corr_factor
    if args.manual_Pstar is not None:
        manual_cutoffK = args.manual_Pstar - 1
    else:
        manual_cutoffK = None
    j2_keys = args.add_currents

    psd_filter_w = args.psd_filterw

    print_elapsed = not args.test_suite_run
    print_cmd = not args.test_suite_run
    fmt = args.savetxt_format

    if DT_FS <= 0.:
        raise ValueError('Time step must be positive')
    if NSTEPS < 0:
        raise ValueError('nsteps must be positive')
    if resample:
        if TSKIP is not None:
            if TSKIP <= 1:
                raise ValueError('Resampling: TSKIP should be > 1')
        elif FSTAR is not None:
            if FSTAR <= 0.:
                raise ValueError('Resampling: FSTAR should be positive')
        else:
            raise ValueError('Resampling: you should specify either TSKIP or FSTAR')
    elif TSKIP is not None:
        raise ValueError('Use flag -r to resample. TSKIP will be ignored')
    elif FSTAR is not None:
        raise ValueError('Use flag -r to resample. FSTAR will be ignored')
    if corr_factor <= 0.:
        raise ValueError('The correction factor must be positive')
    if NSPLIT < 1:
        raise ValueError('The number of splits must be a positive number')

    log.open_file(output + '.log')
    if print_cmd:
        log.write_log('Command:\n ' + ' '.join(argv) + '\n\n')

    # Write some parameters
    if print_cmd:
        log.write_log(' Input file ({}):      {}'.format(input_format, inputfile))
    log.write_log(' Units:      {}'.format(units))
    log.write_log(' Time step:      {} fs'.format(DT_FS))

    # Read data
    selected_keys = [j1_key]
    selected_keys.extend(j2_keys)
    jdata = None
    if input_format == 'table':
        # Table format: data is organized in columns, the selected_keys determines which to read
        #input parameters that are read from file
        for col, pname in parameters_from_input_file:
            selected_keys.append(col)
        jfile = st.i_o.TableFile(inputfile, group_vectors=True, print_elapsed=print_elapsed)
        jfile.read_datalines(start_step=START_STEP, NSTEPS=NSTEPS, select_ckeys=selected_keys)
        jdata = jfile.data
        START_STEP = 0   # reset to zero, as later we will need to read all of jdata

    elif input_format == 'dict':
        # Dict format: data is stored in a binary Numpy file containing a dictionary
        jdata = np.load(inputfile, allow_pickle=True).tolist()

    elif input_format == 'lammps':
        # LAMMPS format: a LAMMPS log file is scanned until the run_keywork is found
        jfile = st.i_o.LAMMPSLogFile(inputfile, run_keyword=run_keyword)
        if args.TEMPERATURE is None:
            selected_keys.append('Temp')
        jfile.read_datalines(NSTEPS, select_ckeys=selected_keys)
        jdata = jfile.data

    else:
        raise NotImplementedError('Input format not implemented.')

    # split data
    if NSPLIT > 1:
        log.write_log('Splitting input data time series into {:d} segments...'.format(NSPLIT))
        data_size = jdata[selected_keys[0]].shape[0]
        if len(jdata[selected_keys[0]].shape) > 1:
            n_proc = jdata[selected_keys[0]].shape[1]
        else:
            n_proc = 1
        rm = data_size % NSPLIT
        steps_start = data_size - rm
        steps_end = data_size / NSPLIT
        if (steps_end % 2) == 1:
            steps_end = steps_end - 1
        for key, value in jdata.items():
            if not key in parameters_from_input_file_key:
                newdata = value[:steps_start].reshape((NSPLIT, data_size / NSPLIT, n_proc)).transpose(
                    (1, 0, 2)).reshape((data_size / NSPLIT, NSPLIT * n_proc))
                jdata[key] = newdata[:steps_end]
        log.write_log('New shape of input data: {}'.format(jdata[selected_keys[0]].shape))

    if NSTEPS == 0:
        NSTEPS = jdata[list(jdata.keys())[0]].shape[0]

    # compute average parameters from input file, if requested
    def average(data, name, units=''):
        ave = np.mean(data)
        std = np.std(data)
        log.write_log(f'Mean {name} (computed): {ave} +/- {std}')
        return ave

    for key, value in parameters.items():
        log.write_log(f'{key} (input): {value}')
    for key, name in parameters_from_input_file:
        parameters[name] = average(jdata[key], name)
        selected_keys.remove(key)

    if structurefile is not None:
        # read volume from LAMMPS data file
        _, volume = st.i_o.read_lammps_datafile.get_box(structurefile)
        log.write_log(' Volume (structure file):    {} A^3'.format(volume))
        #note: here I hardcoded the volume key
        #      nothing guarantees that in the parameter list
        #      of the function that calculates KAPPA_SCALE
        #      you are going to find the VOLUME parameter with this meaning
        parameters['VOLUME'] = volume

    # Time step
    log.write_log(' Time step (input):  {} fs'.format(DT_FS))

    # Define currents
    if jindex is None:
        # read all components
        currents = np.array([jdata[key][START_STEP:(START_STEP + NSTEPS), :] for key in selected_keys])
    else:
        # read only the components jindex
        # NOTE: for multi-current cases, it will select jindex of each current
        if sindex is None:
            currents = np.array([jdata[key][START_STEP:(START_STEP + NSTEPS), jindex] for key in selected_keys])
        else:
            # subtract the components sindex from those jindex
            currents = np.array([
                jdata[key][START_STEP:(START_STEP + NSTEPS), jindex] -
                jdata[key][START_STEP:(START_STEP + NSTEPS), sindex] for key in selected_keys
            ])
    log.write_log('  currents shape is {}'.format(currents.shape))
    log.write_log('snippet:')
    log.write_log(currents)

    # create HeatCurrent object
    j = st.current.all_currents[current_type][0](currents, DT_FS=DT_FS, UNITS=units, **parameters,
                                                 PSD_FILTER_W=psd_filter_w)

    log.write_log(' Number of currents = {}'.format(j.N_CURRENTS))
    log.write_log(' Number of equivalent components = {}'.format(j.N_EQUIV_COMPONENTS))
    log.write_log(' KAPPA_SCALE = {}'.format(j.KAPPA_SCALE))
    log.write_log(' Nyquist_f   = {}  THz'.format(j.Nyquist_f_THz))

    # resample
    if resample:
        if TSKIP is not None:
            jf = j.resample(TSKIP=TSKIP, PSD_FILTER_W=psd_filter_w)
            #FSTAR = j.Nyquist_f_THz / TSKIP   # from st.heatcurrent.resample_current
            FSTAR = jf.Nyquist_f_THz
        else:
            jf = j.resample(fstar_THz=FSTAR, PSD_FILTER_W=psd_filter_w)
        #log.write_log(jf.resample_log)
    else:
        jf = j

    # cepstral analysis
    jf.cepstral_analysis(aic_type='aic', aic_Kmin_corrfactor=corr_factor, manual_cutoffK=manual_cutoffK)
    #log.write_log(jf.cepstral_log)

    ############################################################################
    ## OUTPUT SECTION
    ## TODO: cleanup data files
    ############################################################################

    # DATA OUTPUT
    if binout:
        binoutobj = TCOutput()

        binoutobj.j_DT_FS = j.DT_FS
        binoutobj.j_freqs_THz = j.freqs_THz
        binoutobj.j_fpsd = j.fpsd
        binoutobj.j_flogpsd = j.flogpsd
        binoutobj.j_psd = j.psd
        binoutobj.j_logpsd = j.logpsd
        binoutobj.j_Nyquist_f_THz = j.Nyquist_f_THz
        binoutobj.j_PSD_FILTER_W_THz = psd_filter_w
        if j.MANY_CURRENTS:
            binoutobj.j_cospectrum = j.cospectrum
            binoutobj.j_fcospectrum = j.fcospectrum

        if resample:
            binoutobj.jf_DT_FS = jf.DT_FS
            binoutobj.jf_freqs_THz = jf.freqs_THz
            binoutobj.jf_fpsd = jf.fpsd
            binoutobj.jf_flogpsd = jf.flogpsd
            binoutobj.jf_psd = jf.psd
            binoutobj.jf_logpsd = jf.logpsd
            binoutobj.jf_Nyquist_f_THz = jf.Nyquist_f_THz
            binoutobj.jf_resample_log = jf.resample_log

        binoutobj.kappa = jf.kappa
        binoutobj.kappa_std = jf.kappa_std
        binoutobj.cepstral_log = jf.cepstral_log
        binoutobj.units = jf.UNITS
        binoutobj.KAPPA_SCALE = jf.KAPPA_SCALE
        binoutobj.TEMPERATURE = jf.TEMPERATURE
        binoutobj.VOLUME = jf.VOLUME

        binoutobj.jf_cepf_logpsdK = jf.cepf.logpsdK
        binoutobj.jf_cepf_logpsdK_THEORY_std = jf.cepf.logpsdK_THEORY_std
        binoutobj.jf_cepf_logtau = jf.cepf.logtau
        binoutobj.jf_cepf_logtau_THEORY_std = jf.cepf.logtau_THEORY_std
        binoutobj.jf_cepf_kappa = jf.cepf.tau * jf.KAPPA_SCALE * 0.5
        binoutobj.jf_cepf_kappa_THEORY_std = jf.cepf.tau_THEORY_std * jf.KAPPA_SCALE * 0.5
        binoutobj.jf_cepf_aic_Kmin = jf.cepf.aic_Kmin
        binoutobj.jf_cepf_aic_Kmin_corrfactor = jf.cepf.aic_Kmin_corrfactor
        binoutobj.jf_cepf_cutoffK = jf.cepf.cutoffK
        binoutobj.jf_cepf_psd = jf.cepf.psd
        binoutobj.jf_cepf_logpsd = jf.cepf.logpsd

        if binout_old:
            binoutobj.write_old_binary(output)
        else:
            np.save(output, binoutobj)

    if not no_text_out:
        outfile_name = output + '.psd.dat'
        outarray = np.c_[j.freqs_THz, j.psd, j.fpsd, j.logpsd, j.flogpsd]
        outfile_header = 'freqs_THz  psd  fpsd  logpsd  flogpsd\n'
        np.savetxt(outfile_name, outarray, header=outfile_header, fmt=fmt)
        if j.MANY_CURRENTS:
            outfile_name = output + '.cospectrum.dat'
            outarray = np.c_[j.freqs_THz,
                             j.cospectrum.reshape(
                                 (j.cospectrum.shape[0] * j.cospectrum.shape[1], j.cospectrum.shape[2])).transpose()]
            np.savetxt(outfile_name, np.column_stack([outarray.real, outarray.imag]), fmt=fmt)

            outfile_name = output + '.cospectrum.filt.dat'
            outarray = np.c_[j.freqs_THz,
                             j.fcospectrum.reshape(
                                 (j.fcospectrum.shape[0] * j.fcospectrum.shape[1], j.fcospectrum.shape[2])).transpose()]
            np.savetxt(outfile_name, np.column_stack([outarray.real, outarray.imag]), fmt=fmt)

        if resample:
            outfile_name = output + '.resampled_psd.dat'
            outarray = np.c_[jf.freqs_THz, jf.psd, jf.fpsd, jf.logpsd, jf.flogpsd]
            outfile_header = 'freqs_THz  psd  fpsd  logpsd  flogpsd\n'
            np.savetxt(outfile_name, outarray, header=outfile_header, fmt=fmt)

        outfile_name = output + '.cepstral.dat'
        outarray = np.c_[jf.cepf.logpsdK, jf.cepf.logpsdK_THEORY_std, jf.cepf.logtau, jf.cepf.logtau_THEORY_std,
                         jf.cepf.tau * jf.KAPPA_SCALE * 0.5, jf.cepf.tau_THEORY_std * jf.KAPPA_SCALE * 0.5]
        outfile_header = 'ck  ck_std  L0(P*)  L0_std(P*)  kappa(P*)  kappa_std(P*)\n'
        np.savetxt(outfile_name, outarray, header=outfile_header, fmt=fmt)

        outfile_name = output + '.cepstrumfiltered_psd.dat'
        outarray = np.c_[jf.freqs_THz, jf.cepf.psd, jf.cepf.logpsd]
        outfile_header = 'freqs_THz  cepf_psd cepf_logpsd\n'
        np.savetxt(outfile_name, outarray, header=outfile_header, fmt=fmt)

    ####################################
    # PLOTS
    ####################################

    if do_plot:
        pdf = PdfPages(output + '.plots.pdf')

        addPlotToPdf(j.plot_periodogram, pdf)
        if resample:
            ax = j.plot_resample(jf)
            ax[0].set_xlim([0, 2.5 * FSTAR])
            pdf.savefig()
            plt.close()
        addPlotToPdf(j.plot_psd, pdf, jf, f_THz_max=args.plot_psd_max_THz, k_SI_max=args.plot_psd_max_kappa,
                     k_tick=args.plot_psd_kappa_tick_interval, f_tick=args.plot_psd_THz_tick_interval)
        addPlotToPdf(jf.plot_psd, pdf, f_THz_max=args.plot_psd_max_THz, k_SI_max=args.plot_psd_max_kappa,
                     k_tick=args.plot_psd_kappa_tick_interval, f_tick=args.plot_psd_THz_tick_interval)
        addPlotToPdf(jf.plot_psd, pdf, jf, jf, f_THz_max=args.plot_psd_max_THz, k_SI_max=args.plot_psd_max_kappa,
                     k_tick=args.plot_psd_kappa_tick_interval, f_tick=args.plot_psd_THz_tick_interval)
        try:
            for idx1 in range(j.N_CURRENTS):
                for idx2 in range(idx1, j.N_CURRENTS):
                    addPlotToPdf(j.plot_cospectrum_component, pdf, idx1, idx2)
        except:
            pass

        # plot cepstral coefficients
        ax = jf.plot_ck()
        ax.set_xlim([0, 5 * jf.cepf.cutoffK])
        ax.set_ylim([-1, 15])
        ax.grid()
        pdf.savefig()
        plt.close()

        # plot L0(Pstar)
        ax = jf.plot_L0_Pstar()
        ax.set_xlim([0, 10 * jf.cepf.cutoffK])
        pdf.savefig()
        plt.close()

        # # plot kappa(Pstar)
        # ax = jf.plot_kappa_Pstar()
        # ax.set_xlim([0, 10*jf.cepf.cutoffK])

        addPlotToPdf(jf.plot_kappa_Pstar, pdf, pstar_max=args.plot_conv_max_pstar,
                     pstar_tick=args.plot_conv_pstar_tick_interval, kappa_SI_max=args.plot_conv_max_kappa,
                     kappa_tick=args.plot_conv_kappa_tick_interval)

        # plot cepstral log-PSD
        ax = jf.plot_periodogram()
        jf.plot_cepstral_spectrum(axes=ax, label='cepstrum-filtered')
        ax[0].axvline(x=jf.Nyquist_f_THz, ls='--', c='r')
        ax[1].axvline(x=jf.Nyquist_f_THz, ls='--', c='r')
        #ax[0].set_xlim([0., 2.5*FSTAR_THZ])
        #ax[1].set_ylim([12,18])
        #ax[0].legend(['original', 'resampled', 'cepstrum-filtered'])
        #ax[1].legend(['original', 'resampled', 'cepstrum-filtered']);

        pdf.close()

    log.close_file()

    return 0


#################################


class TCOutput(object):

    # yapf: disable
    def __init__(self):
        # TO BE COMPLETED WIHT ALL PARAMETERS
        self.j_DT_FS            = None
        self.j_freqs_THz        = None
        self.j_fpsd             = None
        self.j_flogpsd          = None
        self.j_psd              = None
        self.j_logpsd           = None
        self.j_Nyquist_f_THz    = None
        self.j_PSD_FILTER_W_THz = None
        self.j_cospectrum       = None

        self.jf_DT_FS         = None
        self.jf_freqs_THz     = None
        self.jf_fpsd          = None
        self.jf_flogpsd       = None
        self.jf_psd           = None
        self.jf_logpsd        = None
        self.jf_Nyquist_f_THz = None
        self.jf_resample_log  = None

        self.jf_cepf_logpsdK            = None
        self.jf_cepf_logpsdK_THEORY_std = None
        self.jf_cepf_logtau             = None
        self.jf_cepf_logtau_THEORY_std  = None
        self.jf_cepf_kappa              = None
        self.jf_cepf_kappa_THEORY_std   = None
        self.jf_cepf_psd                = None
        self.jf_cepf_logpsd             = None
        self.jf_cepf_aic_Kmin           = None
        self.jf_cepf_aic_Kmin_corrfactor = None
        self.jf_cepf_cutoffK            = None

        self.kappa          = None
        self.kappa_std      = None
        self.cepstral_log   = None
        self.UNITS          = None
        self.KAPPA_SCALE    = None
        self.TEMPERATURE    = None
        self.VOLUME         = None
        self.TSKIP          = None

    # yapf: enable
    def write_old_binary(self, output):
        """Write old binary format."""
        opts = {'allow_pickle': False}
        optsa = {'axis': 1}
        outarray = np.c_[self.j_freqs_THz, self.j_fpsd, self.j_flogpsd, self.j_psd, self.j_logpsd]
        np.save(output + '.psd.npy', outarray, **opts)

        if self.j_cospectrum is not None:
            outarray = np.c_[self.j_freqs_THz, self.j_cospectrum.reshape(-1, self.j_cospectrum.shape[-1]).transpose()]
            np.save(output + '.cospectrum.npy', outarray, **opts)

        if self.j_fcospectrum is not None:
            outarray = np.c_[self.j_freqs_THz, self.j_fcospectrum.reshape(-1, self.j_fcospectrum.shape[-1]).transpose()]
            np.save(output + '.cospectrum.filt.npy', outarray, **opts)

        outarray = np.c_[self.jf_freqs_THz, self.jf_psd, self.jf_fpsd, self.jf_logpsd, self.jf_flogpsd]
        np.save(output + '.resampled_psd.npy', outarray, **opts)

        outarray = np.c_[self.jf_cepf_logpsdK, self.jf_cepf_logpsdK_THEORY_std, self.jf_cepf_logtau,
                         self.jf_cepf_logtau_THEORY_std, self.jf_cepf_kappa, self.jf_cepf_kappa_THEORY_std]
        np.save(output + '.cepstral', outarray, **opts)

        outarray = np.c_[self.jf_freqs_THz, self.jf_cepf_psd, self.jf_cepf_logpsd]
        np.save(output + '.cepstrumfiltered_psd', outarray, **opts)


if __name__ == '__main__':
    main()
