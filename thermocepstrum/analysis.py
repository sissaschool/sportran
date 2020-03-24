#!/usr/bin/env python

import os
from sys import path, argv
import argparse
import numpy as np
#import scipy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator

try:
    import thermocepstrum as tc
except ImportError:
    abs_path = os.path.abspath(__file__)
    tc_path = abs_path[:abs_path.rfind('/')]
    path.append(tc_path[:tc_path.rfind('/')])
    try:
        import thermocepstrum as tc
    except ImportError:
        raise ImportError('Cannot locate thermocepstrum.')

# import log-print method
from thermocepstrum.utils.utils import PrintMethod
log = PrintMethod()
log.set_method('bash')

#try to import plotstyle (not fundamental)
try:
    import pkg_resources
    pltstyle_filename = pkg_resources.resource_filename('thermocepstrum.utils', 'plot_style.mplstyle')
except:
    # fallback (maybe it is not installed...)
    try:
        abs_path = os.path.abspath(__file__)
        tc_path = abs_path[:abs_path.rfind('/')]
        path.append(tc_path[:tc_path.rfind('/')])
    except:
        abs_path = '.'
    pltstyle_filename = tc_path + '/utils/plot_style.mplstyle'
try:
    plt.style.use(pltstyle_filename)
except:
    pass

c = plt.rcParams['axes.prop_cycle'].by_key()['color']


def main():
    """
--------------------------------------------------------------------------------
  *** THERMOCEPSTRUM ***  command line interface  (beta)
--------------------------------------------------------------------------------
This script performs the cepstral analysis. It outputs some results in the stdout and log file, and plots in pdf format.

INPUT FORMAT:
 - table  : a column-formatted text file, with a header in the same format of LAMMPS. The name of the LAMMPS compute can start with c_ and end with [#some_number], the code will recognize vectors, and will read automatically all the components.
 - dict   : a Numpy binary file containing a dictionary (e.g. obtained from the script i_o/read_lammps_log.py)
 - LAMMPS : a LAMMPS log file. In this case a --run-keyword  must be provided, that identifies the run to be read (see documentation of i_o/read_lammps_log.py)
The average temperature is computed if a column with the header 'Temp' is found; otherwise you have to specify it.
You must provide the name of the heat flux compute. You can also provide additional currents if your system is a multi-component fluid.
(Notice that the output is the same with any number of components. If you have a lots of components, note that you may want to use more than 3 independent processes -- see theory.)
Units can be metal or real (see LAMMPS documentation at http://lammps.sandia.gov/doc/units.html )

OUTPUT files:
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
  read and analyze "example/Silica.dat" file. The heat-flux columns are called c_flux[1], c_flux[2], c_flux[3]

    ./analysis "example/Silica.dat" -V 3130.431110818 -T 1065.705630 -t 1.0 -k flux1 -u metal -r --FSTAR 28.0 -w 0.5 -o silica_test

-------------------------
"""
    _epilog = """---
Enjoy it!
---
This software was written by Loris Ercole and extended by Riccardo Bertossa to handle the multicomponent stuff, at SISSA, Via Bonomea, 265 - 34136 Trieste ITALY.

Please cite these references:
 - Ercole, Marcolongo, Baroni, Sci. Rep. 7, 15835 (2017), https://doi.org/10.1038/s41598-017-15843-2
 - Bertossa, Grasselli, Ercole, Baroni Phys. Rev. Lett. 122, 255901 (2019), https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.255901
 - Baroni, Bertossa, Ercole, Grasselli, Marcolongo, https://arxiv.org/abs/1802.08006

https://github.com/lorisercole/thermocepstrum
Contact: lercole@sissa.it
"""

    # yapf: disable
    parser = argparse.ArgumentParser(description=main.__doc__, epilog=_epilog, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument( 'inputfile', type=str, help='input file to read (default format: Table)' )
    parser.add_argument( '-t', '--timestep', type=float, required=True, help='Time step of the log.write_loged data (fs)' )
    parser.add_argument( '-k', '--heatfluxkey', type=str, required=True, help='Name of the column keyword that identifies the heat flux' )
    parser.add_argument( '-N', '--nsteps', type=int, default=0, help='Number of steps to read (default: 0=all)' )
    parser.add_argument( '-S', '--start-step', type=int, default=0, help='The first step to read (default: 0=first)' )
    parser.add_argument( '--input-format', default='table', type=str, choices=['table','dict','lammps'], help='Format of the input file' )
    parser.add_argument( '--cindex', nargs='*', type=int, help='Column indexes of the heatflux to read (0,1,2,...)' )
    parser.add_argument( '--sindex', nargs='*', type=int, help='Column indexes of the heatflux to substract from the flux read with --cindex (3,4,5,...)' )
    parser.add_argument( '--run-keyword', type=str, help='Keyword that identifies the run to be read (only for "lammps" format)' )
    parser.add_argument( '--split', type=int, default=1, help='Build a time series with n*m independent processes (n is the number of processes of the original timeseries, m is the number provided with --split). The length of the new time series will be [original length]/m.')

    parser.add_argument( '-o', '--output', type=str, default='output', help='prefix of the output files' )
    parser.add_argument( '-O', '--bin-output', action='store_true', help='save also binary files' )
    parser.add_argument( '--no-text-output', action='store_true', help='do not save text files' )
    parser.add_argument( '--bin-output-old', action='store_true', help='use old format for binary files (compatibility)' )

    parser.add_argument( '-V', '--volume', type=float, help='Volume of the cell (Angstrom). If not set it will be read from structure file or inputfile.' )
    parser.add_argument( '--structure', type=str, help='LAMMPS data file containing the structure. Read to get Volume.' )
    parser.add_argument( '-T', '--temperature', type=float, help='average Temperature (K). If not set it will be read from file' )
    parser.add_argument( '-u', '--units', type=str, default='metal', choices=['metal', 'real', 'qepw', 'gpumd', 'dlpoly'], help='LAMMPS units (default: metal)' )

    parser.add_argument( '-r', '--resample', action='store_true', help='resample the time series (you should define --TSKIP or --FSTAR' )
    resamplearg = parser.add_mutually_exclusive_group()
    resamplearg.add_argument( '--TSKIP', type=int, help='resampling time period (steps)' )
    resamplearg.add_argument( '--FSTAR', type=float, help='resampling target Nyquist frequency (THz)' )
    parser.add_argument( '-c', '--corr-factor', type=float, default=1.0, help='correction factor to the AIC' )
    parser.add_argument( '-j', '--add-currents', type=str, default=[], action='append', help='additional current for multi-component fluids' )

    parser.add_argument( '-w', '--psd-filterw', type=float, help='plot - periodogram - filter window width (THz)' )
    parser.add_argument('--plot-conv-max-pstar',type=int,help='max number of P* in the kappa(P*) plot (x)')
    parser.add_argument('--plot-conv-max-kappa',type=float,help='max kappa in the kappa(P*) plot (y)')
    parser.add_argument('--plot-conv-pstar-tick-interval',type=int,help='tick interval on the x-axis for the kappa(P*) plot')
    parser.add_argument('--plot-conv-kappa-tick-interval',type=float,help='tick interval on the y-axis for the kappa(P*) plot')
    parser.add_argument('--plot-psd-max-THz',type=float, help='max frequency in THz for the psd plot (x)')
    parser.add_argument('--plot-psd-max-kappa',type=float, help='max kappa in W/mK for the psd plot (y)')
    parser.add_argument('--plot-psd-THz-tick-interval',type=float, help='tick interval on the x-axis for the psd plot')
    parser.add_argument('--plot-psd-kappa-tick-interval',type=float, help='tick interval on the y-axis for the psd plot')
    args = parser.parse_args()

    # yapf: enable
    inputfile = args.inputfile
    DT_FS = args.timestep
    j1_key = args.heatfluxkey
    NSTEPS = args.nsteps
    START_STEP = args.start_step
    input_format = args.input_format
    jindex = args.cindex
    sindex = args.sindex
    run_keyword = args.run_keyword
    NSPLIT = args.split

    output = args.output
    binout = args.bin_output
    binout_old = args.bin_output_old
    no_text_out = args.no_text_output

    units = args.units
    temperature = args.temperature
    volume = args.volume
    structurefile = args.structure
    resample = args.resample
    TSKIP = args.TSKIP
    FSTAR = args.FSTAR
    corr_factor = args.corr_factor
    j2_keys = args.add_currents
    psd_filter_w = args.psd_filterw

    if volume is not None:
        if (volume <= 0.):
            raise ValueError('volume must be positive')
    if (DT_FS <= 0.):
        raise ValueError('timestep must be positive')
    if (NSTEPS < 0):
        raise ValueError('nsteps must be positive')
    if temperature is not None:
        if (temperature <= 0.):
            raise ValueError('temperature must be positive')
    if resample:
        if (TSKIP is not None):
            if (TSKIP <= 1):
                raise ValueError('resampling: TSKIP should be > 1')
        elif (FSTAR is not None):
            if (FSTAR <= 0.):
                raise ValueError('resampling: FSTAR should be positive')
        else:
            raise ValueError('resampling: you should specify either TSKIP or FSTAR')
    elif (TSKIP is not None):
        raise ValueError('Use flag -r to resample. TSKIP will be ignored')
    elif (FSTAR is not None):
        raise ValueError('Use flag -r to resample. FSTAR will be ignored')
    if (corr_factor <= 0.):
        raise ValueError('the correction factor must be positive')
    if (NSPLIT < 1):
        raise ValueError('The number of splits must be a positive number')

    ncurrents = len(j2_keys) + 1

    logfile = open(output + '.log', 'w')
    logfile.write('Command:\n ' + ' '.join(argv) + '\n\n')

    selected_keys = [j1_key]
    selected_keys.extend(j2_keys)

    # Write some parameters
    log.write_log(' Input file ({}):      {}'.format(input_format, inputfile))
    logfile.write(' Input file ({}):      {}\n'.format(input_format, inputfile))
    log.write_log(' Units:      {}'.format(units))
    logfile.write(' Units:      {}\n'.format(units))
    log.write_log(' Time step:      {} fs'.format(DT_FS))
    logfile.write(' Time step:      {} fs\n'.format(DT_FS))

    ## Read data
    jdata = None
    if (input_format == 'table'):
        if temperature is None:
            selected_keys.append('Temp')
#      if 'Press' in jfile.ckey:
#         selected_keys.append('Press')
        jfile = tc.i_o.TableFile(inputfile, group_vectors=True)
        jfile.read_datalines(start_step=START_STEP, NSTEPS=NSTEPS, select_ckeys=selected_keys)
        jdata = jfile.data
        START_STEP = 0   # reset to zero, as later we will need to read all of jdata
    elif (input_format == 'dict'):
        jdata = np.load(inputfile, allow_pickle=True).tolist()
    elif (input_format == 'lammps'):
        jfile = tc.i_o.LAMMPSLogFile(inputfile, run_keyword=run_keyword)
        if temperature is None:
            selected_keys.append('Temp')
#      if 'Press' in jfile.ckey:
#         selected_keys.append('Press')
        jfile.read_datalines(NSTEPS, select_ckeys=selected_keys)
        jdata = jfile.data
    else:
        raise NotImplemented('input format not implemented.')

    if (NSPLIT > 1):
        log.write_log('Splitting input data time series into {:d} segments...'.format(NSPLIT))
        logfile.write('Splitting input data time series into {:d} segments...\n'.format(NSPLIT))
        data_size = jdata[selected_keys[0]].shape[0]
        n_proc = 1
        try:   ## HORRIBLE
            n_proc = jdata[selected_keys[0]].shape[1]
        except:
            pass
        rm = data_size % NSPLIT
        steps_start = data_size - rm
        steps_end = data_size / NSPLIT
        if (steps_end % 2 == 1):
            steps_end = steps_end - 1
        for key, value in jdata.items():
            if key != 'Temp':
                newdata = value[:steps_start].reshape((NSPLIT, data_size / NSPLIT, n_proc)).transpose(
                    (1, 0, 2)).reshape((data_size / NSPLIT, NSPLIT * n_proc))
                jdata[key] = newdata[:steps_end]
        log.write_log('New shape of input data: {}'.format(jdata[selected_keys[0]].shape))
        logfile.write('New shape of input data: {}\n'.format(jdata[selected_keys[0]].shape))

    if (NSTEPS == 0):
        NSTEPS = jdata[list(jdata.keys())[0]].shape[0]

    ## Define Temperature
    if temperature is None:
        if 'Temp' in jdata:
            temperature = np.mean(jdata['Temp'])
            temperature_std = np.std(jdata['Temp'])   # this is wrong (needs block average)
            if 'Temp' in selected_keys:
                selected_keys.remove('Temp')
            log.write_log(' Mean Temperature (computed):  {} K  +/-  {}'.format(temperature, temperature_std))
            logfile.write(' Mean Temperature (computed):  {} K  +/-  {}\n'.format(temperature, temperature_std))
        elif 'Temp_ave' in jdata:
            temperature = jdata['Temp_ave']
            if 'Temp_std' in jdata:
                temperature_std = jdata['Temp_std']
                log.write_log(' Mean Temperature (file):      {} K  +/-  {}'.format(temperature, temperature_std))
                logfile.write(' Mean Temperature (file):      {} K  +/-  {}\n'.format(temperature, temperature_std))
            else:
                log.write_log(' Mean Temperature (file):      {} K'.format(temperature))
                logfile.write(' Mean Temperature (file):      {} K\n'.format(temperature))
        else:
            raise RuntimeError('No Temp key found. Please provide Temperature (-T).')
    else:
        log.write_log(' Mean Temperature (input):  {} K'.format(temperature))
        logfile.write(' Mean Temperature (input):  {} K\n'.format(temperature))

    ## Define Volume
    if volume is None:
        if structurefile is not None:
            _, volume = tc.i_o.read_lammps_datafile.get_box(structurefile)
            log.write_log(' Volume (structure file):    {} A^3'.format(volume))
            logfile.write(' Volume (structure file):    {} A^3'.format(volume))
        elif 'Volume' in jdata:
            volume = jdata['Volume']
            log.write_log(' Volume (file):    {} A^3'.format(volume))
            logfile.write(' Volume (file):    {} A^3\n'.format(volume))
        else:
            raise RuntimeError('No Volume key found. Please provide Volume (-V) of structure file (--structure).')
    else:
        log.write_log(' Volume (input):  {} A^3'.format(volume))
        logfile.write(' Volume (input):  {} A^3\n'.format(volume))

    log.write_log(' Time step (input):  {} fs'.format(DT_FS))
    logfile.write(' Time step (input):  {} fs\n'.format(DT_FS))

    ### Compute Pressure (optional)
    #if 'Press' in jdata:
    #   pressure = np.mean(jdata['Press'])
    #   log.write_log ' Mean Pressure (computed):  {} K'.format(pressure)
    #   logfile.write(' Mean Pressure (computed):  {} K'.format(pressure))

    ## Define currents
    log.write_log(selected_keys, jindex)
    if jindex is None:
        currents = np.array([jdata[key][START_STEP:(START_STEP + NSTEPS), :] for key in selected_keys])
    else:
        if sindex is None:
            currents = np.array([jdata[key][START_STEP:(START_STEP + NSTEPS), jindex] for key in selected_keys])
        else:
            currents = np.array([
                jdata[key][START_STEP:(START_STEP + NSTEPS), jindex] -
                jdata[key][START_STEP:(START_STEP + NSTEPS), sindex] for key in selected_keys
            ])
    log.write_log('  currents shape is {}'.format(currents.shape))
    logfile.write('  currents shape is {}\n'.format(currents.shape))
    log.write_log('snippet:')
    log.write_log(currents)

    # create HeatCurrent object
    j = tc.heatcurrent.HeatCurrent(currents, units, DT_FS, temperature, volume, psd_filter_w)

    log.write_log(' Number of currents = {}'.format(ncurrents))
    logfile.write(' Number of currrents = {}\n'.format(ncurrents))
    log.write_log(' Number of components = {}'.format(j.N_COMPONENTS))
    logfile.write(' Number of components = {}\n'.format(j.N_COMPONENTS))
    log.write_log(' kappa_scale = {}'.format(j.kappa_scale))
    logfile.write(' kappa_scale = {}\n'.format(j.kappa_scale))
    log.write_log(' Nyquist_f   = {}  THz'.format(j.Nyquist_f_THz))
    logfile.write(' Nyquist_f   = {}  THz\n'.format(j.Nyquist_f_THz))

    ################################
    ## OUTPUT SECTION   ## TODO: isolate output from computation
    if binout:
        binoutobj = TCOutput()

    with PdfPages(output + '.plots.pdf') as pdf:

        # plot periodogram
        j.plot_periodogram()   #PSD_FILTER_W=psd_filter_w)
        #pdf.attach_note('Nyquist frequency = {} THz'.format(j.Nyquist_f_THz), positionRect=[0, 0 , 5, 1])
        pdf.savefig()
        plt.close()

        if binout:
            binoutobj.j_DT_FS = j.DT_FS
            binoutobj.j_freqs_THz = j.freqs_THz
            binoutobj.j_fpsd = j.fpsd
            binoutobj.j_flogpsd = j.flogpsd
            binoutobj.j_psd = j.psd
            binoutobj.j_logpsd = j.logpsd
            binoutobj.j_Nyquist_f_THz = j.Nyquist_f_THz
            binoutobj.j_PSD_FILTER_W_THz = psd_filter_w
            if j.many_currents:
                binoutobj.j_cospectrum = j.cospectrum
                binoutobj.j_fcospectrum = j.fcospectrum
        #TODO: move all output in one place?
        if not no_text_out:
            outfile_name = output + '.psd.dat'
            outarray = np.c_[j.freqs_THz, j.psd, j.fpsd, j.logpsd, j.flogpsd]
            outfile_header = 'freqs_THz  psd  fpsd  logpsd  flogpsd\n'
            np.savetxt(outfile_name, outarray, header=outfile_header)
            if j.many_currents:
                outfile_name = output + '.cospectrum.dat'
                outarray = np.c_[j.freqs_THz,
                                 j.cospectrum.reshape((j.cospectrum.shape[0] * j.cospectrum.shape[1],
                                                       j.cospectrum.shape[2])).transpose()]
                np.savetxt(outfile_name, outarray)

                outfile_name = output + '.cospectrum.filt.dat'
                outarray = np.c_[j.freqs_THz,
                                 j.fcospectrum.reshape((j.fcospectrum.shape[0] * j.fcospectrum.shape[1],
                                                        j.fcospectrum.shape[2])).transpose()]
                np.savetxt(outfile_name, outarray)

        # resample and plot
        if resample:
            if TSKIP is not None:
                jf, ax = tc.heatcurrent.resample_current(j, TSKIP=TSKIP, plot=True, PSD_FILTER_W=psd_filter_w)
            else:
                jf, ax = tc.heatcurrent.resample_current(j, fstar_THz=FSTAR, plot=True, PSD_FILTER_W=psd_filter_w)
            ax[0].set_xlim([0, 2.5 * FSTAR])
            pdf.savefig()
            plt.close()
            logfile.write(jf.resample_log)

            # plot resampled periodogram
            ax = jf.plot_periodogram()   #PSD_FILTER_W=psd_filter_w)
            pdf.savefig()
            plt.close()

            if binout:
                binoutobj.jf_DT_FS = jf.DT_FS
                binoutobj.jf_freqs_THz = jf.freqs_THz
                binoutobj.jf_fpsd = jf.fpsd
                binoutobj.jf_flogpsd = jf.flogpsd
                binoutobj.jf_psd = jf.psd
                binoutobj.jf_logpsd = jf.logpsd
                binoutobj.jf_Nyquist_f_THz = jf.Nyquist_f_THz
                binoutobj.jf_resample_log = jf.resample_log
            if not no_text_out:
                outfile_name = output + '.resampled_psd.dat'
                outarray = np.c_[jf.freqs_THz, jf.psd, jf.fpsd, jf.logpsd, jf.flogpsd]
                outfile_header = 'freqs_THz  psd  fpsd  logpsd  flogpsd\n'
                np.savetxt(outfile_name, outarray, header=outfile_header)
        else:
            jf = j

        plt_psd(j, jf,\
                  f_THz_max=args.plot_psd_max_THz,\
                  k_SI_max=args.plot_psd_max_kappa,\
                  k_tick=args.plot_psd_kappa_tick_interval,\
                  f_tick=args.plot_psd_THz_tick_interval)
        pdf.savefig()
        plt.close()
        plt_psd(jf,\
                f_THz_max=args.plot_psd_max_THz,\
                k_SI_max=args.plot_psd_max_kappa,\
                k_tick=args.plot_psd_kappa_tick_interval,\
                f_tick=args.plot_psd_THz_tick_interval)
        pdf.savefig()
        plt.close()

        # cepstral analysis
        jf.cepstral_analysis(aic_type='aic', Kmin_corrfactor=corr_factor)
        logfile.write(jf.cepstral_log)
        if binout:
            binoutobj.kappa_Kmin = jf.kappa_Kmin
            binoutobj.kappa_Kmin_std = jf.kappa_Kmin_std
            binoutobj.cepstral_log = jf.cepstral_log
            binoutobj.units = jf.units
            binoutobj.kappa_scale = jf.kappa_scale
            binoutobj.TEMPERATURE = temperature
            binoutobj.VOLUME = volume

        plt_psd(jf, jf, jf,\
                f_THz_max=args.plot_psd_max_THz,\
                k_SI_max=args.plot_psd_max_kappa,\
                k_tick=args.plot_psd_kappa_tick_interval,\
                f_tick=args.plot_psd_THz_tick_interval)
        pdf.savefig()
        plt.close()

        try:
            for idx1 in range(ncurrents):
                for idx2 in range(idx1, ncurrents):
                    plt_other(j, idx1, idx2)
                    pdf.savefig()
                    plt.close()
        except:
            pass

        # plot cepstral coefficients
        ax = jf.plot_ck()
        ax.set_xlim([0, 5 * jf.dct.aic_Kmin])
        #  ax.set_ylim([-0.5, 0.5])
        ax.grid()
        pdf.savefig()
        plt.close()

        # plot L0(Pstar)
        ax = jf.plot_L0_Pstar()
        ax.set_xlim([0, 10 * jf.dct.aic_Kmin])
        pdf.savefig()
        plt.close()

        # plot kappa(Pstar)
        #      ax = jf.plot_kappa_Pstar()
        #      ax.set_xlim([0, 10*jf.dct.aic_Kmin])

        plt_cepstral_conv(jf,\
                          pstar_max=args.plot_conv_max_pstar,\
                          pstar_tick=args.plot_conv_pstar_tick_interval,\
                          k_SI_max=args.plot_conv_max_kappa,\
                          kappa_tick=args.plot_conv_kappa_tick_interval)
        pdf.savefig()
        plt.close()

        if binout:
            binoutobj.jf_dct_logpsdK = jf.dct.logpsdK
            binoutobj.jf_dct_logpsdK_THEORY_std = jf.dct.logpsdK_THEORY_std
            binoutobj.jf_dct_logtau = jf.dct.logtau
            binoutobj.jf_dct_logtau_THEORY_std = jf.dct.logtau_THEORY_std
            binoutobj.jf_dct_kappa = jf.dct.tau * jf.kappa_scale * 0.5
            binoutobj.jf_dct_kappa_THEORY_std = jf.dct.tau_THEORY_std * jf.kappa_scale * 0.5
            binoutobj.jf_dct_aic_Kmin = jf.dct.aic_Kmin
            binoutobj.jf_dct_Kmin_corrfactor = jf.dct.Kmin_corrfactor
        if not no_text_out:
            outfile_name = output + '.cepstral.dat'
            outarray = np.c_[jf.dct.logpsdK, jf.dct.logpsdK_THEORY_std, jf.dct.logtau, jf.dct.logtau_THEORY_std,
                             jf.dct.tau * jf.kappa_scale * 0.5, jf.dct.tau_THEORY_std * jf.kappa_scale * 0.5]
            outfile_header = 'ck  ck_std  L0(P*)  L0_std(P*)  kappa(P*)  kappa_std(P*)\n'
            np.savetxt(outfile_name, outarray, header=outfile_header)

        # plot cepstral log-PSD
        #ax = j.plot_periodogram(()  #PSD_FILTER_W=psd_filter_w)
        ax = jf.plot_periodogram()   #PSD_FILTER_W=psd_filter_w)
        jf.plot_cepstral_spectrum(axes=ax, label='cepstrum-filtered')
        ax[0].axvline(x=jf.Nyquist_f_THz, ls='--', c='r')
        ax[1].axvline(x=jf.Nyquist_f_THz, ls='--', c='r')
        #ax[0].set_xlim([0., 2.5*FSTAR_THZ])
        #ax[1].set_ylim([12,18])
        #ax[0].legend(['original', 'resampled', 'cepstrum-filtered'])
        #ax[1].legend(['original', 'resampled', 'cepstrum-filtered']);

        if binout:
            binoutobj.jf_dct_psd = jf.dct.psd
            binoutobj.jf_dct_logpsd = jf.dct.logpsd
        if not no_text_out:
            outfile_name = output + '.cepstrumfiltered_psd.dat'
            outarray = np.c_[jf.freqs_THz, jf.dct.psd, jf.dct.logpsd]
            outfile_header = 'freqs_THz  cepf_psd cepf_logpsd\n'
            np.savetxt(outfile_name, outarray, header=outfile_header)

        #conv_fact=open(output+'.kappa_scale_aicKmin.dat','w')
        #
        #if units=='metal':
        #    log.write_log 'kappa_scale (with DT_FS, can be used for gk-conversion)= {}'.format(tc.md.scale_kappa_METALtoSI(temperature,volume,DT_FS))
        #    conv_fact.write('{}\n'.format(tc.md.scale_kappa_METALtoSI(temperature,volume,DT_FS)))
        #elif units=='real':
        #    log.write_log 'kappa_scale (with DT_FS, can be used for gk-conversion) = {}'.format(tc.md.scale_kappa_REALtoSI(temperature,volume,DT_FS))
        #    conv_fact.write('{}\n'.format(tc.md.scale_kappa_REALtoSI(temperature,volume,DT_FS)))
        #conv_fact.write('{}\n'.format(jf.dct.aic_Kmin))
        #conv_fact.close()

    logfile.close()

    # write binary output
    if binout:
        if binout_old:
            binoutobj.write_old_binary(output)
        else:
            np.save(output, binoutobj)

    return 0


#################################


def plt_cepstral_conv(jf, pstar_max=None, k_SI_max=None, pstar_tick=None, kappa_tick=None):
    if pstar_max is None:
        pstar_max = (jf.dct.aic_Kmin + 1) * 2.5
    if k_SI_max is None:
        k_SI_max = jf.dct.tau[jf.dct.aic_Kmin] * jf.kappa_scale

    f, ax2 = plt.subplots(1, 1, figsize=(3.8, 2.3))
    ax2.axvline(x=jf.dct.aic_Kmin + 1, ls='--', c='k', dashes=(1.4, 0.6), zorder=-3)
    ax2.fill_between(np.arange(jf.dct.logtau.shape[0]) + 1,\
                     (jf.dct.tau - jf.dct.tau_THEORY_std) * jf.kappa_scale * 0.5,\
                     (jf.dct.tau + jf.dct.tau_THEORY_std) * jf.kappa_scale * 0.5,\
                     alpha=0.3, color=c[4], zorder=-3)
    ax2.plot(np.arange(jf.dct.logtau.shape[0]) + 1, jf.dct.tau * jf.kappa_scale * 0.5,\
             label=r'Cepstral method', marker='o', c=c[4], zorder=-3)
    ax2.set_xlabel(r'$P^*$')
    ax2.set_ylabel(r'$\kappa$ (W/mK)')
    ax2.set_xlim([0, pstar_max])
    ax2.set_ylim([0, k_SI_max])
    #ax2.grid()
    ax2.legend()

    if pstar_tick is None:
        dx1, dx2 = n_tick_in_range(0, pstar_max, 5)
    else:
        dx1 = pstar_tick
        dx2 = dx1 / 2
    if kappa_tick is None:
        dy1, dy2 = n_tick_in_range(0, k_SI_max, 5)
    else:
        dy1 = kappa_tick
        dy2 = dy1 / 2
    ax2.xaxis.set_major_locator(MultipleLocator(dx1))
    ax2.xaxis.set_minor_locator(MultipleLocator(dx2))
    ax2.yaxis.set_major_locator(MultipleLocator(dy1))
    ax2.yaxis.set_minor_locator(MultipleLocator(dy2))


def plt_other(jf, idx1, idx2, f_THz_max=None, k_SI_max=None, k_SI_min=None, k_tick=None, f_tick=None):
    if f_THz_max is None:
        idx_max = index_cumsum(np.abs(jf.fcospectrum[idx1][idx2]), 0.95)
        f_THz_max = jf.freqs_THz[idx_max]
    else:
        maxT = jf.freqs_THz[-1]
        if (maxT < f_THz_max):
            f_THz_max = maxT

    if k_SI_max is None:
        k_SI_max = np.max(
            np.abs(jf.fcospectrum[idx1][idx2])[:int(jf.freqs_THz.shape[0] * f_THz_max / jf.freqs_THz[-1])] *
            jf.kappa_scale * 0.5) * 1.3
    if k_SI_min is None:
        k_SI_min = -k_SI_max

    figure, ax = plt.subplots(1, 1, figsize=(3.8, 2.3))
    ax.plot(jf.freqs_THz, np.real(jf.fcospectrum[idx1][idx2]) * jf.kappa_scale * 0.5, c=c[3], lw=1.0, zorder=1)
    ax.plot(jf.freqs_THz, np.imag(jf.fcospectrum[idx1][idx2]) * jf.kappa_scale * 0.5, c=c[2], lw=1.0, zorder=1)

    ax.set_ylim([k_SI_min, k_SI_max])
    ax.set_xlim([0, f_THz_max])
    ax.set_xlabel(r'$\omega/2\pi$ (THz)')
    ax.set_ylabel(r'$S^{{{}{}}}$'.format(idx1, idx2))

    if f_tick is None:
        dx1, dx2 = n_tick_in_range(0, f_THz_max, 5)
    else:
        dx1 = f_tick
        dx2 = dx1 / 2
    if k_tick is None:
        dy1, dy2 = n_tick_in_range(0, k_SI_max, 5)
    else:
        dy1 = k_tick
        dy2 = dy1 / 2

    ax.xaxis.set_major_locator(MultipleLocator(dx1))
    ax.xaxis.set_minor_locator(MultipleLocator(dx2))
    ax.yaxis.set_major_locator(MultipleLocator(dy1))
    ax.yaxis.set_minor_locator(MultipleLocator(dy2))


def plt_psd(jf, j2=None, j2pl=None, f_THz_max=None, k_SI_max=None, k_tick=None, f_tick=None):
    if f_THz_max is None:
        idx_max = index_cumsum(jf.psd, 0.95)
        f_THz_max = jf.freqs_THz[idx_max]
    else:
        maxT = jf.freqs_THz[-1]
        if j2 is not None:
            if j2.freqs_THz[-1] > maxT:
                maxT = j2.freqs_THz[-1]
        if j2pl is not None:
            if j2pl.freqs_THz[-1] > maxT:
                maxT = j2pl.freqs_THz[-1]
        if maxT < f_THz_max:
            f_THz_max = maxT

    if k_SI_max is None:
        k_SI_max = np.max(
            jf.fpsd[:int(jf.freqs_THz.shape[0] * f_THz_max / jf.freqs_THz[-1])] * jf.kappa_scale * 0.5) * 1.3

    figure, ax = plt.subplots(1, 1, figsize=(3.8, 2.3))
    ax.plot(jf.freqs_THz, jf.psd * jf.kappa_scale * 0.5, lw=0.2, c='0.8', zorder=0)
    ax.plot(jf.freqs_THz, jf.fpsd * jf.kappa_scale * 0.5, c=c[0], zorder=2)
    if j2 is not None:
        plt.axvline(x=j2.Nyquist_f_THz, ls='--', c='k', dashes=(1.4, 0.6), zorder=3)
    if j2pl is not None:
        plt.plot(j2pl.freqs_THz, j2pl.dct.psd * j2pl.kappa_scale * 0.5, c=c[1], zorder=1)
    try:
        plt.plot(jf.freqs_THz, np.real(jf.fcospectrum[0][0]) * jf.kappa_scale * 0.5, c=c[3], lw=1.0, zorder=1)
    except:
        pass

    ax.set_ylim([0, k_SI_max])
    ax.set_xlim([0, f_THz_max])
    ax.set_xlabel(r'$\omega/2\pi$ (THz)')
    ax.set_ylabel(r'${}^{\ell}\hat{S}_{\,k}$ (W/mK)')

    if f_tick is None:
        dx1, dx2 = n_tick_in_range(0, f_THz_max, 5)
    else:
        dx1 = f_tick
        dx2 = dx1 / 2
    if k_tick is None:
        dy1, dy2 = n_tick_in_range(0, k_SI_max, 5)
    else:
        dy1 = k_tick
        dy2 = dy1 / 2

    ax.xaxis.set_major_locator(MultipleLocator(dx1))
    ax.xaxis.set_minor_locator(MultipleLocator(dx2))
    ax.yaxis.set_major_locator(MultipleLocator(dy1))
    ax.yaxis.set_minor_locator(MultipleLocator(dy2))


def n_tick_in_range(beg, end, n):
    import math
    size = end - beg
    n_cifre = math.floor(math.log(size / n, 10.0))
    delta = math.ceil((size / n) / 10**n_cifre) * 10**n_cifre
    return delta, delta / 2


def index_cumsum(arr, p):
    if (p > 1 or p < 0):
        raise ValueError('p must be between 0 and 1')
    arr_int = np.cumsum(arr)
    arr_int = arr_int / arr_int[-1]
    idx = 0
    while (arr_int[idx] < p):
        idx = idx + 1
    return idx


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

        self.jf_dct_logpsdK            = None
        self.jf_dct_logpsdK_THEORY_std = None
        self.jf_dct_logtau             = None
        self.jf_dct_logtau_THEORY_std  = None
        self.jf_dct_kappa              = None
        self.jf_dct_kappa_THEORY_std   = None
        self.jf_dct_psd                = None
        self.jf_dct_logpsd             = None
        self.jf_dct_aic_Kmin           = None
        self.jf_dct_Kmin_corrfactor    = None

        self.kappa_Kmin     = None
        self.kappa_Kmin_std = None
        self.cepstral_log   = None
        self.units          = None
        self.kappa_scale    = None
        self.TEMPERATURE    = None
        self.VOLUME         = None
        self.TSKIP          = None

    # yapf: enable
    def write_old_binary(self, output):
        """Write old binary format."""
        outarray = np.array([self.j_freqs_THz, self.j_fpsd, self.j_flogpsd, self.j_psd, self.j_logpsd])
        np.save(output + '.psd.npy', outarray)

        if self.j_cospectrum is not None:
            outarray = np.array([self.j_freqs_THz, self.j_cospectrum])
            np.save(output + '.cospectrum.npy', outarray)

        if self.j_fcospectrum is not None:
            outarray = np.array([self.j_freqs_THz, self.j_fcospectrum])
            np.save(output + '.cospectrum.filt.npy', outarray)

        outarray = np.array([self.jf_freqs_THz, self.jf_psd, self.jf_fpsd, self.jf_logpsd, self.jf_flogpsd])
        np.save(output + '.resampled_psd.npy', outarray)

        outarray = np.array([
            self.jf_dct_logpsdK, self.jf_dct_logpsdK_THEORY_std, self.jf_dct_logtau, self.jf_dct_logtau_THEORY_std,
            self.jf_dct_kappa, self.jf_dct_kappa_THEORY_std
        ])
        np.save(output + '.cepstral', outarray)

        outarray = np.array([self.jf_freqs_THz, self.jf_dct_psd, self.jf_dct_logpsd])
        np.save(output + '.cepstrumfiltered_psd', outarray)


if __name__ == '__main__':
    main()
