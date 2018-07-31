#!/usr/bin/python

# Needed modules:
#  numpy
#  scipy
#  matplotlib
#  argparse

import argparse
import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['figure.figsize'] = (16, 9)
plt.style.reload_library()
plt.style.use('grafici_belli.mplstyle')
plt.rc('text', usetex=True)
c = plt.rcParams['axes.prop_cycle'].by_key()['color'] 
from matplotlib.ticker import MultipleLocator


from sys import path, argv
tc_path = argv[0][:argv[0].rfind('/')]
path.append(tc_path[:tc_path.rfind('/')])
import thermocepstrum as tc
import math

def main():
   """
--------------------------------------------------------------------------------
  *** THERMOCEPSTRUM ***  command line interface  (beta)
--------------------------------------------------------------------------------
This script performs the cepstral analysis. It outputs some results in the stdout and log file, and plots in pdf format.

INPUT:
For now the input must be a column-formatted text file, with a header in the same format of LAMMPS. The name of the lammps compute can start with c_ and end with [#some_number], the code will recognize vectors, and will read automatically all the components.
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
 - for the multi-component analysis:  Baroni, Bertossa, Ercole, Grasselli, Marcolongo, https://arxiv.org/abs/1802.08006

https://github.com/lorisercole/thermocepstrum
Contact: lercole@sissa.it
"""

   parser = argparse.ArgumentParser(description=main.__doc__, epilog=_epilog, formatter_class=argparse.RawTextHelpFormatter)
   parser.add_argument( 'inputfile', type=str, help='input file to read (default format: Table)' )
   parser.add_argument( '-V', '--volume', type=float, help='Volume of the cell (Angstrom)' )
   parser.add_argument( '-t', '--timestep', type=float, required=True, help='Time step of the printed data (fs)' )
   parser.add_argument( '-k', '--heatfluxkey', type=str, required=True, help='Name of the column keyword that identifies the heat flux' )
   parser.add_argument( '-N', '--nsteps', type=int, default=0, help='Number of steps to read (default: 0=all)' )
   parser.add_argument( '-S', '--start-step', type=int, default=0, help='The first step to read (default: 0=first)' )
   parser.add_argument( '--input-format', default='table', type=str, choices=['table','dict'], help='format of the input file' )
   parser.add_argument( '--cindex', nargs='*', type=int, help='column indexes of the heatflux to read (0,1,2,...)' )

   outarg = parser.add_mutually_exclusive_group()
   outarg.add_argument( '-o', '--output', type=str, default='output', help='prefix of the output files' )
   outarg.add_argument( '-O', '--bin-output', type=str, help='prefix of the output files (use binary file)' )

   parser.add_argument( '-u', '--units', type=str, default='metal', choices=['metal', 'real'], help='LAMMPS units (default: metal)' )
   parser.add_argument( '-T', '--temperature', type=float, help='average Temperature (K). If not set it will be read from file' )

   parser.add_argument( '-r', '--resample', action='store_true', help='resample the time series (you should define --TSKIP or --FSTAR' )
   resamplearg = parser.add_mutually_exclusive_group()
   resamplearg.add_argument( '--TSKIP', type=int, help='resampling time period (steps)' )
   resamplearg.add_argument( '--FSTAR', type=float, help='resampling target Nyquist frequency (THz)' )
   parser.add_argument( '-c', '--corr-factor', type=float, default=1.0, help='correction factor to the AIC' )
   parser.add_argument( '-j', '--add-currents', type=str, default=[], action='append', help='additional current for multi-component fluids' )

   parser.add_argument( '-w', '--psd-filterw', type=float, help='plot - periodogram - filter window width (THz)' )
   args = parser.parse_args()

   inputfile = args.inputfile
   volume = args.volume
   DT_FS = args.timestep
   j1_key = args.heatfluxkey
   NSTEPS = args.nsteps
   START_STEP = args.start_step
   input_format = args.input_format
   jindex = args.cindex
   
   if args.bin_output is not None:
      binout = True
      output = args.bin_output
   else:
      binout = False
      output = args.output

   units = args.units
   temperature = args.temperature
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
   if ((units != 'real') and (units != 'metal')):
      raise ValueError('units must be metal or real')
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
      raise ValueError('Use flag -r to resample. TSKIP will be ignored.')
   elif (FSTAR is not None):
      raise ValueError('Use flag -r to resample. FSTAR will be ignored.')
   if (corr_factor <= 0.):
      raise ValueError('the correction factor must be positive')

   ncurrents = len(j2_keys) + 1
   
   logfile = open(output + '.log', 'w')
   logfile.write('Command:\n ' + ' '.join(argv) + '\n\n')

   selected_keys = [j1_key]
   selected_keys.extend(j2_keys)

   jdata = None
   if (input_format == 'table'):
      if temperature is None:
         selected_keys.append('Temp')
      jfile = tc.i_o.TableFile(inputfile, group_vectors=True)
      jfile.read_datalines(start_step=START_STEP, NSTEPS=NSTEPS, select_ckeys=selected_keys)
      jdata = jfile.data
   elif (input_format == 'dict'):
      jdata = np.load(inputfile)
   else:
      raise NotImplemented('input format not implemented.')

   if NSTEPS==0:
      NSTEPS=jdata[jdata.keys()[0]].shape[0]

   if temperature is None:
      if 'Temp' in jdata:
         temperature = np.mean(jdata['Temp'])
         temperature_std = np.std(jdata['Temp'])
         selected_keys.remove('Temp')
         print ' Mean Temperature (computed):  {} K  +/-  {}'.format(temperature, temperature_std)
         logfile.write(' Mean Temperature (computed):  {} K  +/-  {}\n'.format(temperature, temperature_std))
      elif 'Temp_ave' in jdata:
         temperature = jdata['Temp_ave']
         if 'Temp_std' in jdata:
            temperature_std = jdata['Temp_std']
            print ' Mean Temperature (file):      {} K  +/-  {}'.format(temperature, temperature_std)
            logfile.write(' Mean Temperature (file):      {} K  +/-  {}\n'.format(temperature, temperature_std))
         else:
            print ' Mean Temperature (file):      {} K'.format(temperature)
            logfile.write(' Mean Temperature (file):      {} K\n'.format(temperature))
      else:
         raise RuntimeError('No Temp key found')
   else:
      print ' Mean Temperature (input):  {} K'.format(temperature)
      logfile.write(' Mean Temperature (input):  {} K\n'.format(temperature))

   if volume is None:
      if 'Volume' in jdata:
         volume = jdata['Volume']
         print ' Volume (file):    {} A^3'.format(volume)
         logfile.write(' Volume (file):    {} A^3\n'.format(volume))
      else:
         raise RuntimeError('No Volume key found')
   else:
       print ' Volume (input):  {} A^3'.format(volume)
       logfile.write(' Volume (input):  {} A^3\n'.format(volume))

   print ' Time step (input):  {} fs'.format(DT_FS)
   logfile.write(' Time step (input):  {} fs\n'.format(DT_FS))

   print selected_keys, jindex
   if jindex is None:
      currents = np.array([jdata[key][START_STEP:(START_STEP+NSTEPS),:] for key in selected_keys])
   else:
      currents = np.array([jdata[key][START_STEP:(START_STEP+NSTEPS),jindex] for key in selected_keys])
   print '  currents shape is {}'.format(currents.shape)
   logfile.write('  currents shape is {}\n'.format(currents.shape))
   print 'snippet:'
   print currents

   # create HeatCurrent object
   j = tc.heatcurrent.HeatCurrent(currents, units, DT_FS, temperature, volume, psd_filter_w)

   print ' Number of currents = {}'.format(ncurrents)
   logfile.write(' Number of currrents = {}\n'.format(ncurrents))
   print ' Number of components = {}'.format(j.N_COMPONENTS)
   logfile.write(' Number of components = {}\n'.format(j.N_COMPONENTS))
   print ' kappa_scale = {}'.format(j.kappa_scale)
   logfile.write(' kappa_scale = {}\n'.format(j.kappa_scale))
   print ' Nyquist_f   = {}  THz'.format(j.Nyquist_f_THz)
   logfile.write(' Nyquist_f   = {}  THz\n'.format(j.Nyquist_f_THz))

   with PdfPages(output + ".plots.pdf") as pdf:

      # plot periodogram
      j.plot_periodogram()  #PSD_FILTER_W=psd_filter_w)
      pdf.attach_note('Nyquist frequency = {} THz'.format(j.Nyquist_f_THz), positionRect=[0, 0 , 5, 1])
      pdf.savefig()
      plt.close()

      if binout:
         outarray = np.array([j.freqs_THz, j.fpsd, j.flogpsd, j.psd, j.logpsd])
         np.save(output + ".psd.npy", outarray, allow_pickle=False)
      else:
         outfile = open(output + '.psd.dat', 'w')
         outarray = np.c_[j.freqs_THz, j.psd, j.fpsd, j.logpsd, j.flogpsd]
         outfile.write('freqs_THz  psd  fpsd  logpsd  flogpsd\n')
         np.savetxt(outfile, outarray)
         outfile.close()
  
      # resample and plot
      if resample:
         if TSKIP is not None:
            jf, ax = tc.heatcurrent.resample_current(j, TSKIP=TSKIP, plot=True, PSD_FILTER_W=psd_filter_w)
         else:
            jf, ax = tc.heatcurrent.resample_current(j, fstar_THz=FSTAR, plot=True, PSD_FILTER_W=psd_filter_w)
         ax[0].set_xlim([0, 2.5*FSTAR])
         pdf.savefig()
         plt.close()
         logfile.write(jf.resample_log)

         # plot resampled periodogram
         ax = jf.plot_periodogram()  #PSD_FILTER_W=psd_filter_w)
         pdf.savefig()
         plt.close()

         if binout:
            outarray = np.array([jf.freqs_THz, jf.psd, jf.fpsd, jf.logpsd, jf.flogpsd])
            np.save(output + ".resampled_psd.npy", outarray, allow_pickle=False)
         else:
            outfile = open(output + '.resampled_psd.dat', 'w')
            outarray = np.c_[jf.freqs_THz, jf.psd, jf.fpsd, jf.logpsd, jf.flogpsd]
            outfile.write('freqs_THz  psd  fpsd  logpsd  flogpsd\n')
            np.savetxt(outfile, outarray)
            outfile.close()
      else:
         jf = j

      plt_psd(j,jf)
      pdf.savefig()
      plt.close()
      plt_psd(jf)
      pdf.savefig()
      plt.close()

      # cepstral analysis
      jf.cepstral_analysis(aic_type='aic', Kmin_corrfactor=corr_factor)
      logfile.write(jf.cepstral_log)

      plt_psd(jf,jf,jf)
      pdf.savefig()
      plt.close()
      

      # plot cepstral coefficients
      ax = jf.plot_ck()
      ax.set_xlim([0, 5*jf.dct.aic_Kmin])
   #  ax.set_ylim([-0.5, 0.5])
      ax.grid();
      pdf.savefig()
      plt.close()

      # plot L0(Pstar)
      ax = jf.plot_L0_Pstar()
      ax.set_xlim([0, 10*jf.dct.aic_Kmin])
      pdf.savefig()
      plt.close()
   
      # plot kappa(Pstar)
#      ax = jf.plot_kappa_Pstar()
#      ax.set_xlim([0, 10*jf.dct.aic_Kmin])
      plt_cepstral_conv(jf)
      pdf.savefig()
      plt.close()

      if binout:
         outarray = np.array([jf.dct.logpsdK, jf.dct.logpsdK_THEORY_std, jf.dct.logtau, jf.dct.logtau_THEORY_std, jf.dct.tau*jf.kappa_scale*0.5, jf.dct.tau_THEORY_std*jf.kappa_scale*0.5])
         np.save(output + '.cepstral', outarray, allow_pickle=False)
      else:
         outfile = open(output + '.cepstral.dat', 'w')
         outarray = np.c_[jf.dct.logpsdK, jf.dct.logpsdK_THEORY_std, jf.dct.logtau, jf.dct.logtau_THEORY_std, jf.dct.tau*jf.kappa_scale*0.5, jf.dct.tau_THEORY_std*jf.kappa_scale*0.5]
         outfile.write('ck  ck_std  L0(P*)  L0_std(P*)  kappa(P*)  kappa_std(P*)\n')
         np.savetxt(outfile, outarray)
         outfile.close()
   
      # plot cepstral log-PSD
      #ax = j.plot_periodogram(()  #PSD_FILTER_W=psd_filter_w)
      ax = jf.plot_periodogram()  #PSD_FILTER_W=psd_filter_w)
      jf.plot_cepstral_spectrum(axes=ax, label='cepstrum-filtered')
      ax[0].axvline(x = jf.Nyquist_f_THz, ls='--', c='r')
      ax[1].axvline(x = jf.Nyquist_f_THz, ls='--', c='r')
      #ax[0].set_xlim([0., 2.5*FSTAR_THZ])
      #ax[1].set_ylim([12,18])
      #ax[0].legend(['original', 'resampled', 'cepstrum-filtered'])
      #ax[1].legend(['original', 'resampled', 'cepstrum-filtered']);
   
      if binout:
         outarray = np.array([jf.freqs_THz, jf.dct.psd, jf.dct.logpsd])
         np.save(output+".cepstrumfiltered_psd", outarray, allow_pickle=False)
      else:
         outfile = open(output + '.cepstrumfiltered_psd.dat', 'w')
         outarray = np.c_[jf.freqs_THz, jf.dct.psd, jf.dct.logpsd]
         outfile.write('freqs_THz  cepf_psd cepf_logpsd\n')
         np.savetxt(outfile, outarray)

      conv_fact=open(output+'.kappa_scale_aicKmin.dat','w')
      
      if units=='metal':
          print 'kappa_scale (with DT_FS, can be used for gk-conversion)= {}'.format(tc.md.scale_kappa_METALtoSI(temperature,volume,DT_FS))
          conv_fact.write('{}\n'.format(tc.md.scale_kappa_METALtoSI(temperature,volume,DT_FS)))
      elif units=='real':
          print 'kappa_scale (with DT_FS, can be used for gk-conversion) = {}'.format(tc.md.scale_kappa_REALtoSI(temperature,volume,DT_FS))
          conv_fact.write('{}\n'.format(tc.md.scale_kappa_REALtoSI(temperature,volume,DT_FS)))
      conv_fact.write('{}\n'.format(jf.dct.aic_Kmin))      

      conv_fact.close()
      


   logfile.close()
   return 0


def plt_cepstral_conv(jf,pstar_max=None, k_SI_max=None):
    if pstar_max==None:
       pstar_max=(jf.dct.aic_Kmin+1)*1.618
    if k_SI_max==None:
       k_SI_max=jf.dct.tau[jf.dct.aic_Kmin]*jf.kappa_scale

    f, (ax2) = plt.subplots(1,1,figsize=(3.8,2.3))
    ax2.axvline(x=jf.dct.aic_Kmin+1, ls='--', c='k', dashes=(1.4,0.6), zorder=-3)
    ax2.fill_between(np.arange(jf.dct.logtau.shape[0])+1,\
                     (jf.dct.tau-jf.dct.tau_THEORY_std)*jf.kappa_scale*.5,\
                     (jf.dct.tau+jf.dct.tau_THEORY_std)*jf.kappa_scale*.5,alpha=.3,color=c[4],zorder=-3)#'#3c8da8')
    ax2.plot(np.arange(jf.dct.logtau.shape[0])+1, jf.dct.tau*jf.kappa_scale*.5,\
             label=r'Cepstral method',marker='o',c=c[4],zorder=-3)#'#3c8da8')
    ax2.set_xlabel('$P^*$')
    ax2.set_ylabel('$\kappa$ (W/mK)')
    ax2.set_xlim([0,pstar_max])
    ax2.set_ylim([0, k_SI_max])
#    ax2.grid()
    ax2.legend()
    dx1,dx2=n_tick_in_range(0,pstar_max,5)
    dy1,dy2=n_tick_in_range(0,k_SI_max,5)

    ax2.xaxis.set_major_locator(MultipleLocator(dx1))
    ax2.xaxis.set_minor_locator(MultipleLocator(dx2))
    ax2.yaxis.set_major_locator(MultipleLocator(dy1))
    ax2.yaxis.set_minor_locator(MultipleLocator(dy2))

def plt_psd(jf,j2=None,j2pl=None,f_THz_max=None, k_SI_max=None):

    if f_THz_max==None:
       idx_max=index_cumsum(jf.psd,0.95)
       f_THz_max=jf.freqs_THz[idx_max]

    if k_SI_max==None:
       k_SI_max=np.max(jf.fpsd[:int(jf.freqs_THz.shape[0]*f_THz_max/jf.freqs_THz[-1])]*jf.kappa_scale*.5) *1.3

    plt.figure(figsize=(3.8,2.3))
    plt.plot(jf.freqs_THz,jf.psd*jf.kappa_scale*.5,lw=0.2,c='0.8')
    plt.plot(jf.freqs_THz,jf.fpsd*jf.kappa_scale*.5,c=c[0])
    if j2 != None:
        plt.axvline(x=j2.Nyquist_f_THz,ls='--', c='k', dashes=(1.4,0.6), zorder=3)
    if j2pl != None:
       plt.plot(j2pl.freqs_THz,j2pl.dct.psd*j2pl.kappa_scale*.5,c=c[1])

    plt.ylim([0,k_SI_max])
    plt.xlim([0,f_THz_max])
    plt.xlabel('$\omega/2\pi$ (THz)')
    plt.ylabel('${}^{\ell}\hat{\underline{S}}_{\,k}$ (W/mK)')
    
    dx1,dx2=n_tick_in_range(0,f_THz_max,5)
    dy1,dy2=n_tick_in_range(0,k_SI_max,5)

    plt.axes().xaxis.set_major_locator(MultipleLocator(dx1))
    plt.axes().xaxis.set_minor_locator(MultipleLocator(dx2))
    plt.axes().yaxis.set_major_locator(MultipleLocator(dy1))
    plt.axes().yaxis.set_minor_locator(MultipleLocator(dy2))

def n_tick_in_range(beg,end,n):
    size=end-beg
    n_cifre=math.floor(math.log(size/n,10.0))
    delta=math.ceil((size/n)/10**n_cifre)*10**n_cifre
    return delta,delta/2

def index_cumsum(arr,p):
    if (p>1 or p<0):
        raise ValueError('p must be between 0 and 1')
    arr_int=np.cumsum(arr)
    arr_int=arr_int/arr_int[-1]
    idx=0
    while arr_int[idx]<p:
        idx=idx+1
    return idx

if __name__ == "__main__":
   main()
