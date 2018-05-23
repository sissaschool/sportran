

#on Cineca cluster Marconi please run:
#module load autoload scipy/0.18.1--python--2.7.12

import sys
sys.path.append('/marconi_work/Sis18_baroni')
import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import thermocepstrum as tc

usage=\
"""\
usage:
  {} input_data output units Volume DT_FS F_CUT_THz corr_factor heat_flux_key [other_fluxes ...]

do the cepstral analysis. It outputs some results in the stdout, and some pictures in a pdf called [output].pdf.
The input must be a column-formatted text file, with a header in the same format of LAMMPS. The name of the lammps compute can start with c_ and end with [#some_number], the code will understand automatically if it is a vector, and it will read automatically the correct number of components. In this case you have to provide only the name of the compute (the header without the initial c_ and the final [*]).
There must be a column with the temperature with header 'Temp'. The code computes automatically the mean temperature. The output is in SI units.
You must provide the name of the heat flux compute in the last arguments. The code will understand from the number of the headers that you provide the number of the components, and select the correct procedure. Note that the output is the same with any number of components. If you have a lots of components, note that you may want to use more than 3 independent processes -- see theory.
units can be metal or real (see LAMMPS documentation at http://lammps.sandia.gov/doc/units.html )

Output files (content of the columns):
  [output].periodogram:
      freqs in THz, original periodogram, original log(periodogram)
  [output].resampled_periodogram:
      freqs in THz, resampled periodogram, resampled log(periodogram)
  [output].cepstral:
      jf.dct.logpsdK,jf.dct.logpsdK_THEORY_std,jf.dct.logtau,jf.dct.logtau_THEORY_std,jf.dct.tau*jf.kappa_scale * 0.5,jf.dct.tau_THEORY_std*jf.kappa_scale * 0.5
  [output].cepstrum_filtered_periodogram:
      freqs in THz, cepstral filtered periodogram, cepstral filtered log(periodogram) 

Enjoy it!



This software was written by Loris Ercole and extended by Riccardo Bertossa to handle the multicomponent stuff, at SISSA, via Bonomea, 265 - 34136 Trieste ITALY.

Please see https://arxiv.org/abs/1706.01381 (Loris Ercole, Aris Marcolongo, Stefano Baroni) and https://arxiv.org/abs/1802.08006 (Stefano Baroni, Riccardo Bertossa, Loris Ercole, Federico Grasselli, Aris Marcolongo)
""".format(sys.argv[0])

if len(sys.argv) < 9:
     print usage
     exit(-1)

volume=float(sys.argv[4])
DT_FS=float(sys.argv[5])
FSTAR_THZ=float(sys.argv[6])
corr_factor=float(sys.argv[7])
units=sys.argv[3]
output=sys.argv[2]

if units != 'real' and units != 'metal':
   print usage
   print 'units must be metal or real'
   exit(-1)


selected_keys= sys.argv[8:]
selected_keys.append('Temp')


jfile = tc.i_o.TableFile(sys.argv[1], group_vectors=True)
jfile.read_datalines(start_step=0, NSTEPS=0, select_ckeys=selected_keys)

temperature=np.mean(jfile.data['Temp'])
print "mean temperature: {}K".format(temperature)
currents=[]
if len(selected_keys[:-1]) == 1:
    currents=jfile.data[selected_keys[0]]
else:
    for key in selected_keys[:-1]:
        currents.append(jfile.data[key])

command_file=open(output+'.command','w')
command_file.write(' '.join(sys.argv))
command_file.close()

j = tc.heatcurrent.HeatCurrent(currents,units,DT_FS,temperature,volume)

with PdfPages(output+"_all.pdf") as pdf:
   j.plot_periodogram(PSD_FILTER_W=0.5)
   np.savetxt(output+".periodogram",np.c_[j.freqs_THz,j.fpsd,j.flogpsd])
   pdf.attach_note('Nyquist frequency = {} THz'.format(j.Nyquist_f_THz), positionRect=[0,0 , 100, 50])
   pdf.savefig()
   plt.close()

   jf, ax = tc.heatcurrent.resample_current(j, fstar_THz=FSTAR_THZ, plot=True, freq_units='thz')
   plt.xlim([0, 2.5*FSTAR_THZ])
   ax[1].set_xlim()
   np.savetxt(output+".resampled_periodogram",np.c_[jf.freqs_THz,jf.fpsd,jf.flogpsd])
   pdf.savefig()
   plt.close()

   ax = jf.plot_periodogram(PSD_FILTER_W=0.1)
   pdf.savefig()
   plt.close()

   jf.cepstral_analysis(aic_type='aic', Kmin_corrfactor=corr_factor)
   ax = jf.plot_ck()
   ax.set_xlim([0, 5*jf.dct.aic_Kmin])
#  ax.set_ylim([-0.5, 0.5])
   ax.grid();
   pdf.savefig()
   plt.close()



   # kappa as a function of cutoff K
   ax = jf.plot_L0_Pstar()
   ax.set_xlim([0,10*jf.dct.aic_Kmin])
   pdf.savefig()
   plt.close()
 
   ax = jf.plot_kappa_Pstar()
   ax.set_xlim([0,10*jf.dct.aic_Kmin])
   pdf.savefig()
   plt.close()
   np.savetxt(output+".cepstral",np.c_[jf.dct.logpsdK,jf.dct.logpsdK_THEORY_std,jf.dct.logtau,jf.dct.logtau_THEORY_std,jf.dct.tau*jf.kappa_scale * 0.5,jf.dct.tau_THEORY_std*jf.kappa_scale * 0.5])

   # filtered log-PSD
   ax = j.plot_periodogram(0.5)
   ax = jf.plot_periodogram(0.5, axes=ax)
   ax[0].plot(jf.freqs_THz, jf.dct.psd,    c='b', lw=2, label='filtered')
   ax[1].plot(jf.freqs_THz, jf.dct.logpsd, c='b', lw=2, label='filtered')
   ax[0].axvline(x = jf.Nyquist_f_THz, ls='--', c='r')
   ax[1].axvline(x = jf.Nyquist_f_THz, ls='--', c='r')
   plt.xlim([0., 2.5*FSTAR_THZ])
#   ax[1].set_ylim([12,18])
   ax[0].legend(['original', 'resampled', 'cepstrum-filtered'])
   ax[1].legend(['original', 'resampled', 'cepstrum-filtered']);
   
   np.savetxt(output+".cepstrum_filtered_periodogram",np.c_[jf.freqs_THz,jf.dct.psd,jf.dct.logpsd])

   print '  L_0* = {:15g} +/- {:10f}'.format(jf.dct.logtau_Kmin, jf.dct.logtau_std_Kmin)
   print '  S_0* = {:15f} +/- {:10f}'.format(jf.dct.tau_Kmin, jf.dct.tau_std_Kmin)
   print '-------------------------------------------------'
   print '  kappa* = {:15f} +/- {:10f}  W/mK'.format(jf.kappa_Kmin, jf.kappa_Kmin_std)
   print '-------------------------------------------------'
   
   print 'K of AIC_min = {:d}'.format(jf.dct.aic_Kmin)

