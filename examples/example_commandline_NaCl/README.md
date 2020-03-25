# Command line example

In this example we perform the same analysis that is perfomed in the `example_cepstrum_doublecomp_NaCl.ipynb` jupyter notebook in a much straightforward way.
Symply run, after installing the package, the following command (see `command.sh`)

```
thermocepstrum-analysis ../data/NaCl.dat -k flux -j 'vcm[1]' -t 5.0 -V 65013.301261 -w 0.1 --FSTAR 14.0 -r
```
If the package is not installed you can run this instead:
`../../thermocepstrum/analysis.py ../data/NaCl.dat -k flux -j 'vcm[1]' -t 5.0 -V 65013.301261 -w 0.1 --FSTAR 14.0 -r`

The options have the following meaning:
* `-k flux` use the columns with header `flux` as the energy flux
* `-j 'vcm[1]` use the columns with header `vcm[1]` as the convective flux
* `-t 5.0` set the timestep to 5.0 fs
* `-V 65013.301261` set the volume of the system to 65013.301261 A^3
* `-w 0.1` the width of the moving average filter used only to visualize the spectrum is 0.1 THz
* `--FSTAR 14.0` set the $f^*$ cutoff frequency to 14.0 THz
* `-r` resample according to the value of $f^*$ specified with `--FSTAR`

The output of the program in the terminal is:
```
 Input file (table):      ../data/NaCl.dat
 Units:      metal
 Time step:      5.0 fs
Temp    c_flux[1] c_flux[2] c_flux[3] c_vcm[1][1] c_vcm[1][2] c_vcm[1][3]
 #####################################
  all_ckeys =  {'Temp': [0], 'flux': array([1, 2, 3]), 'vcm[1]': array([4, 5, 6])}
 #####################################
Data length =  20000
  ckey =  {'flux': array([1, 2, 3]), 'vcm[1]': array([4, 5, 6]), 'Temp': [0]}
  ( 20000 ) steps read.
DONE.  Elapsed time:  0.25597500801086426 seconds
 Mean Temperature (computed):  1399.3477811999999 K  +/-  19.318785820942594
 Volume (input):  65013.301261 A^3
 Time step (input):  5.0 fs
['flux', 'vcm[1]'] None
  currents shape is (2, 20000, 3)
snippet:
[[[ 2.5086549e+02  2.0619423e+01  2.0011500e+02]
  [ 1.9622265e+02  8.2667342e+01  2.8433250e+02]
  [ 1.2639441e+02  1.6075472e+02  3.4036769e+02]
  ...
  [ 1.7991856e+02  1.8612706e+01 -1.3265623e+02]
  [ 2.0471193e+02 -4.6643315e-01 -2.0401650e+02]
  [ 2.4123318e+02 -1.8295461e+01 -2.5246475e+02]]

 [[-1.5991832e-01 -7.1370426e-02  2.0687917e-02]
  [-1.3755206e-01 -7.1002931e-02 -1.1279876e-02]
  [-1.0615044e-01 -6.2381243e-02 -4.1568120e-02]
  ...
  [-9.1939899e-02 -8.4778292e-02  6.0011385e-02]
  [-1.3384949e-01 -1.1474530e-01  8.9323546e-02]
  [-1.8385053e-01 -1.3693430e-01  1.1434060e-01]]]
Using multicomponent code.
 Number of currents = 2
 Number of components = 3
 kappa_scale = 1.4604390788939313e-07
 Nyquist_f   = 100.0  THz
Using multicomponent code.
-----------------------------------------------------
  RESAMPLE TIME SERIES
-----------------------------------------------------
 Original Nyquist freq  f_Ny =     100.00000 THz
 Resampling freq          f* =      14.28571 THz
 Sampling time         TSKIP =             7 steps
                             =        35.000 fs
 Original  n. of frequencies =         10001
 Resampled n. of frequencies =          1429
 PSD      @cutoff  (pre-filter) = 443152.37265
                  (post-filter) = 564877.86516
 log(PSD) @cutoff  (pre-filter) =     12.89638
                  (post-filter) =     13.05597
 min(PSD)          (pre-filter) =      0.31536
 min(PSD)         (post-filter) =  22166.11934
 % of original PSD Power f<f* (pre-filter)  = 96.678574
-----------------------------------------------------

-----------------------------------------------------
  CEPSTRAL ANALYSIS
-----------------------------------------------------
  AIC_Kmin  = 3  (P* = 4, corr_factor = 1.000000)
  L_0*   =          15.158757 +/-   0.056227
  S_0*   =     6824108.702608 +/- 383697.095268
-----------------------------------------------------
  kappa* =           0.498310 +/-   0.028018  W/mK
-----------------------------------------------------
```

The program outputs raw data and some pdf plots.\
In this example the output files are called `output.*`.
