
from __future__ import absolute_import

import numpy as np
from scipy.fftpack import fft, ifft, dct
from scipy.signal  import periodogram, lfilter

import matplotlib.pyplot as plt
from statsmodels.tsa.api import acf

EULER_GAMMA = 0.57721566490153286060651209008240243104215933593992  # Euler-Mascheroni constant

