# -*- coding: utf-8 -*-

import numpy as np
from sportran.utils import log


class LowPassFilter(object):

    def __init__(self, *args, **kwargs):
        if (len(args) < 1):
            self.filtertype = kwargs.get('filtertype', None)
            if self.filtertype is None:
                raise ValueError('Not enough input arguments.')
        else:
            self.filtertype = args[0]
        if (self.filtertype == 'exp'):
            if (len(args) == 5):
                self.freqs = args[1]
                self.f0 = args[2]
                self.alpha = args[3]
                self.minatt = args[4]
            else:
                self.freqs = kwargs.get('freqs', None)
                self.f0 = kwargs.get('f0', None)
                self.alpha = kwargs.get('alpha', None)
                self.minatt = kwargs.get('minatt', None)
        return

    def __repr__(self):
        msg = '{} low-pass filter:\n'.format(self.filtertype) + \
              '  f0 = {}\n'.format(self.f0) + \
              '  alpha = {}\n'.format(self.alpha)
        return msg

    def compute_response(self, freqs=None):
        if freqs is not None:
            self.freqs = freqs
        if self.freqs is None:
            raise ValueError('freqs not set.')
        if (self.filtertype == 'exp'):
            self.response = self.exp_filter()
        self.logresponse = np.log(self.response)
        return

    def exp_filter(self):
        if self.f0 is None:
            raise ValueError('f0 not set.')
        if self.alpha is None:
            raise ValueError('alpha not set.')
        if self.minatt is None:
            self.minatt = 1.0e-3
        lf = np.zeros(self.freqs.size)
        for i in range(self.freqs.size):
            if (0. <= self.freqs[i] < 0.5):
                lf[i] = np.exp(-(self.freqs[i] / self.f0)**self.alpha) * (1.0 - self.minatt) + self.minatt
            elif (0.5 <= self.freqs[i] < 1.):
                lf[i] = np.exp(-(np.abs(self.freqs[i] - 1.) / self.f0)**self.alpha) * (1.0 - self.minatt) + self.minatt
            else:
                log.write_log('ERROR: frequency out of range!')
        return lf
