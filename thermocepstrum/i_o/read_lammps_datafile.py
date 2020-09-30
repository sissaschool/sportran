# -*- coding: utf-8 -*-

################################################################################
###
###   ReadLAMMPSdatafile
###
################################################################################

import numpy as np


def get_box(filename):
    """Return the box edges and volume from a LAMMPS Data file."""
    box = np.zeros((2, 3))
    with open(filename, 'r') as f:
        while True:
            line = f.readline().split()
            if 'xlo' in line:
                box[:, 0] = [float(line[0]), float(line[1])]
            if 'ylo' in line:
                box[:, 1] = [float(line[0]), float(line[1])]
            if 'zlo' in line:
                box[:, 2] = [float(line[0]), float(line[1])]
                break
    volume = np.prod(box[1, :] - box[0, :])
    return box, volume


## ** SHOULD ADD ADDITIONAL OF READ FUNCTIONS
