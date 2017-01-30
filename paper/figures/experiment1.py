""" TODO: migration not complete... from Experiment - Validation.ipynb """

# Standard library
from os.path import join, abspath

# Third-party
import corner
import astropy.units as u
import h5py
import matplotlib.pyplot as plt
import numpy as np

# Project
from thejoker.data import RVData
from thejoker.plot import plot_rv_curves
from thejoker.utils import quantity_from_hdf5

from figurehelpers import (plot_units, make_rv_curve_figure, apw_corner, mpl_style,
                           _truth_color, dpi)
plt.style.use(mpl_style)

def main():
    # experiment number
    e_number = 1
    e_name = 'validation'
    filename = abspath("../../cache/experiment{}.h5".format(e_number))

    # read the data
    with h5py.File(filename, 'r') as f:
        data = RVData.from_hdf5(f['data'])

        truth = dict()
        for k in f['truth']:
            truth[k] = quantity_from_hdf5(f['truth'], k)

    # read the samples from fixing the jitter
    with h5py.File(filename, 'r') as f:
        path = 'samples/fixed-jitter'
        samples1 = dict()
        for k in f[path]:
            samples1[k] = quantity_from_hdf5(f[path], k)

        path = 'samples/sample-jitter'
        samples2 = dict()
        for k in f[path]:
            samples2[k] = quantity_from_hdf5(f[path], k)

        # samples1[:,3] = (samples1[:,3] - 360*((data.t_offset/pars1._P) % 1.)) % (360) # HACK

    print("{} samples survived (a)".format(samples1['P'].shape[0]))
    print("{} samples survived (b)".format(samples2['P'].shape[0]))

    fig = make_rv_curve_figure(data, [pars1, pars2], truth_pars=truth_pars,
                               units=samples_units, rv_lim=(21, 67))
    fig.axes[0].set_title("Experiment {}".format(e_number), fontsize=22)

    fig.axes[0].text(55700, 60, "(a)", fontsize=18)
    fig.axes[1].text(55700, 60, "(b)", fontsize=18)

    fig.tight_layout()
    fig.savefig(join('{}-rv-curves.pdf'.format(e_name)), dpi=dpi)

if __name__ == '__main__':
    main()
