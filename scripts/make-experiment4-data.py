"""

Take APOGEE data for the target 2M00110648+6609349 (used in Troup et al.) and
store it in the same format as the other experiment data.

"""

# Standard library
import os

# Third-party
from astropy.io import fits
import astropy.units as u
import h5py
import matplotlib
matplotlib.use('agg')
import numpy as np

# Project
from thejoker.data import RVData

def main():

    if not os.path.exists(os.path.abspath("../scripts")):
        raise RuntimeError("Script must be run from within the scripts directory.")

    root_path = os.path.abspath("..")
    cache_path = os.path.join(root_path, "cache")
    if not os.path.exists(cache_path):
        os.makedirs(cache_path)

    fits_data = fits.getdata(os.path.join(root_path, "data/2M00110648+6609349_visits.fits"))

    # read sub-set data file:
    good_idx = (np.isfinite(fits_data['MJD']) &
                np.isfinite(fits_data['VHELIO']) &
                np.isfinite(fits_data['VRELERR']))
    fits_data = fits_data[good_idx]

    data = RVData(t=fits_data['MJD'],
                  rv=fits_data['VHELIO']*u.km/u.s,
                  stddev=fits_data['VRELERR']*u.km/u.s)
    with h5py.File(os.path.join(cache_path, "experiment4.h5"), "w") as f:
        g = f.create_group('data')
        data.to_hdf5(g)

if __name__ == '__main__':
    # from argparse import ArgumentParser

    # Define parser object
    # parser = ArgumentParser(description="")

    # parser.add_argument("-s", "--seed", dest="seed", default=None,
    #                     type=int, help="Random number seed.")

    # args = parser.parse_args()

    main()
