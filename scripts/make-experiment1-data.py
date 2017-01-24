"""

Generate simulated data for Experiment 1 that match all of our assumptions:
    - A pair of stars, one of which is observed spectroscopically, with radial velocity variations
      that are described purely through the 2-body problem
    - Gaussian uncertainties on the RV measurements

Here, we save the true uncertainties used to generate the data.

"""

# Standard library
import os

# Third-party
import astropy.units as u
import h5py
import matplotlib
matplotlib.use('agg')
import numpy as np

# Project
from thejoker.celestialmechanics import SimulatedRVOrbit
from thejoker.data import RVData
from thejoker.utils import quantity_to_hdf5

def main():

    kms = u.km/u.s
    if not os.path.exists(os.path.abspath("../scripts")):
        raise RuntimeError("Script must be run from within the scripts directory.")

    cache_path = os.path.abspath("../cache")
    if not os.path.exists(cache_path):
        os.makedirs(cache_path)

    # high-eccentricity orbit with reasonable or randomly chosen parameters
    opars = dict(P=103.71*u.day, K=8.134*kms, ecc=0.313,
                 omega=np.random.uniform(0, 2*np.pi)*u.rad,
                 phi0=np.random.uniform(0, 2*np.pi)*u.rad,
                 v0=np.random.normal(0, 30) * kms)
    orbit = SimulatedRVOrbit(**opars)

    n_obs = 5 # number of observations

    bmjd = np.random.uniform(0, 3*365, size=n_obs) + 55555. # 3 year survey
    bmjd.sort()
    rv = orbit.generate_rv_curve(bmjd)
    rv_err = np.random.uniform(100, 200, size=n_obs) * u.m/u.s # apogee-like
    rv = np.random.normal(rv.to(kms).value, rv_err.to(kms).value) * kms

    data = RVData(t=bmjd, rv=rv, stddev=rv_err)
    with h5py.File(os.path.join(cache_path, "experiment1.h5"), "w") as f:
        g = f.create_group('data')
        data.to_hdf5(g)

        g = f.create_group('truth')
        for k in opars:
            quantity_to_hdf5(g, k, opars[k])

if __name__ == '__main__':
    from argparse import ArgumentParser

    # Define parser object
    parser = ArgumentParser(description="")

    parser.add_argument("-s", "--seed", dest="seed", default=None,
                        type=int, help="Random number seed.")

    args = parser.parse_args()

    if args.seed is not None:
        np.random.seed(args.seed)

    main()
