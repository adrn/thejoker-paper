"""

Generate simulated data for Experiment 5:
    - Read in data from Experiment 1
    - Move 2nd data point through 1 cycle in period to show how the sampling depends on when, in
      phase, it is observed

"""

# Standard library
import os

# Third-party
import astropy.units as u
import h5py
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

    opars = dict(P=127.31*u.day, K=8.996045*kms, ecc=0.213,
                 omega=137.234*u.degree,
                 phi0=36.231*u.degree,
                 v0=17.643*kms)
    orbit = SimulatedRVOrbit(**opars)

    EPOCH = 55555. # arbitrary number
    P = opars['P'].to(u.day).value
    f0 = opars['phi0'].to(u.radian).value / (2*np.pi)
    _t = (np.array([0.02, 4.08, 4.45, 4.47]) + f0) * P

    ts = [
        np.concatenate((_t, (np.array([6.46, 7.04]) + f0) * P)) + EPOCH,
        np.concatenate((_t, (np.array([6.04, 6.07]) + f0) * P)) + EPOCH,
        np.concatenate((_t, (np.array([6.62, 6.65]) + f0) * P)) + EPOCH
    ]

    rv_err = np.random.uniform(0.2, 0.3, size=ts[0].size) * kms

    _rnd = np.random.normal(size=ts[0].size)

    rvs = []
    for t in ts:
        rvs.append(orbit.generate_rv_curve(t) + _rnd*rv_err)

    with h5py.File(os.path.join(cache_path, "experiment5.h5"), "w") as root:
        f = root.create_group('data')

        _data = RVData(t=ts[0], rv=rvs[0], stddev=rv_err)
        data0 = _data[:-2]

        g = f.create_group("0")
        data0.to_hdf5(g)
        # g.create_dataset('truth_vector', data=opars.pack())

        for i,t,rv in zip(range(len(ts)), ts, rvs):
            data = RVData(t=t, rv=rv, stddev=rv_err)
            g = f.create_group(str(i+1))
            data.to_hdf5(g)

        g = root.create_group('truth')
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
