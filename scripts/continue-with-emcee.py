"""

TODO: this script doesn't support inferring the jitter and can only handle setting jitter = 0

"""


# Standard library
import os
import sys

# Third-party
import astropy.units as u
import emcee
import h5py
import numpy as np
from schwimmbad import choose_pool

# Project
from thejoker.log import log as logger
from thejoker.sampler import JokerParams
from thejoker.data import RVData
from thejoker.utils import quantity_to_hdf5, quantity_from_hdf5
from thejoker.sampler import mcmc

if not os.path.exists(os.path.abspath("../scripts")):
    raise RuntimeError("Script must be run from within the scripts directory.")

# TODO: single filename, blah blah...
def main(filename, pool, n_steps, overwrite=False, seed=42,
         data_key=None, samples_key=None):

    filename = os.path.abspath(filename)

    # make sure The Joker has already been run
    if not os.path.exists(filename):
        raise IOError("The Joker cache file '{}' can't be found! Are you sure you ran "
                      "'run-sampler.py'?".format(filename))

    with h5py.File(filename, 'a') as root:
        if data_key is not None:
            data_path = 'data/{}'.format(data_key)
        else:
            data_path = 'data'

        if samples_key is not None:
            samples_path = 'samples/{}'.format(samples_key)
        else:
            samples_path = 'samples'

        emcee_path = 'emcee/{}'.format(samples_key)

        if samples_path not in root:
            raise IOError("Key '{}' not in HDF5 file '{}'".format(samples_path, filename))

        # TODO:
        # jitter = f.attrs['fixed_jitter'] * u.m/u.s

        if emcee_path in root:
            if overwrite:
                del root[emcee_path]

            else:
                logger.info("emcee already run on this sample file.")
                pool.close()
                sys.exit(0)

        elif 'emcee' not in root:
            root.create_group('emcee')

    # do this after choosing pool so all processes don't have same seed
    if seed is not None:
        logger.debug("random number seed: {}".format(seed))
        np.random.seed(seed)

    # only accepts HDF5 data formats with units
    logger.debug("Reading data and cached samples from file at '{}'".format(filename))
    with h5py.File(filename, 'r') as f:
        data = RVData.from_hdf5(f[data_path])

        # also read samples
        samples = dict()
        for key in f[samples_path].keys():
            samples[key] = quantity_from_hdf5(f[samples_path], key)

    n_joker_samples = len(samples['P'])

    # TODO: need to load this from HDF5 file...
    # HACK: fixed jitter
    params = JokerParams(P_min=8*u.day, P_max=8192*u.day, anomaly_tol=1E-11,
                         jitter=0.*u.m/u.s)

    # Fire up emcee
    if n_joker_samples > 1:
        P = np.median(samples['P'])
        T = data._t_bmjd.max() - data._t_bmjd.min()
        Delta = 4*P**2 / (2*np.pi*T*u.day)
        P_rms = np.std(samples['P'])
        logger.debug("Period rms for surviving samples: {}, Delta: {}".format(P_rms, Delta))

        if P_rms > Delta:
            logger.error("Period rms > ∆! Re-run The Joker instead - your posterior pdf is "
                         "probably multi-modal")
            pool.close()
            sys.exit(0)

    else:
        logger.debug("Only one surviving sample.")

    # transform samples to the parameters we'll sample using emcee
    _names = 'P', 'phi0', 'ecc', 'omega', 'jitter', 'K', 'v0'
    samples_vec = np.array([samples[k].value for k in _names]).T
    samples_trans = mcmc.to_mcmc_params(samples_vec.T).T

    # TODO HACK: to remove jitter, because it's fixed
    samples_trans = np.delete(samples_trans, 5, axis=1)

    j_max = np.argmax([mcmc.ln_posterior(s, params, data) for s in samples_trans])
    p0 = samples_trans[j_max]

    # HACK: config
    M_min = 128
    n_walkers = M_min
    p0 = emcee.utils.sample_ball(p0, 1E-5*np.abs(p0), size=n_walkers)

    sampler = emcee.EnsembleSampler(n_walkers, p0.shape[1],
                                    lnpostfn=mcmc.ln_posterior, args=(params,data),
                                    pool=pool)

    pos,prob,state = sampler.run_mcmc(p0, n_steps) # MAGIC NUMBER

    pool.close()

    pos = np.hstack((pos, np.zeros((pos.shape[0],1))))
    emcee_samples_vec = mcmc.from_mcmc_params(pos.T)

    import matplotlib.pyplot as plt
    # import corner
    # corner.corner(np.vstack(sampler.chain[:,:-128]))

    nwalkers, nlinks, dim = sampler.chain.shape
    for k in range(dim):
        plt.figure()
        for n in range(nwalkers):
            plt.plot(sampler.chain[n,:,k], marker='', drawstyle='steps', alpha=0.1)

    plt.show()

    # TODO: turn vec into dict with units

    with h5py.File(filename, 'a') as root:
        f = root[emcee_path]
        for key,val in emcee_samples.items():
            quantity_to_hdf5(f, key, val)

if __name__ == "__main__":
    from argparse import ArgumentParser
    import logging

    # Define parser object
    parser = ArgumentParser(description="")

    vq_group = parser.add_mutually_exclusive_group()
    vq_group.add_argument('-v', '--verbose', action='count', default=0, dest='verbosity')
    vq_group.add_argument('-q', '--quiet', action='count', default=0, dest='quietness')

    oc_group = parser.add_mutually_exclusive_group()
    oc_group.add_argument("-o", "--overwrite", dest="overwrite", default=False,
                          action="store_true", help="Overwrite any existing data.")
    oc_group.add_argument("-c", "--continue", dest="_continue", default=False,
                          action="store_true", help="Continue the sampler.")

    parser.add_argument("--seed", dest="seed", default=None, type=int,
                        help="Random number seed")

    # emcee
    parser.add_argument("--nsteps", dest="n_steps", required=True, type=int,
                        help="Number of steps to take.")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--procs", dest="n_procs", default=1,
                       type=int, help="Number of processes.")
    group.add_argument("--mpi", dest="mpi", default=False,
                       action="store_true", help="Run with MPI.")

    parser.add_argument("-f", "--file", dest="filename", default=None, required=True,
                        type=str, help="Path to HDF5 data file to analyze.")
    parser.add_argument("--data-key", dest="data_key", default=None,
                        type=str, help="Path within HDF5 file to read the data.")
    parser.add_argument("--samples-key", dest="samples_key", default=None,
                        type=str, help="Path within HDF5 file to save the samples.")

    # TODO: add other jitter options from run-sampler.py
    # THIS IS IGNORED!
    parser.add_argument("--fixed-jitter", dest="fixed_jitter", default=None, type=str,
                        help="Extra uncertainty to add in quadtrature to the RV measurement "
                             "uncertainties. Must specify a number with units, e.g., '15 m/s'")

    args = parser.parse_args()

    if args.fixed_jitter is not None:
        raise NotImplementedError()

    # Set logger level based on verbose flags
    if args.verbosity != 0:
        if args.verbosity == 1:
            logger.setLevel(logging.DEBUG)
        else: # anything >= 2
            logger.setLevel(1)

    elif args.quietness != 0:
        if args.quietness == 1:
            logger.setLevel(logging.WARNING)
        else: # anything >= 2
            logger.setLevel(logging.ERROR)

    else: # default
        logger.setLevel(logging.INFO)

    pool = choose_pool(mpi=args.mpi, processes=args.n_procs)

    # use a context manager so the prior samples file always gets deleted
    main(filename=args.filename, pool=pool, n_steps=args.n_steps,
         seed=args.seed, overwrite=args.overwrite,
         data_key=args.data_key, samples_key=args.samples_key)
