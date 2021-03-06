# Standard library
import os
import sys

# Third-party
import astropy.units as u
import h5py
import numpy as np
from schwimmbad import choose_pool

# Project
from thejoker.log import log as logger
from thejoker.sampler import JokerParams, TheJoker, save_prior_samples
from thejoker.data import RVData
from thejoker.utils import quantity_to_hdf5, quantity_from_hdf5

if not os.path.exists(os.path.abspath("../scripts")):
    raise RuntimeError("Script must be run from within the scripts directory.")

def main(filename, pool, prior_samples_file, n_samples=1, seed=42,
         data_key=None, samples_key=None,
         overwrite=False, continue_sampling=False, tmp_prior=False,
         P_min=None, P_max=None, fixed_jitter=None,
         log_jitter2_mean=None, log_jitter2_std=None, jitter_unit=None):

    filename = os.path.abspath(filename)
    prior_samples_file = os.path.abspath(prior_samples_file)

    if not os.path.exists(filename):
        raise IOError("Data / cache file '{}' does not exist!".format(filename))

    # Configuration / hyper-parameters
    config = dict()

    if P_min is None:
        config['P_min'] = 8. * u.day
    else:
        config['P_min'] = float(P_min) * u.day

    if P_max is None:
        config['P_max'] = 8192. * u.day
    else:
        config['P_max'] = float(P_max) * u.day

    joker_pars_kw = dict()
    if fixed_jitter is None and log_jitter2_mean is None:
        joker_pars_kw['jitter'] = 0. * u.m/u.s

    elif fixed_jitter is None and log_jitter2_mean is not None:
        assert log_jitter2_std is not None
        assert jitter_unit is not None

        joker_pars_kw['jitter'] = (log_jitter2_mean, log_jitter2_std)
        joker_pars_kw['jitter_unit'] = u.Unit(str(jitter_unit))

    else:
        val,unit = fixed_jitter.split()
        joker_pars_kw['jitter'] = float(val) * u.Unit(unit)

    # paths within HDF5 file
    if data_key is not None:
        data_path = 'data/{}'.format(data_key)
    else:
        data_path = 'data'

    if samples_key is not None:
        samples_path = 'samples/{}'.format(samples_key)
    else:
        samples_path = 'samples'

    rerun = 0
    with h5py.File(filename, 'r+') as f:
        if samples_path in f:
            if not overwrite and not continue_sampling:
                logger.info("Sampling already performed for {}:{} (saved in {}). Use --overwrite "
                            "to redo or --continue to keep sampling.".format(filename,
                                                                             data_path,
                                                                             samples_path))
                pool.close()
                sys.exit(0)

            elif continue_sampling: # we need to increment the random number seed appropriately
                if 'rerun' not in f[samples_path].attrs:
                    rerun = 0
                else:
                    rerun = f[samples_path].attrs['rerun'] + 1

            elif overwrite: # restart rerun counter
                rerun = 0
                del f[samples_path]

            else:
                raise ValueError("Unknown state!")

    # do this after choosing pool so all processes don't have same seed
    if seed is not None:
        if rerun > 0:
            logger.info("This is rerun {} -- incrementing random number seed.".format(rerun))

        logger.debug("random number seed: {} (+ rerun: {})".format(seed, rerun))
        seed = int(str(seed) + str(rerun)) # why do i do this dumb thing?
        rnd = np.random.RandomState(seed=seed)

    else:
        rnd = np.random.RandomState()

    # create TheJoker sampler instance
    params = JokerParams(P_min=config['P_min'], P_max=config['P_max'], anomaly_tol=1E-11,
                         **joker_pars_kw)
    joker = TheJoker(params, random_state=rnd, pool=pool)

    logger.debug("Reading data from input file at '{}'".format(filename))
    with h5py.File(filename, 'r') as f:
        data = RVData.from_hdf5(f[data_path])

    # create prior samples cache, store to file and store filename in DB
    logger.debug("Number of prior samples: {}".format(n_samples))
    if not os.path.exists(prior_samples_file) or tmp_prior:
        logger.debug("Generating prior samples now, caching to '{}'...".format(prior_samples_file))

        prior_samples = joker.sample_prior(n_samples)
        save_prior_samples(prior_samples_file, prior_samples, u.km/u.s) # data is in km/s
        del prior_samples

        logger.debug("...done")

    logger.debug("Running sampler...")
    samples = joker.rejection_sample(data, prior_cache_file=prior_samples_file)

    # save the orbital parameters out to a cache file
    with h5py.File(filename, 'a') as root:
        if samples_path not in root:
            root.create_group(samples_path)

        f = root[samples_path]
        f.attrs['rerun'] = rerun

        for key,val in samples.items():
            if key in f:
                if overwrite: # delete old samples and overwrite
                    del f[key]
                    quantity_to_hdf5(f, key, val)

                elif continue_sampling: # append to existing samples
                    _data = quantity_from_hdf5(f, key)
                    del f[key]
                    _data = np.concatenate((_data.to(val.unit).value, val.value))
                    quantity_to_hdf5(f, key, _data * val.unit)

            else:
                quantity_to_hdf5(f, key, val)

    pool.close()

if __name__ == "__main__":
    from argparse import ArgumentParser
    import logging
    import tempfile

    # Define parser object
    parser = ArgumentParser(description="")

    vq_group = parser.add_mutually_exclusive_group()
    vq_group.add_argument('-v', '--verbose', action='count', default=0, dest='verbosity')
    vq_group.add_argument('-q', '--quiet', action='count', default=0, dest='quietness')

    oc_group = parser.add_mutually_exclusive_group()
    oc_group.add_argument("-o", "--overwrite", dest="overwrite", default=False,
                          action="store_true", help="Overwrite any existing data.")
    oc_group.add_argument("-c", "--continue", dest="continue_", default=False,
                          action="store_true", help="Continue the sampler.")

    parser.add_argument("-s", "--seed", dest="seed", default=None, type=int,
                        help="Random number seed")

    group = parser.add_mutually_exclusive_group()
    group.add_argument("--procs", dest="n_procs", default=1,
                       type=int, help="Number of processes.")
    group.add_argument("--mpi", dest="mpi", default=False,
                       action="store_true", help="Run with MPI.")

    parser.add_argument("-f", "--filename", dest="filename", default=None, required=True,
                        type=str, help="Path to HDF5 data file to analyze.")
    parser.add_argument("--prior-filename", dest="prior_filename", default=None,
                        type=str, help="Path to (store or load) prior samples.")
    parser.add_argument("-n", "--num-samples", dest="n_samples", default=2**20,
                        type=str, help="Number of prior samples to use in rejection sampling.")

    parser.add_argument("--data-key", dest="data_key", default=None,
                        type=str, help="Path within HDF5 file to read the data.")
    parser.add_argument("--samples-key", dest="samples_key", default=None,
                        type=str, help="Path within HDF5 file to save the samples.")

    parser.add_argument("--fixed-jitter", dest="fixed_jitter", default=None, type=str,
                        help="Extra uncertainty to add in quadtrature to the RV measurement "
                             "uncertainties. Must specify a number with units, e.g., '15 m/s'")
    parser.add_argument("--Pmin", dest="P_min", default=8., type=float,
                        help="Minimum period for period prior in days.")
    parser.add_argument("--Pmax", dest="P_max", default=8192, type=float,
                        help="Maximum period for perior prior in days.")

    parser.add_argument("--log-jitter2-mean", dest="log_jitter2_mean", default=None, type=float,
                        help="TODO")
    parser.add_argument("--log-jitter2-std", dest="log_jitter2_std", default=None, type=float,
                        help="TODO")
    parser.add_argument("--jitter-unit", dest="jitter_unit", default=None, type=str,
                        help="TODO")

    args = parser.parse_args()

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

    try:
        n_samples = int(args.n_samples)
    except:
        n_samples = int(eval(args.n_samples)) # LOL what's security?

    pool_kwargs = dict(mpi=args.mpi, processes=args.n_procs)
    pool = choose_pool(**pool_kwargs)

    if args.prior_filename is None:
        with tempfile.NamedTemporaryFile(dir=os.path.abspath("../cache")) as fp:
            main(filename=args.filename, pool=pool, n_samples=n_samples,
                 data_key=args.data_key, samples_key=args.samples_key,
                 prior_samples_file=fp.name, tmp_prior=True,
                 seed=args.seed, overwrite=args.overwrite, continue_sampling=args.continue_,
                 P_min=args.P_min, P_max=args.P_max, fixed_jitter=args.fixed_jitter,
                 log_jitter2_mean=args.log_jitter2_mean, log_jitter2_std=args.log_jitter2_std,
                 jitter_unit=args.jitter_unit)

    else:
        main(filename=args.filename, pool=pool, n_samples=n_samples,
             data_key=args.data_key, samples_key=args.samples_key,
             prior_samples_file=args.prior_filename,
             seed=args.seed, overwrite=args.overwrite, continue_sampling=args.continue_,
             P_min=args.P_min, P_max=args.P_max, fixed_jitter=args.fixed_jitter,
             log_jitter2_mean=args.log_jitter2_mean, log_jitter2_std=args.log_jitter2_std,
             jitter_unit=args.jitter_unit)
