# Standard library
import os
import sys

# Third-party
from astropy import log as logger
import astropy.units as u
import h5py
import numpy as np
from schwimmbad import choose_pool

# Project
from thejoker.sampler import JokerParams, TheJoker, save_prior_samples
from thejoker.data import RVData

if not os.path.exists(os.path.abspath("../scripts")):
    raise RuntimeError("Script must be run from within the scripts directory.")

def main(filename, pool, prior_samples_file, n_samples=1, seed=42, hdf5_key=None,
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
        joker_pars_kw['jitter'] = config['jitter']

    rerun = 0
    with h5py.File(filename, 'r+') as f:
        # get paths within HDF5 file
        if hdf5_key is not None:
            samples_path = 'samples/{}'.format(hdf5_key)
            data_path = 'data/{}'.format(hdf5_key)

        else:
            samples_path = 'samples'
            data_path = 'data'

        if samples_path in f:
            if not overwrite and not continue_sampling:
                logger.info("Sampling already performed for '{}':/data/{}. Use --overwrite "
                            "to redo or --continue to keep sampling.".format(filename, hdf5_key))
                pool.close()
                sys.exit(0)

            elif continue_sampling: # we need to increment the random number seed appropriately
                mode = 'a' # append to output file

                if 'rerun' not in f[samples_path].attrs:
                    rerun = 0
                else:
                    rerun = f[samples_path].attrs['rerun'] + 1

                if rerun != 0:
                    logger.debug("Reading hyperparameters from cache file.")

                    # HACK: this should be better
                    for name,unit in zip(['fixed_jitter', 'P_min', 'P_max'],
                                         [u.m/u.s, u.day, u.day]):
                        if np.isnan(f[samples_path].attrs[name]):
                            config[name] = None
                        else:
                            config[name] = f.attrs[name] * unit

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
        prior_units = save_prior_samples(prior_samples_file, prior_samples, u.km/u.s) # data is in km/s
        del prior_samples

        logger.debug("...done")

    else:
        with h5py.File(prior_samples_file, 'r') as f:
            prior_units = [u.Unit(uu) for uu in f.attrs['units']]

    logger.debug("Running sampler...")
    samples = joker.rejection_sample(data, prior_cache_file=prior_samples_file)

    return

    # save the orbital parameters out to a cache file
    with h5py.File(output_filename, mode) as f:
        f.attrs['rerun'] = rerun
        if hyperpars['fixed_jitter'] is not None:
            f.attrs['fixed_jitter'] = hyperpars['fixed_jitter']
        else:
            f.attrs['fixed_jitter'] = np.nan
        f.attrs['P_min'] = hyperpars['P_min']
        f.attrs['P_max'] = hyperpars['P_max']

        for i,(name,unit) in enumerate(OrbitalParams._name_to_unit.items()):
            if name in f:
                if overwrite: # delete old samples and overwrite
                    del f[name]
                    f.create_dataset(name, data=orbital_params.T[i])

                elif continue_sampling: # append to existing samples
                    _data = f[name][:]
                    del f[name]
                    f[name] = np.concatenate((_data, orbital_params.T[i]))
            else:
                f.create_dataset(name, data=orbital_params.T[i])

            f[name].attrs['unit'] = str(unit)

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
    parser.add_argument("--hdf5-key", dest="hdf5_key", default=None,
                        type=str, help="Path within an HDF5 file to the data.")
    parser.add_argument("-n", "--num-samples", dest="n_samples", default=2**20,
                        type=str, help="Number of prior samples to use in rejection sampling.")

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
            main(filename=args.filename, pool=pool, n_samples=n_samples, hdf5_key=args.hdf5_key,
                 prior_samples_file=fp.name, tmp_prior=True,
                 seed=args.seed, overwrite=args.overwrite, continue_sampling=args.continue_,
                 P_min=args.P_min, P_max=args.P_max, fixed_jitter=args.fixed_jitter,
                 log_jitter2_mean=args.log_jitter2_mean, log_jitter2_std=args.log_jitter2_std,
                 jitter_unit=args.jitter_unit)

    else:
        main(filename=args.filename, pool=pool, n_samples=n_samples, hdf5_key=args.hdf5_key,
             prior_samples_file=args.prior_filename,
             seed=args.seed, overwrite=args.overwrite, continue_sampling=args.continue_,
             P_min=args.P_min, P_max=args.P_max, fixed_jitter=args.fixed_jitter,
             log_jitter2_mean=args.log_jitter2_mean, log_jitter2_std=args.log_jitter2_std,
             jitter_unit=args.jitter_unit)
