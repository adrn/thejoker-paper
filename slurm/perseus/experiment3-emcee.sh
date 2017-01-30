#!/bin/bash
#SBATCH -J exp3-emcee           # job name
#SBATCH -o exp3-emcee.o%j             # output file name (%j expands to jobID)
#SBATCH -e exp3-emcee.e%j             # error file name (%j expands to jobID)
#SBATCH -n 130                   # total number of mpi tasks requested
#SBATCH -t 02:00:00             # run time (hh:mm:ss)
#SBATCH --mail-user=adrn@princeton.edu
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job finishes

cd /tigress/adrianp/projects/thejoker-paper/scripts/

module load openmpi/gcc/1.10.2/64

source activate thejoker-paper

export NSTEPS=16384
export SEED=42

# Run emcee on output from experiment 3
srun python continue-with-emcee.py -v --mpi -o \
--nsteps=$NSTEPS \
-f ../cache/experiment3.h5 \
--data-key="11" \
--samples-key="11" \
--seed=$SEED

srun python continue-with-emcee.py -v --mpi -o \
--nsteps=$NSTEPS \
-f ../cache/experiment3.h5 \
--data-key="9" \
--samples-key="9" \
--seed=$SEED

# Also continue this one:
srun python run-sampler.py -v --mpi -c \
-n 2**28 \
-f ../cache/experiment3.h5 \
--data-key="7" \
--samples-key="7" \
--seed=42 \
--fixed-jitter='0 m/s'
