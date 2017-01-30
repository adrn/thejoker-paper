#!/bin/bash
#SBATCH -J exp5           # job name
#SBATCH -o exp5.o%j             # output file name (%j expands to jobID)
#SBATCH -e exp5.e%j             # error file name (%j expands to jobID)
#SBATCH -n 260                   # total number of mpi tasks requested
#SBATCH -t 00:45:00             # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=adrn@princeton.edu
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job finishes

cd /tigress/adrianp/projects/thejoker-paper/scripts/

module load openmpi/gcc/1.10.2/64

source activate thejoker-paper

# Run experiment 5!
python make-experiment5-data.py -s 42

export NSAMPLES="2**28"
export SEED=42

# Run experiment 5, moving data point index=1 through 1 period in 8 steps:
srun python run-sampler.py -v --mpi -o \
-n $NSAMPLES \
-f ../cache/experiment5.h5 \
--data-key='0' --samples-key='0' \
--seed=$SEED

srun python run-sampler.py -v --mpi -o \
-n $NSAMPLES \
-f ../cache/experiment5.h5 \
--data-key='1' --samples-key='1' \
--seed=$SEED

srun python run-sampler.py -v --mpi -o \
-n $NSAMPLES \
-f ../cache/experiment5.h5 \
--data-key='2' --samples-key='2' \
--seed=$SEED

srun python run-sampler.py -v --mpi -o \
-n $NSAMPLES \
-f ../cache/experiment5.h5 \
--data-key='3' --samples-key='3' \
--seed=$SEED
