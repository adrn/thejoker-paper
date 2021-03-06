#!/bin/bash
#SBATCH -J exp3           # job name
#SBATCH -o exp3.o%j             # output file name (%j expands to jobID)
#SBATCH -e exp3.e%j             # error file name (%j expands to jobID)
#SBATCH -n 260                   # total number of mpi tasks requested
#SBATCH -t 01:00:00             # run time (hh:mm:ss)
#SBATCH --mail-user=adrn@princeton.edu
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job finishes

cd /tigress/adrianp/projects/thejoker-paper/scripts/

module load openmpi/gcc/1.10.2/64

source activate thejoker-paper

python make-experiment3-data.py -s 42

export NSAMPLES="2**28"
export SEED=42

srun python run-sampler.py -v --mpi -o \
-n $NSAMPLES \
-f ../cache/experiment3.h5 \
--data-key="11" \
--samples-key="11" \
--seed=$SEED \
--fixed-jitter='0 m/s'

srun python run-sampler.py -v --mpi -o \
-n $NSAMPLES \
-f ../cache/experiment3.h5 \
--data-key="9" \
--samples-key="9" \
--seed=$SEED \
--fixed-jitter='0 m/s'

srun python run-sampler.py -v --mpi -o \
-n $NSAMPLES \
-f ../cache/experiment3.h5 \
--data-key="7" \
--samples-key="7" \
--seed=$SEED \
--fixed-jitter='0 m/s'

srun python run-sampler.py -v --mpi -o \
-n $NSAMPLES \
-f ../cache/experiment3.h5 \
--data-key="5" \
--samples-key="5" \
--seed=$SEED \
--fixed-jitter='0 m/s'

srun python run-sampler.py -v --mpi -o \
-n $NSAMPLES \
-f ../cache/experiment3.h5 \
--data-key="3" \
--samples-key="3" \
--seed=$SEED \
--fixed-jitter='0 m/s'
