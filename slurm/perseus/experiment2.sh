#!/bin/bash
#SBATCH -J exp2           # job name
#SBATCH -o exp2.o%j             # output file name (%j expands to jobID)
#SBATCH -e exp2.e%j             # error file name (%j expands to jobID)
#SBATCH -n 260                   # total number of mpi tasks requested
#SBATCH -t 00:30:00             # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=adrn@princeton.edu
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job finishes

cd /tigress/adrianp/projects/thejoker-paper/scripts/

module load openmpi/gcc/1.10.2/64

source activate thejoker-paper

# Run experiment 2!
python make-experiment2-data.py -s 123

srun python run-sampler.py -v --mpi -o \
-n 2**28 -s 42 \
-f ../cache/experiment2.h5 \
--save-key='fixed-jitter-a' \
--read-key="a" \
--fixed-jitter='0 m/s'

srun python run-sampler.py -v --mpi -o \
-n 2**28 -s 42 \
-f ../cache/experiment2.h5 \
--save-key='fixed-jitter-b' \
--read-key="b" \
--fixed-jitter='0 m/s'

srun python run-sampler.py -v --mpi -o \
-n 2**28 -s 42 \
-f ../cache/experiment2.h5 \
--save-key='sample-jitter-b' \
--read-key='b' \
--log-jitter2-mean=10.5 --log-jitter2-std=0.5 --jitter-unit='m/s'
