#!/bin/bash
PBS -N test
PBS -q small
PBS -l nodes=3:ppn=4
PBS -l walltime=00:05:00
cd $PBS_O_WORKDIR
source /opt/software/intel/initpaths intel64
export I_MPI_FABRICS=shm:dapl
mpirun -np 12 ./exe 4 16777216 10 101