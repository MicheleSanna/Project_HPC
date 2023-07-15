#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=mklblasrun
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 1
#SBATCH --time=01:59:00
#SBATCH --nodelist=epyc[001]
#SBATCH --exclusive


module load architecture/AMD
module load openMPI/4.1.4/gnu/12.2.1

srun make 

export OMP_NUM_THREADS=1

srun project -i -k 20000 -e 1 -f snapshot_00000.pgm -n 100 -s 0
