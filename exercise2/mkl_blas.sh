#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --job-name=mklblasrun
#SBATCH --nodes=1
#SBATCH --cpus-per-task 64
#SBATCH --time=01:59:00
#SBATCH --nodelist=epyc[008]

module load architecture/AMD
module load mkl
module load openBLAS/0.3.21-omp
srun -n1 make cpu

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
srun -n1 --cpus-per-task=12 ./gemm_mkl.x 10000 10000 10000 >> mkl.txt
srun -n1 --cpus-per-task=12 ./gemm_oblas.x 10000 10000 10000 >> oblas.txt
export OMP_NUM_THREADS=2
srun -n1 --cpus-per-task=12 ./gemm_mkl.x 10000 10000 10000 >> mkl.txt
srun -n1 --cpus-per-task=12 ./gemm_oblas.x 10000 10000 10000 >> oblas.txt
export OMP_NUM_THREADS=4
srun -n1 --cpus-per-task=12 ./gemm_mkl.x 10000 10000 10000 >> mkl.txt
srun -n1 --cpus-per-task=12 ./gemm_oblas.x 10000 10000 10000 >> oblas.txt
export OMP_NUM_THREADS=8
srun -n1 --cpus-per-task=12 ./gemm_mkl.x 10000 10000 10000 >> mkl.txt
srun -n1 --cpus-per-task=12 ./gemm_oblas.x 10000 10000 10000 >> oblas.txt
export OMP_NUM_THREADS=16
srun -n1 --cpus-per-task=12 ./gemm_mkl.x 10000 10000 10000 >> mkl.txt
srun -n1 --cpus-per-task=12 ./gemm_oblas.x 10000 10000 10000 >> oblas.txt
export OMP_NUM_THREADS=32
srun -n1 --cpus-per-task=12 ./gemm_mkl.x 10000 10000 10000 >> mkl.txt
srun -n1 --cpus-per-task=12 ./gemm_oblas.x 10000 10000 10000 >> oblas.txt
export OMP_NUM_THREADS=64
srun -n1 --cpus-per-task=12 ./gemm_mkl.x 10000 10000 10000 >> mkl.txt
srun -n1 --cpus-per-task=12 ./gemm_oblas.x 10000 10000 10000 >> oblas.txt
export OMP_NUM_THREADS=128
srun -n1 --cpus-per-task=12 ./gemm_mkl.x 10000 10000 10000 >> mkl.txt
srun -n1 --cpus-per-task=12 ./gemm_oblas.x 10000 10000 10000 >> oblas.txt
