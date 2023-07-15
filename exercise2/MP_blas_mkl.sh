#!/bin/bash
#SBATCH --partition=THIN
#SBATCH --job-name=mklblasrun
#SBATCH --nodes=1
#SBATCH --cpus-per-task 12
#SBATCH --time=01:59:00
#SBATCH --exclusive
#SBATCH --nodelist=thin[008]

module load mkl
module load openBLAS/0.3.23-omp
srun -n1 make cpu

for j in {1..9..1}
do
    srun -n1 --cpus-per-task=1 ./gemm_mkl.x 8000 8000 8000 >> mpi/mkl_float.txt
    srun -n1 --cpus-per-task=1 ./gemm_oblas.x 8000 8000 8000 >> mpi/oblas_float.txt
done

for j in {1..9..1}
do
    srun -n1 --cpus-per-task=2 ./gemm_mkl.x 8000 8000 8000 >> mpi/mkl_float.txt
    srun -n1 --cpus-per-task=2 ./gemm_oblas.x 8000 8000 8000 >> mpi/oblas_float.txt
done

for j in {1..9..1}
do
    srun -n1 --cpus-per-task=4 ./gemm_mkl.x 8000 8000 8000 >> mpi/mkl_float.txt
    srun -n1 --cpus-per-task=4 ./gemm_oblas.x 8000 8000 8000 >> mpi/oblas_float.txt
done

for j in {1..9..1}
do
    srun -n1 --cpus-per-task=6 ./gemm_mkl.x 8000 8000 8000 >> mpi/mkl_float.txt
    srun -n1 --cpus-per-task=6 ./gemm_oblas.x 8000 8000 8000 >> mpi/oblas_float.txt
done

for j in {1..9..1}
do
    srun -n1 --cpus-per-task=8 ./gemm_mkl.x 8000 8000 8000 >> mpi/mkl_float.txt
    srun -n1 --cpus-per-task=8 ./gemm_oblas.x 8000 8000 8000 >> mpi/oblas_float.txt
done

for j in {1..9..1}
do
    srun -n1 --cpus-per-task=10 ./gemm_mkl.x 8000 8000 8000 >> mpi/mkl_float.txt
    srun -n1 --cpus-per-task=10 ./gemm_oblas.x 8000 8000 8000 >> mpi/oblas_float.txt
done

for j in {1..9..1}
do
    srun -n1 --cpus-per-task=12 ./gemm_mkl.x 8000 8000 8000 >> mpi/mkl_float.txt
    srun -n1 --cpus-per-task=12 ./gemm_oblas.x 8000 8000 8000 >> mpi/oblas_float.txt
done

