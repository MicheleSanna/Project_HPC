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

export OMP_PLACES=sockets
export OMP_PROC_BIND=master
export OMP_NUM_THREADS=12

for i in {1000..5000..500}
do
    for j in {1..15..1}
    do
        srun -n1 --cpus-per-task=12 ./gemm_mkl.x $i $i $i >> master_sockets/mkl_double.txt
        srun -n1 --cpus-per-task=12 ./gemm_oblas.x $i $i $i >> master_sockets/oblas_double.txt
    done
done

for i in {6000..10000..1000}
do
    for j in {1..15..1}
    do
        srun -n1 --cpus-per-task=12 ./gemm_mkl.x $i $i $i >> master_sockets/mkl_double.txt
        srun -n1 --cpus-per-task=12 ./gemm_oblas.x $i $i $i >> master_sockets/oblas_double.txt
    done
done

for i in {12000..20000..2000}
do
    for j in {1..15..1}
    do
        srun -n1 --cpus-per-task=12 ./gemm_mkl.x $i $i $i >> master_sockets/mkl_double.txt
        srun -n1 --cpus-per-task=12 ./gemm_oblas.x $i $i $i >> master_sockets/oblas_double.txt
    done
done
