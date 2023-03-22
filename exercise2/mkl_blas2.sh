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

export OMP_NUM_THREADS=64

for i in {1000..5000..500}
do
    for j in {1..5..1}
    do
        srun -n1 --cpus-per-task=12 ./gemm_mkl.x $i $i $i >> mkl.txt
        srun -n1 --cpus-per-task=12 ./gemm_oblas.x $i $i $i >> oblas.txt
    done
done

for i in {6000..10000..1000}
do
    for j in {1..5..1}
    do
        srun -n1 --cpus-per-task=12 ./gemm_mkl.x $i $i $i >> mkl.txt
        srun -n1 --cpus-per-task=12 ./gemm_oblas.x $i $i $i >> oblas.txt
    done
done

for i in {12000..20000..2000}
do
    for j in {1..5..1}
    do
        srun -n1 --cpus-per-task=12 ./gemm_mkl.x $i $i $i >> mkl.txt
        srun -n1 --cpus-per-task=12 ./gemm_oblas.x $i $i $i >> oblas.txt
    done
done
