### information about exercise 2
All the runs were done on the epyc node\\

#Explanation of the files

-*_double.txt contains all the results of of a script where, with 64 OMP threads, the benchmark were runned with different sizes of matrices from 1000 to 20000, in double precision

-*_experiments.txt contains results of runs with different combination of OMP_PLACES and OMP_PROC_BIND, in this order: (cores,master), (cores,spread), (cores,close), (threads,master), (threads,spread), (threads,close)

-oblas.txt and mkl.txt contains all the results of the runs for testing the scaling at the increasing of openMP threads

-mkl_blas.sh and mkl_blas2.sh are scripts used to run the tests, mkl_blas2 for the scaling at the increasing of the matrix dimension and mkl_blas for the scaling at the increasing of the OpenMP threads.

The files containing the results for the matrices scaling test but with float were unfortunately overwritten, but the results of the tests are in the report and in this folder there is the script that generated them (mkl_blas2.sh).

