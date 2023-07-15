
# Exercise 1 - Parallel Programming

### Foundations of HPC @ UniTS-DSSC 2022

This project should run just by loading the openMPI module and running the make command. Currently it's compiled on epyc nodes. This project follows the rules specified in the assignment pdf, the only different thing is that the initialized playground is not void, it contains a bunch of gliders to check if the program is running correctly

The rules of the project were different with respect to the original one, so i made a macro to switch between them. Right now the project is compiled with the official version of the rules

The macro IMAGE_TEST is useful to create game fields that make easier the debug