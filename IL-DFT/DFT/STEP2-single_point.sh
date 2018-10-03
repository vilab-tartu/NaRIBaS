#!/bin/bash


### 2. Calculate single point energies with different functionals

# User input is specified below in the form of input lists stored in a specific folder.
# Don't change the variable names, only add/rename lists given within the brackets.

inputlists_folder=inputlists

# Note that this step uses a different list for functionals to iterate over multiple functionals
# while an optimization is only performed with one functional

inputlists=(anion.list cation.list functional_sp.list basisset.list version.list)

####### Add here the file that contains the task definition
tasks=task_STEP2-single_point.sh

####### Add here the paths to topologies, files and programs
currentdir=$(pwd)

dir_systempreparation=$currentdir/sysprep
dir_experiments=$currentdir/../Experiments_DFT/
dir_orca=/home/master/Downloads/Orca3/
