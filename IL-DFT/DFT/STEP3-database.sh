#!/bin/bash


### 3. Calculate SP properties

# User input is specified below in the form of input lists stored in a specific folder.
# Don't change the variable names, only add/rename lists given within the brackets.

inputlists_folder=inputlists

inputlists=(anion.list cation.list functional_sp.list basisset.list version.list)

####### Add here the file that contains the task definition
tasks=task_STEP3-database.sh

####### Add here the paths to topologies, files and programs
currentdir=$(pwd)

dir_systempreparation=$currentdir/sysprep
dir_experiments=$currentdir/../Experiments_DFT/
dir_orca=/home/master/Downloads/Orca3/
