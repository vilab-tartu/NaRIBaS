#!/bin/bash


### 2. Run singlepoint calculation of the model

# User input is specified below in the form of input lists stored in a specific folder.
# Don't change the variable names, only add/rename lists given within the brackets.

inputlists_folder=inputlists

inputlists=(ion.list slab.list method.list version.list)

####### Add here the file that contains the task definition
tasks=task_STEP2-singlepoint.sh

####### Add here the paths to scripts
currentdir=$(pwd)

dir_systempreparation=$currentdir/sysprep
dir_experiments=$currentdir/../Experiments_PMC/
