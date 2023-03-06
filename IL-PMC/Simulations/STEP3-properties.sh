#!/bin/bash

### 3. Calculate properties

# User input is specified below in the form of input lists stored in a specific folder.
# Don't change the variable names, only add/rename lists given within the brackets.

inputlists_folder=inputlists

inputlists=(ion.list slab.list method.list version.list)

####### Add here the file that contains the task definition
tasks=task_STEP3-properties.sh

####### Add here the paths to scripts
currentdir=$(pwd)

dir_systempreparation=$currentdir/sysprep
dir_experiments=$currentdir/../Experiments_PMC/
