#!/bin/bash


### 2. Perform production run


# User input is specified below in the form of input lists stored in a specific folder.
# Don't change the variable names, only add/rename lists given within the brackets.

inputlists_folder=inputlists_Kappa

inputlists=(electrode.list ion.list temperature.list version.list replica.list)

####### Add here the file that contains the task definition
tasks=task_Kappa_STEP5-copy.sh

####### Add here the paths to topologies, files and programs
currentdir=$(pwd)

dir_systempreparation=$currentdir/sysprep_Kappa
dir_experiments=$currentdir/../NaRIBaS_Data/Experiments_Kappa
dir_analysis=$currentdir/../NaRIBaS_Data/Analysis_Kappa
dir_temp=$currentdir/../NaRIBaS_Data/Submit_Kappa
dir_packmol=$currentdir/sysprep_Kappa
