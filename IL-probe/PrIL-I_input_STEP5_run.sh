#!/bin/bash


### 5. Perform production run


# User input is specified below in the form of input lists stored in a specific folder.
# Don't change the variable names, only add/rename lists given within the brackets.

inputlists_folder=inputlists

inputlists=(electrode.list anion.list cation.list additive.list probe.list orientation.list numberofionpairs.list mol_fraction.list temperature.list version.list surfacecharge.list distance.list replica.list)

####### Add here the file that contains the task definition
tasks=task_PrIL-I_STEP5-run.sh

####### Add here the paths to topologies, files and programs
currentdir=$(pwd)

dir_systempreparation=~/NaRiBaS/PrIL-I/Init
dir_experiments=~/NaRiBaS/PrIL-I/Experiments
dir_analysis=~/NaRiBaS/PrIL-I/Analysis
dir_temp=~/NaRiBaS/PrIL-I/Submit
dir_packmol=~/NaRiBaS/PrIL-I/Init
