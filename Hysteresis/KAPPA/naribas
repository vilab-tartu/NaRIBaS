#!/bin/bash

starttime=$(date +%s.%N)

currentdir=$(pwd)

####### User input will be read in here. The file should contain variable definition of 'inputlists_folder', 'inputlists', 'tasks'
source $currentdir/$1
#######

Which_dir_to_read_in=$currentdir/$inputlists_folder

source $currentdir/source/basicfunction_read.sh
source $currentdir/source/basicfunction_iterator.sh
source $currentdir/source/basicfunction_fileediting.sh
source $currentdir/source/basicfunction_slab.sh

# Reading in all data, done by a function stored in basicfunction_read.sh
read_in_all_input_lists
#read -p "Press enter to continue..." # allows to generate a break (e.g. to check the output of a former step) and proceed after enter is pressed).

# Now use the iterator pattern to access all elements of all lists, the basic functions are stored in the file basicfunction_iterator.sh

echo 'Execute task list.'

proceed=1
start=1

while [ $proceed -eq 1 ]; do
if [ $start -eq 1 ]; then
	start=0
# Initialize the array that holds for every list the index
	generate_initial_index_array
	generate_max_index_array
	totalnumberofsetups=0
	changed_list=1
	listcount=0
else
	if [ $changed_list -eq 1 ]; then
		totalnumberofsetups=$((totalnumberofsetups+1))
#		echo 'Next setup is'
#		print_current_setup $variablename $current_index_all
#		echo '------'

		get_current_setup

##############
		eval source $currentdir/$tasks
##############

		loop_through_all_elements
	else
		proceed=0
	fi
fi
done

echo ''
echo 'Total number of setups equals ' $totalnumberofsetups '.'

stoptime=$(date +%s.%N)
echo ''
printf "Script finished after:    %.3F seconds.\n"  $(echo "$stoptime - $starttime"|bc )

exit
