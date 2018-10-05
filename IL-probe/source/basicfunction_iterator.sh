#!/bin/bash

function generate_initial_index_array()
{
	count=0
	while [ $count -lt ${#inputlists[@]} ]; do
		current_index_all[$count]=0
		count=$((count+1))
	done

#	echo ${current_index_all[*]}
}


function generate_max_index_array()
{
	count=0
	while [ $count -lt ${#inputlists[@]} ]; do
		variablename=${inputlists[$count]%.list}
		eval numberofitems=\${#$variablename[*]}
		max_index_all[$count]=$numberofitems
		count=$((count+1))
	done

#	echo ${max_index_all[*]}
}

function print_current_setup()
{
   local  variable=$1
   local  listindex_all=$2

	count=0

	while [ $count -lt ${#inputlists[@]} ]; do
		variablename=${inputlists[$count]%.list}
		i=${current_index_all[$count]}

#		test=$(eval echo \${$variablename[$i]})
#		echo $variablename $i $(eval echo \${$variablename[$i]}) #${!variablename[$i]}
		eval echo \${$variablename[$i]}

		count=$((count+1))
	done
}

function get_current_setup()
{
	count=0

	while [ $count -lt ${#inputlists[@]} ]; do
		variablename=${inputlists[$count]%.list}
		i=${current_index_all[$count]}

#		test=$(eval echo \${$variablename[$i]})
#		echo $variablename $i $(eval echo \${$variablename[$i]}) #${!variablename[$i]}
		eval current_$variablename=\${$variablename[$i]}
#		eval echo \$current_$variablename
		count=$((count+1))
	done
#   echo $variablename $i $(eval echo \${$variablename[$i]}) #${!variablename[$i]}
}

function hasnext()
{
   local  variable=$1
   local  listindex=$2

   eval numberofitems=\${#$variable[*]}
#   echo $numberofitems

   if [ $((listindex+1)) -lt $numberofitems ]; then
      echo 1
   else
      echo 0
   fi
}

function next()
{
   local  variable=$1
   local  listindex=$2

   echo $((listindex+1))
}


function loop_through_all_elements()
{
listcount=0
while [ $listcount -lt ${#inputlists[@]} ]; do
   variablename=${inputlists[$listcount]%.list}
   current_index=${current_index_all[$listcount]}

   if [ $(hasnext $variablename $current_index) -eq 1 ]; then
		current_index_all[$listcount]=$(next $variablename $current_index)
		current_index=${current_index_all[$listcount]} 
		changed_list=1
		break
   elif [ $(hasnext $variablename $current_index) -eq 0 ]; then
		changed_list=0
		current_index_all[$listcount]=0
#		echo ${current_index_all[$listcount]}
#		listcount=$((listcount-1))
   fi
   listcount=$((listcount+1))
done
}

