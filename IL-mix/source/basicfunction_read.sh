#!/bin/bash

function read_sinlge_list()
{
   local  _resultvar=$1
   local  numofrows=$(cat $2 | wc -l)
   local  numofcolumns=$(awk '{print NF}' $2 | sort -nu | tail -n 1)
	local  row=$3
	local  column=0
	while [ $column -lt $numofcolumns ]; do
	  column=$((column + 1))
     local myresult[$column]=$(sed -n "$row p" $2 | cut -d' ' -f$column)
	done
   eval $_resultvar="'${myresult[*]}'"
}

function read_in_all_input_lists()
{
#   echo ${inputlists[*]}

#   input=${inputlists[*]}
#   path_to_lists=$Which_dir_to_read_in

   echo ''
   echo 'Starting the read in of the input data ...'
   echo ''
   echo 'Accessing lists in '$Which_dir_to_read_in
   echo 'Reading in the following lists:'
   echo ''

   listcount=0
   while [ $listcount -lt ${#inputlists[@]} ]; do
      variablename=${inputlists[$listcount]%.list}
      numrow=1
      while [ $numrow -le $(cat $Which_dir_to_read_in/${inputlists[$listcount]} | wc -l) ]; do
         read_sinlge_list $variablename[$((numrow-1))] $Which_dir_to_read_in/${inputlists[$listcount]} $numrow
         numrow=$((numrow+1))
      done

      eval numberofitems=\${#$variablename[*]}
      echo $variablename.list 'containing ' $numberofitems ' item(s).'

      listcount=$((listcount+1))
   done

   echo ''
   echo 'Finished data input.'
}
