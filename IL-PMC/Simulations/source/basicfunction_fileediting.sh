#!/bin/bash

function edit_delete_last_lines()
{
	# Delete the last lines of a file until number of rows is equal for all columns
	lownum=$(awk '{print NF}' $1 | sort -nu | head -n 1) #gives the lowest column number
	highnum=$(awk '{print NF}' $1 | sort -nu | tail -n 1) #gives the highest column number
	while [ $lownum -lt $highnum ]; do
		head -n -1 $1 > tmp.dat
		mv tmp.dat $1
		lownum=$(awk '{print NF}' $1 | sort -nu | head -n 1)
	done
}
