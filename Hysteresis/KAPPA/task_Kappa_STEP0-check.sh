#!/bin/bash

# Place to add code

# Some basic functions to access the data

#echo ${concentration[*]} # All elements of the array
#echo ${concentration[9]} # Element number 10 (counter starts at 0)
#numberofitems=${#concentration[*]} # calculate the number of elements in array
#echo $numberofitems
#echo ${concentration[1]} | awk '{print $1}' # Access the data that is stored within an element array
#echo '-----'
#print_current_setup

# Transfer list entries to bash variables
electrode_name=$(echo ${current_electrode[0]} | awk '{print $1}')
xbox_nm=$(echo ${current_electrode[0]} | awk '{print $2}') #in nm
ybox_nm=$(echo ${current_electrode[0]} | awk '{print $3}') #in nm
zbox_nm=$(echo ${current_electrode[0]} | awk '{print $4}') #in nm
r_wall_nm=$(echo ${current_electrode[0]} | awk '{print $5}') #in nm
electrodeatoms=$(echo ${current_electrode[0]} | awk '{print $6}')

ion_name=$(echo ${current_ion[0]} | awk '{print $1}')
dens_ion_nm2=$(echo ${current_ion[0]} | awk '{print $2}') #in nm
ion_r_name=$(echo ${current_ion[0]} | awk '{print $3}')
ion_charge_name=$(echo ${current_ion[0]} | awk '{print $4}')
ion_centre_name=$(echo ${current_ion[0]} | awk '{print $5}')

numberofions_name=$(echo ${current_numberofions[0]} | awk '{print $1}')

####################################
# Calculations
xbox=$(awk "BEGIN {print "$xbox_nm*10.0"}" /dev/null)  #in A
ybox=$(awk "BEGIN {print "$ybox_nm*10.0"}" /dev/null)  #in A

A_box=$(awk "BEGIN {print "$xbox*$ybox"}" /dev/null)    #box area in sqA
A_ion=$(awk "BEGIN {print "100/$dens_ion_nm2"}" /dev/null)      #ion area in sqA
maxnumberofions_name=$(awk "BEGIN {print int("$A_box/$A_ion")}" /dev/null)
varnumberofions_name=$(awk "BEGIN {print "$maxnumberofions_name+1.0*$numberofions_name"}" /dev/null)

echo $ion_name $maxnumberofions_name $varnumberofions_name
