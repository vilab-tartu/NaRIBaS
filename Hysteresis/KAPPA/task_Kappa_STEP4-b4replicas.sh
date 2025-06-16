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
r_ion_nm=$(echo ${current_ion[0]} | awk '{print $2}') #in nm
ion_r_name=$(echo ${current_ion[0]} | awk '{print $3}')
ion_charge_name=$(echo ${current_ion[0]} | awk '{print $4}')

#numberofions_name=$(echo ${current_numberofions[0]} | awk '{print $1}') 
temperature_name=$(echo ${current_temperature[0]} | awk '{print $1}')
version_name=$(echo ${current_version[0]} | awk '{print $1}')
replica_name=$(echo ${current_replica[0]} | awk '{print $1}')

####################################
# Calculations
xbox=$(awk "BEGIN {print "$xbox_nm*10.0"}" /dev/null)  #in A
ybox=$(awk "BEGIN {print "$ybox_nm*10.0"}" /dev/null)  #in A
zbox=$(awk "BEGIN {print "$zbox_nm*10.0"}" /dev/null)  #in A
xbox12=$(awk "BEGIN {print "$xbox*0.5"}" /dev/null)
ybox12=$(awk "BEGIN {print "$ybox*0.5"}" /dev/null)
zbox12=$(awk "BEGIN {print "$zbox*0.5"}" /dev/null)
r_wall=$(awk "BEGIN {print "$r_wall_nm*10.0"}" /dev/null)   #in A
z1_ion=$(awk "BEGIN {print "$zbox12+$r_wall_nm*10.0"}" /dev/null)   #in A
z2_ion=$(awk "BEGIN {print "$zbox12+$r_wall_nm*10.0+$r_ion_nm*10.0"}" /dev/null) #in A

A_box=$(awk "BEGIN {print "$xbox*$ybox"}" /dev/null)	#box area in sqA
#A_ion=$(awk "BEGIN {print "4*$r_ion_nm*$r_ion_nm*100"}" /dev/null)	#ion area in sqA
#maxnumberofions_name=$(awk "BEGIN {print int("$A_box/$A_ion")}" /dev/null)
#varnumberofions_name=$(awk "BEGIN {print "$maxnumberofions_name+1.0*$numberofions_name"}" /dev/null)

echo $varnumberofions_name



# Define path for storing simulation data
fullpath_start=$electrode_name/$ion_name/$temperature_name/$version_name/$replica_name/


# Take new n(mono), calculate the r_ion_nm and replace it in ion.list
#if [ -a $dir_experiments/$fullpath_start/nmono.txt ]
#then
#nmono=$(< $dir_experiments/$fullpath_start/nmono.txt)
#r2_ion=$(awk "BEGIN {print "$A_box/$nmono/400"}" /dev/null)
#r_ion=$(python -c "print round(float($r2_ion)**(0.5),3)")
#sed -i 's/'$ion_name'	'$r_ion_nm'/'$ion_name'	'$r_ion'/g' $inputlists_folder/ion.list
#fi

