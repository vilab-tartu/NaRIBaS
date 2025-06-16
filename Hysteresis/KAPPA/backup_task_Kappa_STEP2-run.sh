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
zbox12=$(awk "BEGIN {print "$zbox*0.1"}" /dev/null)
r_wall=$(awk "BEGIN {print "$r_wall_nm*10.0"}" /dev/null)   #in A
z1_ion=$(awk "BEGIN {print "$zbox12+$r_wall"}" /dev/null)   #in A
z2_ion=$(awk "BEGIN {print "$zbox12+5*2+$r_wall"}" /dev/null) #in A

A_box=$(awk "BEGIN {print "$xbox*$ybox"}" /dev/null)    #box area in sqA
A_ion=$(awk "BEGIN {print "100/$dens_ion_nm2"}" /dev/null)      #ion area in sqA
maxnumberofions_name=$(awk "BEGIN {print int("$A_box/$A_ion")}" /dev/null)
varnumberofions_name=$(awk "BEGIN {print "$maxnumberofions_name+1.0*$numberofions_name"}" /dev/null)
electrode_atom_charge=$(awk "BEGIN {printf \"%.20f\n\","-1.0*$ion_charge_name*$varnumberofions_name/$electrodeatoms"}" /dev/null)

####################################

# Define path for storing simulation data
fullpath_start=$electrode_name/$ion_name/$temperature_name/$version_name/$replica_name/$varnumberofions_name

cd $dir_systempreparation

#Add all necessary .mdp files
sed 's/SED_temperature_name_SED/'$temperature_name'/g' 2_NVT_production_varTemp.mdp > 2_NVT_run.mdp
mv 2_NVT_run.mdp $dir_experiments/$fullpath_start/

####################################
cd $dir_experiments/$fullpath_start/

#prepare the index file
gmx make_ndx -f NVT_lowtimestep.gro -o index.ndx << EOF
keep 0
r $electrode_name
name 1 Electrode
r $ion_r_name
name 2 Ion
a $ion_centre_name
name 3 Centre
q
EOF

if [ -f NVT.gro ]; then
    echo "Calculation already done"
else
    #grompp and initiate production run
    gmx grompp -f 2_NVT_run.mdp -c NVT_lowtimestep.gro -p topol_local.top -n index.ndx -o NVT -maxwarn 1
    rm mdout.mdp
    gmx mdrun -ntmpi 1 -ntomp 8 -pin off -deffnm NVT
    rm *#
    echo "Done"
fi

#get electrostatic potential
gmx potential -f NVT.xtc -s NVT.tpr -n index.ndx -sl 2000 -o potential_$varnumberofions_name.xvg << EOF
0
q
EOF
#get number density profile
gmx density -f NVT.xtc -n index.ndx -s NVT.tpr -sl 1000 -dens number -b 10 -o density_$varnumberofions_name.xvg << EOF
3
EOF

rm *#
rm step*

pwd
