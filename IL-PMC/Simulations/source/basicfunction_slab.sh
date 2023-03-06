#!/bin/bash

function calc_electrode_positions_in_slab()
{

# Function to calculate the electrode positions

# Read the number of ion pairs
while read inputline; do
numberofionpairs_name="$(echo $inputline | cut -d' ' -f1)"
done < $currentdir/$inputlists_folder/numberofionpairs.list

# Get the particle density
particledensity=0.0
while read inputline
do
   currenttemperature="$(echo $inputline | cut -d' ' -f1)"
   currentdens="$(echo $inputline | cut -d' ' -f2)"
   if [ $currenttemperature -eq $temperature_name ]; then
      particledensity=$currentdens
   fi
done < $currentdir/$inputlists_folder/particledensity_$cation_name\_$anion_name\_$additive_name\_$mol_fraction.list

# Complain and exit if particle density is not given in the list
if [ $(echo $particledensity 0.0 | awk '{if ($1 > $2) print "1"; else print "0";}') -eq 0 ]; then
   echo 'Error: Particle density not defined for the given temperature '$temperature'K.'
   exit
fi

# Calculate the volume the IL needs and divide by the x- and y-dimension which are already fixed by the electrode geometry.
exp="$numberofionpairs_name/($particledensity*$xbox_nm*$ybox_nm)*10.0" # Gromacs thinks in nm, packmol in A !!
zbox=$(awk "BEGIN {print $exp}" /dev/null) # in A

###################################
# To calculate the position of electrodes use the following algorithm:
#
# The |-lines are positions of the walls, ()-brackets show the size of the atoms/molecules that have to be taken into account when packing ions in the box.
#
# distances:              r_wall       2*r_ion                2*r_ion      r_wall
#                     |            )(     .     )..........(     .     )(           |
# abs. pos.:  pos_left_electrode   z1   z1_ion                z2_ion   z2     pos_right_electrode
#
#
# pos_left_electrode = r_wall*2*num_of_walls-rwall
# z1_accesiblesurface=z1 = pos_left_electrode+r_wall
# z1_ion = z1+r_ion
# z2_accesiblesurface = z2 = z1+zbox
# pos_right_electrode = z2+r_wall
# z2_ion = z2-r_ion
#
####################################

# Calculate the wall position
exp="$xbox_nm*10.0"
xbox=$(awk "BEGIN {print $exp}" /dev/null)  #in A

exp="$ybox_nm*10.0"
ybox=$(awk "BEGIN {print $exp}" /dev/null)  #in A

exp="$r_wall_nm*10.0"
r_wall=$(awk "BEGIN {print $exp}" /dev/null)  #in A
exp="$r_wall*2.0"
pos_left_electrode=$(awk "BEGIN {print $exp}" /dev/null)  #in A

if [ $(echo $r_cation_nm $r_anion_nm | awk '{if ($1 > $2) print "1"; else print "0";}') -eq 1 ]; then
   r_ion_nm=$r_cation_nm
else
   r_ion_nm=$r_anion_nm
fi
exp="$r_ion_nm*10.0"
r_ion=$(awk "BEGIN {print $exp}" /dev/null) #in A

exp="$pos_left_electrode+$r_wall"
z1=$(awk "BEGIN {print $exp}" /dev/null)
exp="$z1+$zbox"
z2=$(awk "BEGIN {print $exp}" /dev/null)
exp="$z2+$r_wall"
pos_right_electrode=$(awk "BEGIN {print $exp}" /dev/null)

exp="$z1+$r_ion"
z1_ion=$(awk "BEGIN {print $exp}" /dev/null)
exp="$z2-$r_ion"
z2_ion=$(awk "BEGIN {print $exp}" /dev/null)

# Some unit conversion
exp="$pos_left_electrode/10.0"
pos_left_electrode_nm=$(awk "BEGIN {print $exp}" /dev/null)  #in nm

exp="$pos_right_electrode/10.0"
pos_right_electrode_nm=$(awk "BEGIN {print $exp}" /dev/null)  #in nm

####################################

# A crosscheck if the dimensions of the box are physical. We have to account for the Yeh/Berkowitz correction for electrostatics in Slab geometry and add a vaccum slab to our simulation box. This slab should be larger than the maximum of box size in (x,y)-dimension. Moreover the total length of the box in z-direction should be larger than three times the maximum of box size in (x,y)-dimension.

if [ $(echo $xbox_nm $ybox_nm | awk '{if ($1 >= $2) print "1"; else print "0";}') -eq 1 ]; then
   z_vac_min=$xbox_nm
else
   z_vac_min=$ybox_nm
fi

exp="$pos_right_electrode/10.0+$r_wall_nm+$z_vac_min*1.5" #in nm
zbox_vacuum=$(awk "BEGIN {print $exp}" /dev/null)

}


function calc_slicing_of_slab()
{
# How many slices to take into account? 
# Call calc_electrode_positions_in_slab first to obtain zbox_vacuum. 
# Provide deltaz in nm before calling the function.
exp="$zbox_vacuum/$deltaz"
number=$(awk "BEGIN {print $exp}" /dev/null)
# Round the number as numberofslices need to be an integer
numberofslices=$(echo $number | awk '{printf("%d\n",$1 + 0.5)}')
}
