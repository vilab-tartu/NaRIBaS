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
# Note that in this step uses a different list for functionals to iterate over multiple functionals
# while an optimization is only performed with one functional
anion_name=$(echo ${current_anion[0]} | awk '{print $1}')
cation_name=$(echo ${current_cation[0]} | awk '{print $1}')
functional_name=$(echo ${current_functional_sp[0]} | awk '{print $1}' | sed 's,\/,\\/,g')
basisset_name=$(echo ${current_basisset[0]} | awk '{print $1}' | sed 's,\/,\\/,g')
version_name=$(echo ${current_version[0]} | awk '{print $1}')

# This time another variable is read in from the cation_list: the cation atom count
# It is defined in the lists as second element on each row and it is necessary for
# separating anions and cations for calculating interaction energies
cation_atom_num=$(echo ${current_cation[0]} | awk '{print $2}')

#These variables are now accessed as in the following example:
echo SYSMEM PARAMETERS
echo Anion $anion_name
echo Cation $cation_name
echo Functional $functional_name
echo Basis-set $basisset_name

# Define path for storing computational data
# Single point results are placed in different folders named after the method
fullpath=$(echo $cation_name$anion_name/$functional_name/ | sed 's,\\\/,,g')
mkdir -p $dir_experiments/$fullpath

# Preparing files for calculation
# The geometry is taken from the previous optimization step folder
cd $dir_experiments/$fullpath
cp $dir_experiments/$cation_name$anion_name/opt/$cation_name$anion_name.xyz ./

# Assigning ion pair name and ion lengths as variables
Natom=`awk 'NR==1 {print $0}' $cation_name$anion_name".xyz"` 
Natom=$(awk "BEGIN {print "$Natom+1"}" /dev/null)
Ncat=$(awk "BEGIN {print "$cation_atom_num+2"}" /dev/null)
k=2

# Assigning atoms of the xyz file to temporary files corresponding to pair, anion or cation
# This is necessary for the automatic separation of anion and cation geometries
# For this tho work, the cation atoms need to preceed the anion atoms in the xyz file
while [ "$k" -le "$Natom" ]
do
        
	k=$(awk "BEGIN {print "$k+1"}" /dev/null)
	if [ "$k" -le "$Ncat" ]
	then
		at=`awk 'NR=='$k' {print $1}' $cation_name$anion_name".xyz"`
		x=`awk 'NR=='$k' {print $2}' $cation_name$anion_name".xyz"`
		y=`awk 'NR=='$k' {print $3}' $cation_name$anion_name".xyz"`
		z=`awk 'NR=='$k' {print $4}' $cation_name$anion_name".xyz"`
		echo $at " " $x " " $y " " $z >> $cation_name$anion_name'_Cat'
		echo $at " " $x " " $y " " $z >> $cation_name$anion_name'_IP'
	else 
                at=`awk 'NR=='$k' {print $1}' $cation_name$anion_name".xyz"`
                x=`awk 'NR=='$k' {print $2}' $cation_name$anion_name".xyz"`
                y=`awk 'NR=='$k' {print $3}' $cation_name$anion_name".xyz"`
                z=`awk 'NR=='$k' {print $4}' $cation_name$anion_name".xyz"`
		echo $at " " $x " " $y " " $z >> $cation_name$anion_name'_An'
		echo $at " " $x " " $y " " $z >> $cation_name$anion_name'_IP'
	fi 
done

# As an example, the input files are created from scratch this time
for i in An Cat IP
do
	# Different fragments have different total charges
	if [ "$i" = "An" ] 
	then
		ch=$(awk "BEGIN {print "-1"}" /dev/null)
	fi 
	if [ "$i" = "Cat" ] 
	then
		ch=$(awk "BEGIN {print "1"}" /dev/null)
	fi
	if [ "$i" = "IP" ] 
	then
		ch=$(awk "BEGIN {print "0"}" /dev/null)
	fi

	# Inserting functional, basis set and other keywords to a new empty input file
	echo "!$functional_name $basisset_name TIGHTSCF CHELPG" | sed 's/___/ /g' | sed 's,\\\/,\/,g' >> $i.inp
	echo "!PAL4" >> $i.inp

	# Inserting fragment charge to input file
	echo "* xyz" $ch "1" >> $i.inp

        # Concatenating the fragment geometry from the previously created temporary file to input file
	cat $cation_name$anion_name'_'$i >> $i.inp
        rm $cation_name$anion_name'_'$i
	echo "*" >> $i.inp
	echo "" >> $i.inp
	echo "%scf MaxIter 300 end" >> $i.inp
	echo "% MaxCore 2000" >> $i.inp

	# Giving the input a systematic name
	mv $i.inp $cation_name$anion_name"_"$i.inp

	# Executing single point calculation
	$dir_orca/orca $cation_name$anion_name"_"$i.inp > $cation_name$anion_name"_"$i.out	
done




