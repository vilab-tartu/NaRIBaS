#!/bin/bash

# ----------------------------
# 1. Variable Setup & Data Extraction
# ----------------------------

# (Your existing setup: get values from arrays)
electrode_name=$(echo ${current_electrode[0]} | awk '{print $1}')
xbox_nm=$(echo ${current_electrode[0]} | awk '{print $2}')
ybox_nm=$(echo ${current_electrode[0]} | awk '{print $3}')
zbox_nm=$(echo ${current_electrode[0]} | awk '{print $4}')
r_wall_nm=$(echo ${current_electrode[0]} | awk '{print $5}')
electrodeatoms=$(echo ${current_electrode[0]} | awk '{print $6}')
ion_name=$(echo ${current_ion[0]} | awk '{print $1}')
dens_ion_nm2=$(echo ${current_ion[0]} | awk '{print $2}')
ion_r_name=$(echo ${current_ion[0]} | awk '{print $3}')
ion_charge_name=$(echo ${current_ion[0]} | awk '{print $4}')
ion_centre_name=$(echo ${current_ion[0]} | awk '{print $5}')
numberofions_name=$(echo ${current_numberofions[0]} | awk '{print $1}')
temperature_name=$(echo ${current_temperature[0]} | awk '{print $1}')
version_name=$(echo ${current_version[0]} | awk '{print $1}')
replica_name=$(echo ${current_replica[0]} | awk '{print $1}')

# Calculations
xbox=$(awk "BEGIN {print $xbox_nm*10.0}" /dev/null)
ybox=$(awk "BEGIN {print $ybox_nm*10.0}" /dev/null)
zbox=$(awk "BEGIN {print $zbox_nm*10.0}" /dev/null)
xbox12=$(awk "BEGIN {print $xbox*0.5}" /dev/null)
ybox12=$(awk "BEGIN {print $ybox*0.5}" /dev/null)
zbox12=$(awk "BEGIN {print $zbox*0.1}" /dev/null)
r_wall=$(awk "BEGIN {print $r_wall_nm*10.0}" /dev/null)
z1_ion=$(awk "BEGIN {print $zbox12+$r_wall}" /dev/null)
z2_ion=$(awk "BEGIN {print $zbox12+20+$r_wall}" /dev/null)
A_box=$(awk "BEGIN {print $xbox*$ybox}" /dev/null)
A_ion=$(awk "BEGIN {print 100/$dens_ion_nm2}" /dev/null)
maxnumberofions_name=$(awk "BEGIN {print int($A_box/$A_ion)}" /dev/null)
varnumberofions_name=$(awk "BEGIN {print $maxnumberofions_name+1.0*$numberofions_name}" /dev/null)
electrode_atom_charge=$(awk "BEGIN {printf \"%.20f\n\",-1.0*$ion_charge_name*$varnumberofions_name/$electrodeatoms}" /dev/null)

# Output directory
fullpath_start=$electrode_name/$ion_name/$temperature_name/$version_name/$replica_name/$varnumberofions_name
mkdir -p $dir_experiments/$fullpath_start

# ----------------------------
# 2. Main Replica Loop
# ----------------------------

cd $dir_experiments/$fullpath_start/

success=0
while [ $success -eq 0 ]; do
    # Generate a random 7-digit seed
    seed=$(( (RANDOM * 32768 + RANDOM) % 9000000 + 1000000 ))
    echo $seed > seed.txt

    # Clean up all intermediate files but not NVT_lowtimestep.gro
    rm -f packmol.gro STEEP* mdout.mdp \#* *.trr index.ndx 0_STEEP.mdp 1_NVT_lowtimestep.mdp electrode_local.itp topol_local.top

    # ----------------------
    # System Preparation
    # ----------------------
    cd $dir_systempreparation

    sed 's+SED_dir_systempreparation_SED+'$dir_systempreparation'+g' packmol_Kappa.inp > packmol.inp
    sed -i 's/SED_seed_SED/'$seed'/g' packmol.inp
    sed -i 's/SED_ion_name_SED/'$ion_name'/g' packmol.inp
    sed -i 's/SED_ion_num_SED/'$varnumberofions_name'/g' packmol.inp
    sed -i 's/SED_xbox_SED/'$xbox'/g' packmol.inp
    sed -i 's/SED_ybox_SED/'$ybox'/g' packmol.inp
    sed -i 's/SED_zbox_left_SED/'$z1_ion'/g' packmol.inp
    sed -i 's/SED_zbox_right_SED/'$z2_ion'/g' packmol.inp
    sed -i 's/SED_xbox12_SED/'$xbox12'/g' packmol.inp
    sed -i 's/SED_ybox12_SED/'$ybox12'/g' packmol.inp
    sed -i 's/SED_zbox12_SED/'$zbox12'/g' packmol.inp
    sed -i 's/SED_electrode_name_SED/'$electrode_name'/g' packmol.inp

    # Packmol run
    packmol < packmol.inp
    gmx editconf -f packmol.pdb -o packmol.gro
    rm -f packmol.pdb packmol.inp

    # Edit box size for vacuum
    sed -i '$d' packmol.gro
    echo $xbox_nm $ybox_nm $zbox_nm >> packmol.gro

    # Index file
    gmx make_ndx -f packmol.gro -o index.ndx <<-EOF
    keep 0
    r $electrode_name
    name 1 Electrode
    r $ion_r_name
    name 2 Ion
    a $ion_centre_name
    name 3 Centre
    q
EOF

    # MDP files
    sed 's/SED_temperature_name_SED/'$temperature_name'/g' 0_STEEP_varTemp.mdp > 0_STEEP.mdp
    sed 's/SED_temperature_name_SED/'$temperature_name'/g' 1_NVT_lowtimestep_varTemp.mdp > 1_NVT_lowtimestep.mdp

    # Topology
    cp $dir_systempreparation/top/$electrode_name.itp electrode_local.itp
    sed -i 's/SED_electrode_charge_SED/'$electrode_atom_charge'/g' electrode_local.itp
    sed 's+SED_dir_systempreparation_SED+'$dir_systempreparation'+g' $dir_systempreparation/topol_local_Kappa.top > topol_local.top
    sed -i 's/SED_ion_name_SED/'$ion_name'/g' topol_local.top
    sed -i 's/SED_ion_r_name_SED/'$ion_r_name'/g' topol_local.top
    sed -i 's/SED_ion_num_SED/'$varnumberofions_name'/g' topol_local.top
    sed -i 's/SED_electrode_SED/'$electrode_name'/g' topol_local.top
    sed -i 's/SED_electrodeatoms_SED/'$electrodeatoms'/g' topol_local.top

    # Copy all prepared files to working directory
    cp 0_STEEP.mdp 1_NVT_lowtimestep.mdp electrode_local.itp topol_local.top packmol.gro index.ndx $dir_experiments/$fullpath_start/

    cd $dir_experiments/$fullpath_start/

    # ----------------------
    # Simulation Run Steps
    # ----------------------

    # 5x STEEP
    prevfile=packmol.gro
    for i in {1..5}; do
        gmx grompp -f 0_STEEP.mdp -c $prevfile -p topol_local.top -n index.ndx -o STEEP -maxwarn 1
        if [ ! -f STEEP.tpr ]; then
            echo "GROMPP failed for STEEP at attempt $i"
            break
        fi
        rm -f mdout.mdp
        gmx mdrun -ntmpi 1 -ntomp 8 -pin off -deffnm STEEP
        if [ ! -f STEEP.gro ]; then
            echo "MDRUN failed for STEEP at attempt $i"
            break
        fi
        prevfile=STEEP.gro
    done
    rm -f \#*

    # NVT equilibration
    gmx grompp -f 1_NVT_lowtimestep.mdp -c STEEP.gro -p topol_local.top -n index.ndx -o NVT_lowtimestep -maxwarn 1
    rm -f mdout.mdp
    gmx mdrun -ntmpi 1 -ntomp 8 -pin off -deffnm NVT_lowtimestep
    rm -f *trr \#*

    # ----------------------
    # Success Check
    # ----------------------

    if [ -f "NVT_lowtimestep.gro" ]; then
        echo "################################"
        echo "############ DONE! #############"
        echo "################################"
        success=1
    else
        echo "FAILED with seed $seed. Retrying with new seed..."
    fi

done

echo "Replica complete in $dir_experiments/$fullpath_start with seed $(cat seed.txt)"

