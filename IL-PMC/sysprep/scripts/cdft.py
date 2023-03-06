### Load modules
from ase import Atom, Atoms
from ase.calculators.mixing import SumCalculator
from ase.constraints import FixAtoms
from ase.dft.bee import BEEFEnsemble
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.optimize import QuasiNewton
from ase.parallel import parprint
from ase.spacegroup.symmetrize import FixSymmetry
from ase.units import Bohr
from dftd4.ase import DFTD4
from gpaw import GPAW, PW, restart
from gpaw.utilities import h2gpts
from gpaw.cdft.cdft import CDFT
import csv
import matplotlib.pyplot as plt
import numpy as np
import os

### Set variables
# Variables are substituted by NaRIBaS
ion  = 'SED_ION_SED'
dist =  SED_DISTANCE_SED
met  = 'SED_SLAB_SED'
size = [SED_XY_SED,SED_Z_SED]
face = 'SED_FACE_SED'
vac  =  SED_VACUUM_SED
mode = 'SED_MODE_SED'
func = 'SED_FUNCTIONAL_SED'
disp = 'SED_VDW_SED'
ctff =  SED_CUTOFF_SED
grid =  SED_GRIDSPACE_SED
klen =  SED_KPDENSITY_SED
name =f'{met}_{face}_{size[0]},{size[1]},{size[2]}_{ion}'
div  =  16
lat  =  grid*div*(2**0.5)
left =  6.0
sufx = 'opt'
char = [1,2,3,4]
lcao = 'lcao'
type = 'sys'

########################
### Run calculations ###
########################

for val in char:
    cdft = 'cdft_'+str(val)

    # Restart from traj file
    if os.path.exists(f'{name}_{sufx}.gpw') == True and os.stat(f'{name}_{sufx}.gpw').st_size > 0:
        if os.path.exists(f'{name}_{sufx}.xyz') == True and os.stat(f'{name}_{sufx}.xyz').st_size > 0:
            atoms = read(f'{name}_{sufx}.xyz')
            parprint(f'Start with optimized geometry.')
        else:
            parprint(f'Can not find optimized xyz file.')
    # Notify about incomplete optimization
    else:
        parprint(f'Optimization is not done.')

    ########################
    ### CDFT calculation ###
    ########################

    # Get atoms and calc
    atoms, calc = restart(f'{name}_{lcao}_{type}.gpw')

    # Set CDFT calculator
    ccalc = CDFT(atoms   =atoms,
                 calc    =calc,
                 charges =[float(val)],
                 charge_regions=[[atom.index for atom in atoms if atom.symbol != f'{met}']],
                 method  ='L-BFGS-B',
                 minimizer_options={'gtol':0.02},
                 txt     =f'{name}_{cdft}_{type}.txt',
                )

    # Get energy
    atoms.calc = ccalc
    atoms.get_potential_energy()    

    # Write xyz, gpw, and cube files
    calc.write(f'{name}_{cdft}_{type}.gpw')
    write(f'{name}_{cdft}_{type}.cube', atoms, data=calc.get_all_electron_density(gridrefinement=2)*Bohr**3)

    ########################
    ### El.  potential   ###
    ########################

    ### White and plot potential profile
    efermi = calc.get_fermi_level()
    v = calc.get_electrostatic_potential().mean(1).mean(0)
    z = np.linspace(0, calc.atoms.cell[2, 2], len(v), endpoint=False)
    plt.figure(figsize=(3.25, 3.25))
    plt.plot(z, v, label='xy-averaged potential')
    plt.plot([0, z[-1]], [efermi, efermi], label='Fermi level')
    plt.xlabel('Distance [Ang]')
    plt.ylabel('Electrostatic potential [V]')
    plt.xlim([0., z[-1]])
    plt.savefig(f'{name}_{cdft}_{type}.png', bbox_inches='tight')
    zv = list(zip(z,v))
    with open(f'{name}_{cdft}_{type}.pot', 'w', encoding='UTF8', newline='') as f:
        f.write(f'# EF/eV = {efermi}')
        f.write('\n')
        f.write(f'# Distance [A], Electrostatic potential [V]')
        f.write('\n')
        writer = csv.writer(f)
        writer.writerows(zv)
        f.close()

    del atoms
    # Delete CDFT calc
    del calc
