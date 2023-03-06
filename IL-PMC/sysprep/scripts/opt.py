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
import csv
import matplotlib.pyplot as plt
import os
import numpy as np

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

def getkpts(klen,size):
     kptx = -(-klen//size[0])
     kpty = -(-klen//size[1])
     kptz = 1
     return tuple([kptx,kpty,kptz])

def getgpts(grid,cell):
     gptx = round(cell[0][0]//grid/div)*div
     gpty = round(cell[0][0]//grid/div)*div
     gptz = round(cell[2][2]//grid/div)*div
     return tuple([gptx,gpty,gptz])


########################
### Run optimization ###
########################

### Check if optimization is done
if os.path.exists(f'{name}_{sufx}.gpw') == True and os.stat(f'{name}_{sufx}.gpw').st_size > 0:
    parprint(f'{name} optimization done -- gpw file exists.')
### Run optimization
else:
    ### Read atoms
    if os.path.exists(f'{name}_{sufx}_step.gpw') == True and os.stat(f'{name}_{sufx}_step.gpw').st_size > 0:
        atoms,calc = restart(f'{name}_{sufx}_step.gpw', txt=None)
    elif os.path.exists(f'{name}_{sufx}_step.traj') == True and os.stat(f'{name}_{sufx}_step.traj').st_size > 0:
        atoms = read(f'{name}_{sufx}_step.traj',-1)
        parprint(f'Start with traj geometry.')
    else:
        atoms = read(f'{name}.xyz')
        parprint(f'Start with xyz geometry.')

    ### Set constraints
    c = []
    c.append(FixAtoms(mask=[atom.tag > 2 for atom in atoms]))
    #c.append(FixSymmetry(atoms))
    atoms.set_constraint(c)

    ### Set calculator
    calc = GPAW(basis='dzp',
                gpts=getgpts(grid,atoms.get_cell()),
                kpts=getkpts(klen,size),
                mode=PW(ctff),
                parallel={'augment_grids':True,'sl_auto':True,'use_elpa':True},
                poissonsolver={'dipolelayer':'xy'},
                txt=f'{name}_{sufx}.txt',
                xc=func,
               )

    ### Add dispersion correction
    if disp == 'D4':
        atoms.calc = SumCalculator([DFTD4(method=func), calc])
    else:
        atoms.calc = calc

    ### Optimize
    opt = QuasiNewton(atoms, trajectory=f'{name}_{sufx}.traj', logfile=f'{name}_{sufx}.log')
    traj= Trajectory(f'{name}_{sufx}_step.traj', 'a', atoms)
    opt.attach(traj.write, interval=1)
    def writegpw():
        calc.write(f'{name}_{sufx}_step.gpw')
    opt.attach(writegpw, interval=1)
    opt.run(fmax=0.05, steps=23)

    ### Write xyz, gpw, and cube files
    atoms.write(f'{name}_{sufx}.xyz')
    calc.write(f'{name}_{sufx}.gpw')
