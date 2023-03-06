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
#from gpaw.cdft.cdft import CDFT
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
frag = ['sys','ads','met']
cdft = 'cdft'
lcao = 'lcao'

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

def xyzh(atoms, filename):
    from gpaw.analyse.hirshfeld import HirshfeldPartitioning
    f = open('{0}.xyz'.format(filename), 'w') # set xmol format
    a = atoms.get_chemical_symbols()          # get chemical symbols
    p = atoms.get_positions()                 # get positions of the atoms
    h = HirshfeldPartitioning(atoms.calc).get_charges()
    f.write('{0}\n'.format(len(a)))
    f.write('Properties=species:S:1:pos:R:3:charge:R:1\n') # ensure compatibility with ASE
    for i in range(len(a)):                   # print symbol, positions, and charge
        s = '{0}'.format(a[i])
        x = '{0:.6f}'.format(round(p[i][0],6))
        y = '{0:.6f}'.format(round(p[i][1],6))
        z = '{0:.6f}'.format(round(p[i][2],6))
        c = '{0:.3f}'.format(round(h[i],3))
        f.write('{0:<4}{1:>16}{2:>16}{3:>16}{4:>10}\n'.format(s,x,y,z,c))
    f.close()

########################
### Run calculations ###
########################

for type in frag:
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

    if type == 'ads':
        # for adsorbate  only remove all metals
        del atoms[[atom.index for atom in atoms if atom.symbol == f'{met}']]
    elif type == 'met':
        # for metal slab only remove all non-metals
        del atoms[[atom.index for atom in atoms if atom.symbol != f'{met}']]
    else:
        pass

    ########################
    ### PW   calculation ###
    ########################

    # Set PW calculator
    calc = GPAW(basis   ='dzp',
                gpts    =getgpts(grid,atoms.get_cell()),
                kpts    =getkpts(klen,size),
                mode    =PW(ctff),
                parallel={'augment_grids':True,'sl_auto':True,'use_elpa':True},
                poissonsolver={'dipolelayer':'xy'},
                txt     =f'{name}_{type}.txt',
                xc      =func,
               )

    # Set calculator and get energies
    atoms.calc = calc
    pw_e = atoms.get_potential_energy()
    bee  = BEEFEnsemble(calc).get_ensemble_energies()
    
    ### Write xyz, gpw, and cube files
    atoms.write(f'{name}_{type}.xyz')
    calc.write(f'{name}_{type}.gpw')
    write(f'{name}_{type}.cube', atoms, data=calc.get_all_electron_density(gridrefinement=2)*Bohr**3)

    if type == 'sys':
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
        plt.savefig(f'{name}_{type}.png', bbox_inches='tight')
        zv = list(zip(z,v))
        with open(f'{name}_{type}.pot', 'w', encoding='UTF8', newline='') as f:
            f.write(f'# EF/eV = {efermi}')
            f.write('\n')
            f.write(f'# Distance [A], Electrostatic potential [V]')
            f.write('\n')
            writer = csv.writer(f)
            writer.writerows(zv)
            f.close()
    else:
        pass
    # Delete PW calc
    del calc

    ########################
    ### LCAO calculation ###
    ########################

    # Set LCAO calculator
    calc = GPAW(basis   ='dzp',
                gpts    =getgpts(grid,atoms.get_cell()),
                kpts    =getkpts(klen,size),
                mode    ='lcao',
                parallel={'augment_grids':True,'sl_auto':True,'use_elpa':True},
                poissonsolver={'dipolelayer':'xy'},
                spinpol =True,
                txt     =f'{name}_{lcao}_{type}.txt',
                xc      =func,
               )

    # Set calculator and get energies
    atoms.calc = calc
    lcao_e = atoms.get_potential_energy()

    # Write xyz, gpw, and cube files
    calc.write(f'{name}_{lcao}_{type}.gpw')
    write(f'{name}_{lcao}_{type}.cube', atoms, data=calc.get_all_electron_density(gridrefinement=2)*Bohr**3)

    ########################
    ### El.  potential   ###
    ########################

    if type == 'sys':
        ### Write Hirshfeld charges
        xyzh(atoms, f'{name}_{lcao}_{type}_Hirshfeld')

        ### Write and plot potential profile
        efermi = calc.get_fermi_level()
        v = calc.get_electrostatic_potential().mean(1).mean(0)
        z = np.linspace(0, calc.atoms.cell[2, 2], len(v), endpoint=False)
        plt.figure(figsize=(3.25, 3.25))
        plt.plot(z, v, label='xy-averaged potential')
        plt.plot([0, z[-1]], [efermi, efermi], label='Fermi level')
        plt.xlabel('Distance [Ang]')
        plt.ylabel('Electrostatic potential [V]')
        plt.xlim([0., z[-1]])
        plt.savefig(f'{name}_{lcao}_{type}.png', bbox_inches='tight')
        zv = list(zip(z,v))
        with open(f'{name}_{lcao}_{type}.pot', 'w', encoding='UTF8', newline='') as f:
            f.write(f'# EF/eV = {efermi}')
            f.write('\n')
            f.write(f'# Distance [A], Electrostatic potential [V]')
            f.write('\n')
            writer = csv.writer(f)
            writer.writerows(zv)
            f.close()
    else:
        pass
    # Delete LCAO calc
    del calc

    ########################
    ### D4  calculation  ###
    ########################

    atoms.calc = DFTD4(method=func)
    d4_e = atoms.get_potential_energy()

    ########################
    ### Write  energies  ###
    ########################

    # Write energies
    with open(f'{name}_{type}.nrg', 'w', encoding='UTF8', newline='') as f:
        f.write(f'# Potential  energy for {func} [eV]')
        f.write('\n')
        f.write(f'{pw_e}')
        f.write('\n')
        f.write(f'# Dispersion energy for {disp} [eV]')
        f.write('\n')
        f.write(f'{d4_e}')
        f.write('\n')
        f.write(f'# LCAO energy for {func} [eV]')
        f.write('\n')
        f.write(f'{lcao_e}')
        f.write('\n')
        f.write(f'# BEE ensemple for {func} (eV)')
        f.write('\n')
        np.savetxt(f, bee)
        f.write('\n')
        f.close
