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
import cube_tools as ct
import matplotlib.pyplot as plt
import numpy as np
import os
import math

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

###Vector operations
def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))
def length(v):
    return math.sqrt(dotproduct(v, v))
def angle(v1, v2):
    return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))
def area(v1, v2):
    return np.sin(angle(v1, v2))*length(v1)*length(v2)

for type in ['sys','ads','met']:
    for work in ['', '_lcao', '_cdft']:
        # Write electron density profile
        if os.path.exists(f'{name}{work}_{type}.cube') == True and os.stat(f'{name}{work}_{type}.cube').st_size > 0:
            cb = ct.cube(f'{name}{work}_{type}.cube')
            dz = cb.Z[2]
            #ct.planar_average returns sum over ED of a plane at a given z in e/Angs, to get average in units of Angs3 additional division with unit cell plane area is necessary

            z,d= zip(*[(z,d) for z,d in cb.planar_average('z')])
            d  = np.array(d)/(area(np.asarray(cb.X)*cb.NX*Bohr,np.asarray(cb.Y)*cb.NY*Bohr)) #e/Angs3
            zd = [["# Distance [A]", "Electronic density [e/A^3]"]]
            zd += list(zip(z,d))
            with open(f'{name}{work}_{type}.den', 'w', encoding='UTF8', newline='') as f:
                writer = csv.writer(f)
                writer.writerows(zd)
                f.close()
        else:
            parprint(f'Can not find {name}{work}_{type}.cube file')
