# Load modules
from ase import Atom, Atoms
from ase.db import connect
from ase.io import read, write
from ase.parallel import parprint
from ase.spacegroup.symmetrize import FixSymmetry
from ase.units import Bohr
from gpaw import GPAW, restart
from sqlite3 import OperationalError
from time import sleep
import csv
#import cube_tools as ct
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import math

###Vector operations
def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))

def length(v):
    return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
    return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def estimateArea(v1, v2):
    return np.sin(angle(v1, v2))*length(v1)*length(v2)

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
path = 'SED_DBPATH_SED'
ver  = 'SED_VERSION_SED'
name =f'{met}_{face}_{size[0]},{size[1]},{size[2]}_{ion}'
mthd =f'{mode}_{ctff}_{func}_{disp}_{grid}_{klen}'
div  =  16
lat  =  grid*div*(2**0.5)
left =  6.0
sufx = 'sys'
frag = ['sys','ads','met']

for type in frag:

    if type=='sys':

        atoms,calc = restart(f'{name}_{type}.gpw', txt=None)
        atoms.calc = calc        # Read atoms

        # Potential profile for system in LCAO
        potential_profile = pd.read_csv(f'{name}_lcao_{type}.pot', sep=",", comment='#', skiprows=2)
        sys_lcao_pot_zax  = potential_profile.iloc[:, 0].values.tolist()
        sys_lcao_pot_val  = potential_profile.iloc[:, 1].values.tolist()

        # Potential profile for system in CDFT
        potential_profile = pd.read_csv(f'{name}_cdft_{type}.pot', sep=",", comment='#', skiprows=2)
        sys_cdft_pot_zax  = potential_profile.iloc[:, 0].values.tolist()
        sys_cdft_pot_val  = potential_profile.iloc[:, 1].values.tolist()

        # Potential profile for system in PW
        potential_profile = pd.read_csv(f'{name}_{type}.pot', sep=",", comment='#', skiprows=2)
        sys_pot_zax  = potential_profile.iloc[:, 0].values.tolist()
        sys_pot_val  = potential_profile.iloc[:, 1].values.tolist()

        # Density profile for system in LCAO
        density_profile = pd.read_csv(f'{name}_lcao_{type}.den', sep=",", comment='#', skiprows=1)
        sys_lcao_den_zax  = density_profile.iloc[:, 0].values.tolist()
        sys_lcao_den_val  = density_profile.iloc[:, 1].values.tolist()

        # Density profile for system in CDFT
        density_profile = pd.read_csv(f'{name}_cdft_{type}.den', sep=",", comment='#', skiprows=1)
        sys_cdft_den_zax  = density_profile.iloc[:, 0].values.tolist()
        sys_cdft_den_val  = density_profile.iloc[:, 1].values.tolist()

        # Density profile for system in PW
        density_profile = pd.read_csv(f'{name}_{type}.den', sep=",", comment='#', skiprows=1)
        sys_den_zax  = density_profile.iloc[:, 0].values.tolist()
        sys_den_val  = density_profile.iloc[:, 1].values.tolist()

        # Energies for system
        with open(f'{name}_{type}.nrg', "r") as nrgs_file:
            nrgs = nrgs_file.readlines()
        sys_epot_val = nrgs[1]
        sys_disp_val = nrgs[3]
        sys_bee_vals = nrgs[7:]

        # DDEC analysis in LCAO
        num_atoms    = open(f'{name}_lcao_ddec.xyz', 'r').readlines()[0]
        # ! skipping 18 lines of text plus the number of atoms
        ddf = pd.read_csv(f'{name}_lcao_ddec.xyz', delimiter="\s+", skiprows=18+int(num_atoms), nrows=int(num_atoms), header=None)
        ddf.columns  = ["numner", "symbol", "xaxis", "yaxis", "zaxis", "charge", "xdipole", "ydipole", "zdipole", "dipole", "Qxy", "Qxz", "Qyz", "Qx2y2", "Qz2", "", "", ""]       
        sys_lcao_chg_val  = ddf["charge"].tolist()
        sys_lcao_dip_val  = ddf["zdipole"].tolist()

        # DDEC analysis in CDFT
        num_atoms    = open(f'{name}_cdft_ddec.xyz', 'r').readlines()[0]
        # ! skipping 18 lines of text plus the number of atoms
        ddf = pd.read_csv(f'{name}_cdft_ddec.xyz', delimiter="\s+", skiprows=18+int(num_atoms), nrows=int(num_atoms), header=None)
        ddf.columns  = ["numner", "symbol", "xaxis", "yaxis", "zaxis", "charge", "xdipole", "ydipole", "zdipole", "dipole", "Qxy", "Qxz", "Qyz", "Qx2y2", "Qz2", "", "", ""]       
        sys_cdft_chg_val  = ddf["charge"].tolist()
        sys_cdft_dip_val  = ddf["zdipole"].tolist()

        # DDEC analysis in PW
        num_atoms    = open(f'{name}_ddec.xyz', 'r').readlines()[0]
        # ! skipping 18 lines of text plus the number of atoms
        ddf = pd.read_csv(f'{name}_ddec.xyz', delimiter="\s+", skiprows=18+int(num_atoms), nrows=int(num_atoms), header=None)
        ddf.columns  = ["numner", "symbol", "xaxis", "yaxis", "zaxis", "charge", "xdipole", "ydipole", "zdipole", "dipole", "Qxy", "Qxz", "Qyz", "Qx2y2", "Qz2", "", "", ""]       
        sys_chg_val  = ddf["charge"].tolist()
        sys_dip_val  = ddf["zdipole"].tolist()

        # Estimate surface charge density
        cell     = atoms.get_cell()
        ddf.sort_values(by=['zaxis'], inplace=True)
        metdf    = ddf.loc[ddf['symbol'].isin([met])]
        metchg   = sum(metdf['charge'].tolist())
        area     = estimateArea(cell[0],cell[1])
        sigma    = metchg/area #Keeping units to e/Angs2

        # Estimate distance
        adsdf    = ddf.loc[ddf['symbol'].isin(['B','C','N','S'])]
        ads_dist = adsdf['zaxis'].min() - metdf['zaxis'].max()

    elif type=='ads':

        # Energies for ion
        with open(f'{name}_{type}.nrg', "r") as nrgs_file:
            nrgs = nrgs_file.readlines()
        ads_epot_val  = nrgs[1]
        ads_disp_val  = nrgs[3]
        ads_bee_vals  = nrgs[7:]

        # Density profile for ion in LCAO
        density_profile = pd.read_csv(f'{name}_lcao_{type}.den', sep=",", comment='#', skiprows=1)
        ads_lcao_den_zax   = density_profile.iloc[:, 0].values.tolist()
        ads_lcao_den_val   = density_profile.iloc[:, 1].values.tolist()

        # Density profile for ion in PW
        density_profile = pd.read_csv(f'{name}_{type}.den', sep=",", comment='#', skiprows=1)
        ads_den_zax   = density_profile.iloc[:, 0].values.tolist()
        ads_den_val   = density_profile.iloc[:, 1].values.tolist()

    elif type=='met':

        # Energies for ion
        with open(f'{name}_{type}.nrg', "r") as nrgs_file:
            nrgs = nrgs_file.readlines()
        met_epot_val  = nrgs[1]
        met_disp_val  = nrgs[3]
        met_bee_vals  = nrgs[7:]

        # Density profile for ion in LCAO
        density_profile = pd.read_csv(f'{name}_lcao_{type}.den', sep=",", comment='#', skiprows=1)
        met_lcao_den_zax   = density_profile.iloc[:, 0].values.tolist()
        met_lcao_den_val   = density_profile.iloc[:, 1].values.tolist()

        # Density profile for ion in PW
        density_profile = pd.read_csv(f'{name}_{type}.den', sep=",", comment='#', skiprows=1)
        met_den_zax   = density_profile.iloc[:, 0].values.tolist()
        met_den_val   = density_profile.iloc[:, 1].values.tolist()

    else:
        pass

ads_enrg = float(sys_epot_val) - float(ads_epot_val) - float(met_epot_val)
ads_disp = float(sys_disp_val) - float(ads_disp_val) - float(met_disp_val)

# Data agregation
data_set ={'sys_pot_zax':  sys_pot_zax,
           'sys_pot_val':  sys_pot_val,
           'sys_chg_val':  sys_chg_val,
           'sys_dip_val':  sys_dip_val,
           'sys_den_zax':  sys_den_zax,
           'sys_den_val':  sys_den_val,
           'sys_bee_vals':sys_bee_vals,
           'ads_den_zax':  ads_den_zax,
           'ads_den_val':  ads_den_val,
           'ads_bee_vals':ads_bee_vals,
           'met_den_zax':  met_den_zax,
           'met_den_val':  met_den_val,
           'met_bee_vals':met_bee_vals,
           'sys_lcao_pot_zax':  sys_lcao_pot_zax,
           'sys_lcao_pot_val':  sys_lcao_pot_val,
           'sys_lcao_chg_val':  sys_lcao_chg_val,
           'sys_lcao_dip_val':  sys_lcao_dip_val,
           'sys_lcao_den_zax':  sys_lcao_den_zax,
           'sys_lcao_den_val':  sys_lcao_den_val,
           'ads_lcao_den_zax':  ads_lcao_den_zax,
           'ads_lcao_den_val':  ads_lcao_den_val,
           'met_lcao_den_zax':  met_lcao_den_zax,
           'met_lcao_den_val':  met_lcao_den_val,
           'sys_cdft_pot_zax':  sys_cdft_pot_zax,
           'sys_cdft_pot_val':  sys_cdft_pot_val,
           'sys_cdft_chg_val':  sys_cdft_chg_val,
           'sys_cdft_dip_val':  sys_cdft_dip_val,
           'sys_cdft_den_zax':  sys_cdft_den_zax,
           'sys_cdft_den_val':  sys_cdft_den_val,
          }

# Write to database
while True:
    try:
        with connect(f'{path}/psc.db') as db:
            checkstr=f'{name}_{mthd}'
            parprint(checkstr)
            # Check id
            try:
                id = db.get(check=checkstr).id
                db.update(id,
                          atoms,
                          version   = str(ver),
                          electrode = str(met),
                          ads_name  = str(ion),
                          ads_enrg  = round(ads_enrg,3),
                          ads_disp  = round(ads_disp,3),
                          ads_chg   = round(-metchg, 3),
                          ads_dist  = round(ads_dist,2),
                          check     = str(checkstr),
                          area      = round(area,2),
                          face      = str(face),
                          mode      = str(mode),
                          size      = str(size),
                          sigma     = sigma,
                          data      = data_set,
                         )
            except KeyError:
                db.write(atoms,
                         version   = str(ver),
                         electrode = str(met),
                         ads_name  = str(ion),
                         ads_enrg  = round(ads_enrg,3),
                         ads_disp  = round(ads_disp,3),
                         ads_chg   = round(-metchg, 3),
                         ads_dist  = round(ads_dist,2),
                         check     = str(checkstr),
                         area      = round(area,2),
                         face      = str(face),
                         mode      = str(mode),
                         size      = str(size),
                         sigma     = sigma,
                         data      = data_set,
                        )
        break
    except OperationalError:
        parprint('Can not connect to the database!')
        sleep(11)
