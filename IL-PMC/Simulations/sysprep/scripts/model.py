# Load modules
from ase import Atom, Atoms
from ase.io import write, read
from ase.build import add_adsorbate, fcc100, fcc110, fcc111, fcc211, molecule

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
kden =  SED_KPDENSITY_SED
ads  =  read(f'{ion}.xyz')
name =f'{met}_{face}_{size[0]},{size[1]},{size[2]}_{ion}'
div  =  16
lat  =  grid*div*(2**0.5)
left =  6.0

# Create slab
atoms= fcc111(met, size=(size[0],size[1],size[2]), a=lat, orthogonal=False, vacuum=left)

# Add ion
add_adsorbate(atoms, ads, dist, 'ontop')

# Adjust model box to fit grid spacing
zmax = max([i[2] for i in atoms.get_positions()])
cell = atoms.get_cell()
cell[2][2] = grid*div*round((zmax+vac)/grid/div)
atoms.set_cell(cell)

# Save xyz
write(f'{name}.xyz', atoms)
