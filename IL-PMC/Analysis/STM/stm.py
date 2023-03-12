# Creates: 2d.png, 2d_I.png, line.png, dIdV.png
from ase.dft.stm import STM
from gpaw import GPAW, PW, FermiDirac
from ase.io import write, read
import matplotlib.pyplot as plt
import numpy as np
from ase.visualize.plot import plot_atoms


slab=read('Au_fcc111_6,6,4_C53N_sys.xyz')
slab.set_pbc([True, True, False])


calc = GPAW(mode='pw',
#            h=0.16,
            kpts=(1, 1, 1),
            symmetry='off',
            txt='au111.txt')
slab.calc = calc
energy = slab.get_potential_energy()
"""
calc.write('test.gpw', 'all')
calc = GPAW('test.gpw')
slab = calc.get_atoms()
"""

### Run STM calculation
stm = STM(slab)
bias = -2
z = 18.5#5
c = stm.get_averaged_current(bias, z)
x, y, h = stm.scan(bias, c, repeat=(11, 9))

#### Surface image
elSubsMulti = slab.repeat((12,10, 1))
f, (ax1) = plt.subplots(1, 1, sharey=False,sharex=False, figsize=(6, 6))
plot_atoms(elSubsMulti, ax1, radii=1.0)
ax1.set_xlim(120, 220)
ax1.set_ylim(27, 127)
ax1.set_xticks([])
ax1.set_yticks([])
f.savefig("./Surf.png", 
          dpi=600, bbox_inches="tight", pad_inches=0)

#### STM image
f, (ax2) = plt.subplots(1, 1, sharey=False,sharex=False, figsize=(3.33, 3.33), dpi=100)
cont = ax2.contourf(x, y, h, 40, cmap="gray")
ax2.set_xlabel(r"x / nm", fontsize=10)
ax2.set_ylabel(r"y / nm", fontsize=10)
cont.set_clim(17,19)
#10x10
ax2.set_xlim(60, 160)
ax2.set_ylim(0, 100)
ax2.set_xticks(np.arange(60, 161, 20))
ax2.set_xticklabels(np.arange(0,11, 2), fontsize=8)
ax2.set_yticks(np.arange(0, 101, 20))
ax2.set_yticklabels(np.arange(0,11, 2), fontsize=8)
ax2.set_aspect('equal')
f.savefig("STM.png", dpi=600)


