#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.Physics (v. 0.5):
    Provides methods to generate physical parameters, functions etc.
    
    Requirements: Python 2.2(?)->
                  numpy

    TODO: - add comments
          - generalise WritePales()
    
    Author: Martti Louhivuori (martti.louhivuori@helsinki.fi)

    Date: 25.04.2006

    ---

    Copyright (C) 2006-2008  Martti Louhivuori

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    The full license text can be found in the file LICENSE.
"""
from numpy import zeros, float, sum, array
import math
import os
from Mara.Math import InertiaMatrix

atomic_mass_unit = 1.66053886 * 10**-27 # kg
avogadros_number = 6.0221415 * 10**23 # mol^-1

# Gyromagnetic ratios according to Cavanagh et al. Protein NMR
# Spectroscopy (Academic Press, 1996)
# Note: H = 1H, N = 14N
gyromagnetic_constants = {
    'H' : 2.6752e8,
    'C' : 6.728e7,
    'N' : 1.934e7,
    'O' : -3.628e7,
    'F' : 2.5181e8,
    'Na' : 7.080e7,
    'P' : 1.0841e8,
    'Cd' : 5.934e7
    }

def calculate_rg(filename):
    """
    Calculates the radius of gyration of a molecule in a PDB (Cyana-format) 
    file using only alpha carbons.
    """
    s = os.system
    s("grep 'ATOM[ \t1-9]*CA[ \t]' %s > %s"%(filename,filename+'-rg'))
    file = open(filename+'-rg')
    lines = file.readlines()
    file.close()
    s('rm %s'%filename+'-rg')
    if not len(lines):
        return None
    atoms = zeros((len(lines),3), float)
    for i in range(len(lines)):
        p = lines[i].strip().split()
        atoms[i] = [float(p[5]),float(p[6]),float(p[7])]
    CoM = sum(atoms,0)/len(atoms)
    q = 0
    for atom1 in atoms:
        for atom2 in atoms:
            r = atom1 - atom2
            q += (r[0]**2 + r[1]**2 + r[2]**2)
    return q / (2*(len(atoms)+1)**2)

def CoM(molecule, whitelist=None):
    """
    Calculate the center of mass of a molecule.
    """
    try:
        atoms = molecule.getAtoms()
    except AttributeError:
        atoms = molecule
    if whitelist is None:
        accept = lambda x: True
    else:
        accept = lambda x: x in whitelist
    r = zeros((3), float)
    m = 0.0
    for atom in atoms:
        if accept(atom.name):
            r += array(atom.position) * float(atom)
            m += float(atom)
    return r / m

def InertiaTensor(molecule, unit_length=10**-10, unit_weight=atomic_mass_unit, \
        scale=False):
    """
    Calculate the inertia tensor of a molecule.
    """
    try:
        atoms = molecule.getAtoms()
    except AttributeError:
        atoms = molecule
    com = CoM(atoms)
    itensor = zeros((3,3), float)
    mass = 0.0
    for atom in atoms:
        itensor += InertiaMatrix(*((atom.position - com) * unit_length)) \
                * float(atom) * unit_weight
        mass += float(atom)
    if scale:
        return itensor / mass
    else:
        return itensor


class Potentials:
    """
    Electromagnetic potentials.
    """

    def __init__(self):
        pass

    def count_debye(self, ion, surface=0.47e18, temperature = 25):
        """
        Counts the debye length (Å) in a given ion concentration (mM). Returns 
        also ...
        """
        ion *= 1e-3
        e = -1.6022e-19
        n = 6.023e26
        E = 80
        E0 = 8.8542e-12
        kB = 1.38065e-23
        T = 273.15 + temperature
        K = math.sqrt(e**2*2*n*ion/(E*E0*kB*T))
        surf = surface*e
        A = 2*E*E0*kB*T/(e*surf)
        gamma = -K*A+math.sqrt(1+K**2*A**2)
        return (1/K)*10**10, K, gamma
    
    def debye(self,ion,cutoff,step=0.1e-10,surface=0.47e18,temperature=25):
        """
        Calculate the Debye potential for a given environment.
        """
        if surface > 0:
            sign = 1
        else:
            sign = -1
        r = 0
        (D,K,g) = self.count_debye(ion,sign*surface,temperature)
        positions = []
        potential = []
        while r < cutoff+step:
            p = sign*2*math.log((1+g*math.exp(-r*K))/(1-g*math.exp(-r*K)))
            positions.append(r)
            potential.append(p)
            r += step
        return potential, positions

    def WritePales(self,filename,positions,potential):
        """
        Write a Pales potential file.
        """
        file = open(filename,'w')
        for pos,pot in zip(positions,potential):
            str = " %.5E  %.5E\n" % (pos*1e9,pot)
            file.write(str)
        file.close()
        return True

