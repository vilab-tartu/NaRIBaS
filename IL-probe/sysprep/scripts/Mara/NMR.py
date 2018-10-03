#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.NMR (v. 0.5):
    NMR related methods and classes.
    
    Requirements: Python 2.2(?)->
                  numpy

    TODO: - add comments
          - clean obsolete stuff
          - fix relative maximums of dipoles
    
    Author: Martti Louhivuori (martti.louhivuori@helsinki.fi)

    Date: 25.04.2008

    ---

    Copyright (C) 2006-2008  Martti Louhivuori

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    The full license text can be found in the file LICENSE.
"""
import re
from numpy import array, float, inner, zeros
from string import center
from Mara.aux import BasicIterator
from Mara import loki

__dipole_types__ = ['HACA'] # OBSOLETE? 

def get_dipoles(pdb, saupe, types=[]): # OBSOLETE?
    from myPDB import PDB # FIX THIS!!!
    if not types:
        types = __dipole_types__
    molecule = PDB(pdb)
    chain = molecule.getChain()
    cosines = chain.getCosines(types)
    dipoles = zeros((len(cosines),5),float)
    for i in range(len(cosines)):
        x,y,z = cosines[i].values
        dipoles[i] = [y**2-x**2, z**2-x**2, 2*x*y, 2*x*z, 2*y*z]
    return sum(dipoles * array(saupe),1)

"""
# should generalise for ALL dipoles in this fashion...
    ro = re.compile('([A-Z]+)(\[([+-][0-9]+)\])?')
    for type in types:
        name, x, shift = ro.search(type)
        angles = chain.getDipoleAngles(type=name)
"""

def reconstructSaupe(saupe):
    """
    Construct Saupe matrix from the five independent components given as
    Szz, Syy-Sxx, Sxy, Sxz, Syz.
    """
    szz, sxy, sxz, syz = saupe[0], saupe[2], saupe[3], saupe[4]
    syy = -0.5*(saupe[1]+szz)
    sxx = -syy-szz
    return array([[sxx,sxy,sxz],[sxy,syy,syz],[sxz,syz,szz]],float)

def deconstructSaupe(saupe):
    """
    Reduce Saupe matrix to the five independent components in the form
    Szz, Syy-Sxx, Sxy, Sxz, Syz.
    """
    return [saupe[2,2], saupe[0,0]-saupe[1,1], saupe[0,1], saupe[0,2], \
            saupe[1,2]]

def rotateSaupe(rotation, saupe):
    """
    Rotate Saupe matrix using a rotation matrix.
    """
    if type(saupe) != type(array([])):
        saupe = reconstructSaupe(saupe)
    return inner(saupe,rotation)

def makeOrderMatrix(cosines):
    """
    Construct an Order matrix from a list of cosines.
    """
    matrix = zeros((len(cosines),5),float)
    for i,c in zip(range(len(cosines)), cosines):
        matrix[i] = c.values[1]**2 - c.values[0]**2, \
                c.values[2]**2 - c.values[0]**2, 2*c.values[0]*c.values[1], \
                2*c.values[0]*c.values[2], 2*c.values[1]*c.values[2]
    return matrix

def csum(couplings):
    """
    Prune a Coupling list by summing all similar couplings together.
    """
    s = Couplings()
    boycot = []
    for i in range(len(couplings)):
        if i not in boycot:
            if couplings.count(couplings[i].tag) > 1:
                x = float(couplings[i])
                for j in range(i+1,len(couplings)):
                    if couplings[i].tag == couplings[j].tag:
                        x += float(couplings[j])
                        boycot.append(j)
                s.append(couplings[i].replica(x))
            else:
                s.append(couplings[i])
    return s

def tags2couplings(tags, values=0.0):
    """
    Convert a tag list to Couplings.
    """
    s = Couplings()
    if type(values) in [int, float]:
        values = [values] * len(tags)
    for tag, value in zip(tags, values):
        i, name = tag
        c = Coupling(value, DCouplings[name], i)
        s.append(c)
    return s

def fuse_seq2couplings(seq, couplings):
    """
    Fuse sequence info into Couplings.
    """
    curr_id = None
    i = -1
    for c in couplings:
        if curr_id != int(c):
            i += 1
            curr_id = int(c)
        t = c.getType()
        c.residues = (seq[i + t.shift[0]], seq[i + t.shift[1]])
    return couplings

def renormalise_couplings(couplings):
    """
    Renormalise Couplings, i.e. multiply with coupling-specific gammas.
    """
    for c in couplings:
        c.value = float(c) * float(c.getType())
    return couplings


class DCoupling:
    def __init__(self, name, gamma, atoms, shift):
        self.__name__ = str(name)
        self.__gamma__ = float(gamma)
        self.atoms = tuple(atoms[0:2])
        self.shift = tuple(shift[0:2])
    def __str__(self):
        return self.__name__
    def __float__(self):
        return self.__gamma__
    def __int__(self):
        return int(self.__gamma__)
    def __repr__(self):
        main = 'Dipolar coupling %s.\n'%self.__name__
        gamma = 'gamma: %f\t'%self.__gamma__
        pos = 'position shifts: %s'%repr(self.shift)
        return '%s'%(main + gamma + center('atoms: %s\t'%repr(self.atoms), \
                70-len(gamma)-len(pos)) + pos)
    def __call__(self, x):
        return x/self.__gamma__

### FIX !!! ###
# Main chain, intra-residue dipoles
DCouplings = {
    'DHaCa': DCoupling('DHaCa', 44539.47, ('HA','CA'), (0,0)),
    'DHaCo': DCoupling('DHaCo', 1.0, ('HA','C'), (0,0)),
    'DCaCo': DCoupling('DCaCo', 4284.77, ('CA','C'), (0,0)),
    'DCaCb': DCoupling('DCaCb', 4175.99, ('CA','CB'), (0,0)),
    'DCbCg': DCoupling('DCbCg', 4097.46, ('CB','CG'), (0,0)),
    'DCgCd': DCoupling('DCgCd', 4014.15, ('CG','CD'), (0,0)),
    'DCdCe': DCoupling('DCdCe', 4193.01, ('CD','CE'), (0,0)),
    'DCeNz': DCoupling('DCeNz', -1916.37, ('CE','NZ'), (0,0)),
    'DHaN':  DCoupling('DHaN', 1.0, ('HA','N'), (0,0)),
    'DNCa':  DCoupling('DNCa', 1.0, ('N','CA'), (0,0)),
    'DHnN':  DCoupling('DHnN', -21585.2, ('HN','N'), (0,0)),
# Main chain dipoles btw adjacent residues
    'DCoN':  DCoupling('DCoN', -2609.05, ('C','N'), (0,1)),
    'DHnCo':  DCoupling('DHnCo', 6666.07, ('HN','C'), (1,0))
}
### FIX !!! ###
# Provide direct instances of each dipolar coupling.
DHaCa = DCouplings['DHaCa']
DHaCo = DCouplings['DHaCo']
DCaCo = DCouplings['DCaCo']
DCaCb = DCouplings['DCaCb']
DCbCg = DCouplings['DCbCg']
DCgCd = DCouplings['DCgCd']
DCdCe = DCouplings['DCdCe']
DCeNz = DCouplings['DCeNz']
DHaN = DCouplings['DHaN']
DNCa = DCouplings['DNCa']
DHnN = DCouplings['DHnN']
DCoN = DCouplings['DCoN']
DHnCo = DCouplings['DHnCo']

class CouplingIterator:
    def __init__(self, data):
        self.data = data
        self.index = -1
        self.max = len(data)-1
    def __iter__(self):
        return self
    def __next__(self):
        if self.index < self.max:
            self.index += 1
            return self.data[self.index]
        else:
            raise StopIteration

class Coupling:

    def __init__(self, value, type, id, residues=('???','???'), max=None, \
            error=0.0):
        from Mara.Library import omniTranslator
        self.value = float(value)
        self.__id__ = int(id)
        if type is str and DCouplings.has_key(type):
            self.__type__ = DCouplings[type]
        elif isinstance(type, DCoupling):
            self.__type__ = type
        else:
            raise TypeError, "Incompatible coupling type '%s'." % repr(type)
        if (isinstance(residues, tuple) or isinstance(residues, list)) \
                and len(residues) == 2 and \
                isinstance(residues[0], str) and isinstance(residues[1], str):
            self.residues = (
                    omniTranslator(residues[0], '3-letter-aminoacids'),
                    omniTranslator(residues[1], '3-letter-aminoacids'))
        else:
            raise TypeError, "Incompatible residue names '%s'." % repr(residues)
        self.error = float(error)
        self.tag = (self.__id__, str(self.__type__))
        if max:
            self.__max__ = float(max)
        else:
            self.__max__ = None
        self.__name__ = 'Coupling'

    def __div__(self, x):
        return self.replica(float(self)/x)

    def __idiv__(self, x):
        self.value /= x

    def __rdiv__(self, x):
        return self.replica(x/float(self))

    def __mul__(self, x):
        return self.replica(float(self)*x)

    def __imul__(self, x):
        self.value *= x

    def __rmul__(self, x):
        return self.replica(x*float(self))

    def __float__(self):
        return self.value

    def __int__(self):
        return self.__id__

    def __str__(self):
        return str(self.__type__)

    def __repr__(self):
        return 'Coupling(%f, %s, %d, %s, %s, %f)' % \
                (self.value, str(self.__type__), self.__id__,
                repr(self.residues), repr(self.__max__),
                self.error)

    def match(self, tag):
        return self.tag == tag

    def normalise(self):
        return float(self) / (self.__max__ or float(self.__type__))

    def type(self): #TODO: get rid of this...
        return self.__type__
    def getType(self):
        return self.__type__

    def replica(self, value):
        return Coupling(value, self.__type__, self.__id__, \
                residues=self.residues, max=self.__max__, error=self.error)

    def getResidues(self):
        return self.residues


class Couplings:

    def __init__(self, couplings=[]):
        self.__data__ = []
        self += couplings
        self.__name__ = 'Couplings'

    def __translate__(self, couplings):
        if '__name__' not in dir(couplings):
            if type(couplings) != type([]):
                couplings = [Coupling(*couplings)]
            else:
                new = []
                for c in couplings:
                    if '__name__' not in dir(c) or c.__name__ != 'Coupling':
                        new.append(Coupling(*c))
                    else:
                        new.append(c)
                couplings = new
        elif couplings.__name__ == 'Coupling':
            couplings = [couplings]
        else:
            couplings = list(couplings)
        return couplings

    def __iadd__(self, couplings):
        self.__data__ += self.__translate__(couplings)
        return self

    def __add__(self, couplings):
        return Couplings(self.__data__ + self.__translate__(couplings))

    def __div__(self, x):
        new = Couplings()
        for c in self:
            new.append(c/x)
        return new

    def __rdiv__(self, x):
        new = Couplings()
        for c in self:
            new.append(x/c)
        return new

    def __idiv__(self, x):
        for c in self:
            c /= x

    def __mul__(self, x):
        new = Couplings()
        for c in self:
            new.append(c*x)
        return new

    def __rmul__(self, x):
        new = Couplings()
        for c in self:
            new.append(x*c)
        return new

    def __imul__(self, x):
        for c in self:
            c *= x

    def __str__(self):
        return '\n'.join([str(x) for x in self.__data__])

    def __repr__(self):
        glue = ',\n'+' '*11
        return 'Couplings([%s])'%(glue.join([repr(x) for x in self.__data__]))

    def __iter__(self):
        return BasicIterator(self.__data__)

    def __getitem__(self, i):
        return self.__data__[i]

    def __len__(self):
        return len(self.__data__)

    def append(self, coupling):
        self.__iadd__(coupling)

    def extend(self, other):
        for c in other:
            self.__iadd__(c)

    def count(self, tag):
        i = 0
        for c in self:
            if c.tag == tag:
                i += 1
        return i

    def countType(self, ctype):
        i = 0
        for c in self:
            if c.type() == ctype:
                i += 1
        return i

    def flip(self, targets):
        new = Couplings()
        for c in self:
            if str(c.type()) in targets: s = -1
            else: s = 1
            new.append(c*s)
        return new

    def getTypes(self):
        types = []
        for c in self:
            if c.type() not in types:
                types.append(c.type())
        return types

    def index(self, coupling):
        return self.__data__.index(coupling)

    def insert(self): pass
    def pop(self): pass
    def remove(self, x):
        if type(x) == tuple and len(x) == 2:
            i = self.index(self.find(x))
        elif isinstance(Coupling, x):
            i = self.index(x)
        else:
            raise ArgumentError, "%s is neither a tag nor a Coupling."%repr(x)
        del self.__data__[i]
        
    def reverse(self): pass
    def sort(self): pass

    def find(self, tag):
        for c in self:
            if c.match(tag):
                return c
        return None

    def getTags(self):
        tags = []
        for c in self:
            tags.append(c.tag)
        return tags

    def makeUnique(self):
        tags = []
        kill = []
        for c in self:
            if c.tag not in tags:
                tags.append(c.tag)
            else:
                kill.append(c)
        for c in kill:
            del self.__data__[self.index(c)]

    def getSequence(self, unknown='???'):
        mapping = {}
        for coupling in self:
            shift1, shift2 = coupling.getType().shift
            if not mapping.has_key(int(coupling) + shift1):
                mapping[int(coupling) + shift1] = coupling.residues[0]
            if not mapping.has_key(int(coupling) + shift2):
                mapping[int(coupling) + shift2] = coupling.residues[1]
        keys = mapping.keys()
        if unknown != None:
            missing = 0
            for i in range(min(keys) + 1, max(keys)):
                if not mapping.has_key(i):
                    mapping[i] = unknown
                    missing += 1
            if missing:
                loki.info('Missing %d residue(s) in the sequence.' % missing)
        sequence = []
        for i in range(min(keys), max(keys) + 1):
            if mapping.has_key(i):
                sequence.append(mapping[i])
        return sequence


from Mara.Physics import gyromagnetic_constants as __gyros__

CouplingConstants = {
    'HnN' : __gyros__['H'] * __gyros__['N'] / 1.02**3,
    'CaCo' : __gyros__['C']**2 / 1.51**3
}
del __gyros__

