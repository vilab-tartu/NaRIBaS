#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.Atom (v. 0.6):
    Quantitise the universe.
    
    Methods and classes to create physical objects (atoms) and to group them 
    into larger constructs. A brief explanation of each class follows:

    Nucleus -- as the name implies, the nucleus of a atom, i.e. the invariant
               information of one atom type
    Atom -- an actual, physical atom w/ position, charge etc.
    Molecule -- a group of Atoms
    AminoAcid -- a Molecule that has extra amino acid specific information
    Residue -- an AminoAcid w/ sequence id (i.e. is part of a Polypeptide)
    Polypeptide -- a group of Residues
    Ensemble -- a group of Atoms, Molecules, Polypeptides or whatever

    Requirements: Python 2.4->
                  numpy
    
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

__all__ = ['Nucleus', 'Atom', 'PseudoAtom', 'Molecule', 'AminoAcid', \
        'Residue', 'Polypeptide', 'Ensemble', 'allNuclei', 'allAminoAcids']

from Mara.aux import *
from Mara.NMR import DCouplings
from Mara.Library import omniTranslator
#from Mara.Library.Natura import libChemistry
from Mara import loki
from numpy import dot, array, sum, float, zeros, sqrt
from copy import deepcopy

__atomic_mass__ = 1.66053886e-27
__electron_mass__ = 9.10938188e-31 / __atomic_mass__

class Nucleus:
    def __init__(self, name, protons, neutrons):
        self.name = str(name)
        self.protons = int(protons)
        self.neutrons = int(neutrons)
        self.mass = self.protons + self.neutrons
    def __str__(self):
        return self.name
    def __repr__(self):
        return "Nucleus('%s', %d, %d)"%(self.name, self.protons, \
                self.neutrons)
    def __int__(self):
        return self.protons

class Atom:
    def __init__(self, nucleus, name=None, charge=0, position=None):
        if isinstance(nucleus, Nucleus): self.nucleus = nucleus
        else: raise ArgumentError, 'Nucleus needs to be a Nucleus instance.'
        self.charge = int(charge)
        self.position = [float(x) for x in position[:3]]
        if name: self.name = str(name)
        else: self.name = str(self.nucleus)
        self.mass = self.nucleus.mass + (self.nucleus.protons - self.charge) \
                * __electron_mass__
    def __str__(self):
        return self.name
    def __repr__(self):
        return "Atom(%s, name='%s', charge=%d, position=%s)"%(
                repr(self.nucleus), self.name, self.charge, \
                repr(self.position))
    def __int__(self):
        return self.charge
    def __float__(self):
        return self.mass
    def rotate(self, matrix):
        try:
            self.position = list(dot(array(matrix), array(self.position)))
        except TypeError:
            raise ArgumentError, 'Rotation matrix must be a numpy array of shape (3,3).'
    def translate(self, vector):
        try:
            self.position = [x + float(y) for x,y in \
                    zip(self.position, vector[:3])]
        except TypeError:
            raise ArgumentError, 'Translation vector must a subscriptable object of (at least) three numerical (or convertable) components.'

class PseudoAtom(Atom):
    def __init__(self, name=None, position=None):
        self.nucleus = Nucleus('Q', 0, 0)
        self.charge = 0
        self.position = [float(x) for x in position[:3]]
        if name: self.name = str(name)
        else: self.name = str(self.nucleus)
        self.mass = 0.0

class Molecule:
    def __init__(self, name, atoms=[], nicknames=[]):
        self.name = str(name)
        self.nicknames = [str(x) for x in nicknames]
        self.atoms = []
        self.extend(atoms)
    def __str__(self):
        return self.name
    def __repr__(self):
        return "Molecule('%s', %s, %s)"%(self.name, self.atoms, self.nicknames)
    def __len__(self):
        return len(self.atoms)
    def __getitem__(self, i):
        return self.atoms[i]
    def __iter__(self):
        return BasicIterator(self.atoms)
    def append(self, atom):
        if isinstance(atom, Atom): self.atoms.append(atom)
        else: print "Discarding non-Atom '%s'."%repr(atom)
    def extend(self, atoms):
        for atom in atoms: self.append(atom)
    def rotate(self, matrix):
        for atom in self.atoms: atom.rotate(matrix)
    def translate(self, vector):
        for atom in self.atoms: atom.translate(vector)
    def getAtoms(self):
        return self.atoms
    def getBond(self, source, target):
        src, tgt = None, None
        for atom in self.atoms:
            if str(atom) == source:
                src = atom
            elif str(atom) == target:
                tgt = atom
        bond = None
        if None not in [src, tgt]:
            bond = Bond((src, tgt))
        return bond


class AminoAcid(Molecule):
    def __init__(self, name, abbrev='', symbol='', atoms=[]):
        self.name = str(name).capitalize()
        self.abbrev = str(abbrev).capitalize()
        self.symbol = str(symbol).capitalize()
        self.atoms = []
        self.extend(atoms)
    def __str__(self):
        return self.abbrev
    def __repr__(self):
        return "AminoAcid('%s', '%s', '%s', %s)"%(self.name, self.abbrev, \
                self.symbol, repr(self.atoms))
    def new(self, atoms):
        return AminoAcid(self.name, self.abbrev, self.symbol, atoms)

class Residue(Molecule):
    def __init__(self, name, id, abbrev='', symbol='', atoms=[], chain=None):
        if name in allAminoAcids:
            aa = allAminoAcids[name]
            self.name = aa.name
            self.abbrev = aa.abbrev
            self.symbol = aa.symbol
            self.atoms = aa.atoms
        else:
            self.name = name
            self.abbrev = abbrev
            self.symbol = symbol
            self.atoms = []
        if atoms and type(atoms) in [list, tuple]:
            self.atoms = list(atoms)
        self.id = int(id)
        self.chain = chain
    def __str__(self):
        return self.abbrev
    def __repr__(self):
        return "Residue('%s', %d, atoms=%s)"%(self.name, self.id, \
                repr(self.atoms))
    def __int__(self):
        return self.id
    def __float__(self):
        return float(self.id)
    def new(self, atoms):
        return Residue(self.name, self.id, self.abbrev, self.symbol, atoms)

class Polypeptide:
    def __init__(self, name, residues=[]):
        self.name = str(name)
        self.residues = []
        self.extend(residues)
    def __str__(self):
        return self.name
    def __repr__(self):
        return "Polypeptide('%s', %s)"%(self.name, repr(self.residues))
    def __len__(self):
        return len(self.residues)
    def __getitem__(self, i):
        return self.residues[i]
    def __iter__(self):
        return BasicIterator(self.residues)
    def append(self, residue):
        if isinstance(residue, Residue): self.residues.append(residue)
        else: loki.warn("Discarding non-Residue '%s'.", repr(residue))
    def extend(self, residues):
        for residue in residues: self.append(residue)
    def rotate(self, matrix):
        for residue in self.residues: residue.rotate(matrix)
    def translate(self, vector):
        for residue in self.residues: residue.translate(vector)
    def getDCouplingVectors(self, dcouplings, include='all', exclude=[]):
        if include == 'all':
            inc = range(1, len(self) + 1)
        else:
            inc = []
            for i in include:
                if type(i) is int:
                    inc.append(i)
                else:
                    loki.warn("discarding invalid residue '%s'", repr(i))
            inc = intersection(inc, range(1, len(self) + 1))
        exc = []
        for i in exclude:
            if type(i) is int:
                exc.append(i)
            else:
                loki.warn("discarding invalid residue '%s'", repr(i))
        include = sublist(inc, exc)
        loki.debug('include=' + repr(include))
        dcouplings = deepcopy(dcouplings)
        loki.debug('dcouplings=' + repr(dcouplings)) # REMOVE
        record = {}
        for i in range(len(dcouplings)):
            if type(dcouplings[i]) == str:
                dcouplings[i] = DCouplings[dcouplings[i]]
            dc = dcouplings[i]
            record.setdefault(dc.atoms[0], {})
            record.setdefault(dc.atoms[1], {})
        if 'HN' in record.keys():
            record.setdefault('H', {})
        wanted = record.keys()
        for i in include:
            for atom in self[i-1]:
                if str(atom) in wanted:
                    record[str(atom)][int(self[i-1])] = atom.position
        if 'H' in wanted and len(record['H'].keys()):
            for k in record['H'].keys():
                if not record['HN'].has_key(k):
                    record['HN'][k] = record['H'][k]
                    loki.info('using H for HN in residue %d.', k)
        loki.debug('record=%s', repr(record))
        table = {}
        for dc in dcouplings:
            x, y = dc.atoms
            residues = intersection(record[x].keys(), record[y].keys())
            for i in residues:
                tag = (i, str(dc))
                table[tag] = [m-n for m,n in zip(record[x][i], record[y][i])]
        tags = table.keys()
        tags.sort()
        vectors = [table[x] for x in tags]
        return tags, vectors

    def getAtoms(self):
        atoms = []
        for residue in self.residues:
            atoms.extend(residue.getAtoms())
        return atoms

    def getBonds(self, source, target):
        bonds = Bonds()
        for residue in self.residues:
            bond = residue.getBond(source, target)
            if bond is not None:
                bonds.append(bond)
        return bonds

    def getSequence(self, format='3-letter-aminoacids'):
        if not omniTranslator.recognise_language(format):
            raise ValueError, "Unknown format '%s'." % repr(format)
        sequence = []
        for residue in self:
            name = omniTranslator(str(residue), format)
            sequence.append(name)
        return sequence
    def getChainIDs(self):
        chains = []
        for residue in self.residues:
            if residue.chain not in chains:
                chains.append(residue.chain)
        if len(chains) == 1 and None in chains:
            return []
        else:
            return chains


class Ensemble:
    def __init__(self, name, members=[]):
        self.name = str(name)
        self.members = []
        self.extend(members)
    def __str__(self):
        return self.name
    def __repr__(self):
        return "Ensemble('%s', %s)"%(self.name, repr(self.members))
    def __len__(self):
        return len(self.members)
    def __getitem__(self, i):
        return self.members[i]
    def __delitem__(self, i):
        del self.members[i]
    def __iter__(self):
        return BasicIterator(self.members)
    def append(self, member):
        if isinstance(member, Atom) or isinstance(member, Molecule) or \
                isinstance(member, AminoAcid) or isinstance(member, Residue) \
                or isinstance(member, Polypeptide):
            self.members.append(member)
        else: loki.warn("Discarding unrecognisable '%s'.", repr(member))
    def extend(self, members):
        for member in members: self.append(member)
    def rotate(self, matrix):
        for member in self.members: member.rotate(matrix)
    def translate(self, vector):
        for member in self.members: member.translate(vector)
    def remove(self, name):
        kill = []
        for i in range(len(self)):
            if str(self[i]) == name:
                kill.append(i)
        while kill:
            del self[kill.pop()]
    def getAtoms(self):
        atoms = []
        for member in self.members:
            if not isinstance(member, Atom):
                atoms.extend(member.getAtoms())
            else:
                atoms.append(member)
        return atoms
    def getBonds(self, source, target):
        bonds = Bonds()
        for member in self.members:
            bonds.extend(member.getBonds(source, target))
        return bonds


class Bond:
    def __init__(self, atoms, vector=None, type='pseudo'):
        self.setAtoms(atoms)
        self.setVector(vector)
        self.setType(type)
    def __str__(self):
        return str(atoms[0]) + str(atoms[1])
    def __repr__(self):
        return 'Bond(%s, vector=%s, type=%s)' % (repr(self.__atoms__), \
                repr(self.__vector__), repr(self.__type__))
    def setAtoms(self, atoms):
        if atoms is not None \
        and (isinstance(atoms, list) or isinstance(atoms, tuple)) \
        and len(atoms) == 2:
            self.__atoms__ = (atoms[0], atoms[1])
        else:
            self.__atoms__ = None
    def setVector(self, vector=None):
        self.__vector__ = None
        if vector is None:
            try:
                self.__vector__ = list(\
                        array(self.__atoms__[1].position) - \
                        array(self.__atoms__[0].position)\
                        )
            except AttributeError, msg:
                loki.debug('Bond.setVector(): ' + repr(msg))
        else:
            if (isinstance(vector, list) or isinstance(vector, tuple)) \
            and len(vector) == 3:
                self.__vector__ = list(vector)
    def setType(self, type):
        if type in ['pseudo']: # TODO: add more types
            self.__type__ = type
        else:
            self.__type__ = 'pseudo'
    def getAtoms(self):
        return self.__atoms__
    def getVector(self):
        return self.__vector__
    def getType(self):
        return self.__type__

class Bonds:
    def __init__(self, bonds=None):
        self.set(bonds)
    def __str__(self):
        return str(self.__bonds__)
    def __repr__(self):
        return 'Bonds(%s)' % str(self.__bonds__)
    def __len__(self):
        return len(self.__bonds__)
    def __getitem__(self, i):
        return self.__bonds__[i]
    def __delitem__(self, i):
        del self.__bonds__[i]
    def __iter__(self):
        return BasicIterator(self.__bonds__)
    def set(self, bonds):
        self.reset()
        if bonds is not None:
            self.extend(bonds)
    def reset(self):
        self.__bonds__ = []
    def append(self, bond):
        if isinstance(bond, Bond):
            self.__bonds__.append(bond)
        else:
            try:
                self.__bonds__.append(Bond(bond))
            except:
                loki.warn("Invalid bond '%s'. Ignoring" % repr(bond))
    def extend(self, bonds):
        if isinstance(bonds, Bonds) or isinstance(bonds, list) \
        or isinstance(bonds, tuple):
            for bond in bonds:
                self.append(bond)
        else:
            loki.warn("Unknown format '%s' for a bond list. Ignoring." % \
                    type(bonds))
    def getAverage(self, normalise=False):
        r = zeros((3), float)
        for bond in self.__bonds__:
            v = bond.getVector()
            if v:
                loki.debug('v=' + repr(v))
                r += array(v)
        if normalise:
            r /= sqrt(sum(r**2))
        else:
            r /= float(len(self))
        return r


allNuclei = {
    'H': Nucleus('H',1,0),
    'He': Nucleus('He',2,2),
    'B': Nucleus('B',5,5),
    'C': Nucleus('C',6,6),
    'N': Nucleus('N',7,7),
    'O': Nucleus('O',8,8),
    'F': Nucleus('F',9,9),
    'P': Nucleus('P',15,15),
    'S': Nucleus('S',16,16),
    'SE': Nucleus('SE',34,34),
    'Q': Nucleus('Q',0,0)
}

alanine = AminoAcid('Alanine','Ala','A')
arginine = AminoAcid('Arginine','Arg','R')
asparagine = AminoAcid('Asparagine','Asn','N')
aspartate = AminoAcid('Aspartate','Asp','D')
cysteine = AminoAcid('Cysteine','Cys','C')
glutamate = AminoAcid('Glutamate','Glu','E')
glutamine = AminoAcid('Glutamine','Gln','Q')
glycine = AminoAcid('Glycine','Gly','G')
histidine = AminoAcid('Histidine', 'His','H')
isoleucine = AminoAcid('Isoleucine','Ile','I')
leucine = AminoAcid('Leucine','Leu','L')
lysine = AminoAcid('Lysine','Lys','K')
methionine = AminoAcid('Methionine','Met','M')
phenylalanine = AminoAcid('Phenylalanine','Phe','F')
proline = AminoAcid('Proline','Pro','P')
serine = AminoAcid('Serine','Ser','S')
threonine = AminoAcid('Threonine','Thr','T')
tryptophan = AminoAcid('Tryptophan','Trp','W')
tyrosine = AminoAcid('Tyrosine','Tyr','Y')
valine = AminoAcid('Valine','Val','V')
selenomethionine = AminoAcid('Selenomethionine','Mse','Ms')

def create_allAminoAcids():
    all = [alanine, arginine, asparagine, aspartate, cysteine, glutamate, \
            glutamine, glycine, histidine, isoleucine, leucine, lysine, \
            methionine, phenylalanine, proline, serine, threonine, tryptophan, \
            tyrosine, valine, selenomethionine]
    allAminoAcids = {}
    for aa in all:
        for key in [aa.name, aa.abbrev, aa.symbol]:
            allAminoAcids[key.lower()] = aa
            allAminoAcids[key.upper()] = aa
            allAminoAcids[key.capitalize()] = aa
    return allAminoAcids

allAminoAcids = create_allAminoAcids()

