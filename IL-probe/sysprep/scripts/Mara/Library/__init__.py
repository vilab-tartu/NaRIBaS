#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.Library (v. 0.8):
    A collection of dictionaries, translation codeces etc.

    This file includes class declarations and initialisation of common 
    Translators (see below for more info).

    Requirements: Python 2.2->

    TODO: - add more specialised translators
          - find a way around the vicious-circle in imports

    Author: Martti Louhivuori (martti.louhivuori@helsinki.fi), 8.2.2006

    Date: 20.10.2010

    ---

    Copyright (C) 2006  Martti Louhivuori

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    The full license text can be found in the file LICENSE.
"""

from Mara.Linguistics import *

class CommonLib:
    """
    A very light common library infrastructure.
    """
    def __init__(self, description='lib???'):
        self.description = str(description)
    def __str__(self):
        return self.description

class ChemLib(CommonLib):
    """
    A specialised library for chemical entities.
    """
    def __init__(self, description='libChem???'):
        self.description = str(description)
        self.nuclei = {}
        self.aminoacids = {}
    def getNucleus(self, name):
        if name in self.nuclei.keys():
            return self.nuclei[name]
        else:
            raise ArgumentError, "Unknown nucleus '%s'."%str(name)
    def getAminoAcid(self, name):
        realname = omniTranslator(name, 'english')
        if realname in self.aminoacids.keys():
            return self.aminoacids[realname]
        else:
            raise ArgumentError, "Unknown amino acid '%s'."%str(name)

from Lingua import libLingua

# the all-knowing, albeit erring, super translator
omniTranslator = Translator(libLingua.allDictionaries, libLingua.allCodices)
# NMR specific translator, good when working with couplings
NMRTranslator = Translator(\
        [libLingua.dictionaryMara, libLingua.dictionaryNMR], \
        [libLingua.codexNMR], 'mara')
# atom-name translator w/ both atomistic & CG names
atomTranslator = Translator(
        [libLingua.dictionaryAtom, libLingua.dictionaryMartini], 
        [], 'atom')
# aminoacid translator
aaTranslator = Translator(
        [libLingua.dictionaryAmino1, libLingua.dictionaryAmino3, \
                libLingua.dictionaryEnglish], \
        [libLingua.codexAmino, libLingua.codexAminoSymbols], \
        '3-letter-aminoacids')

