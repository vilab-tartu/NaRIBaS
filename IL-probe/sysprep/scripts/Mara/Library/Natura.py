#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.Library.Natura (v. 0.8):
    A collection of libraries of physical entities etc.

    Requirements: Python 2.2->

    Author: Martti Louhivuori (martti.louhivuori@helsinki.fi)

    Date: 8.2.2006

    ---

    Copyright (C) 2006  Martti Louhivuori

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    The full license text can be found in the file LICENSE.
"""

from Mara.Atom import Nucleus, AminoAcid
from Mara.Linguistics import *
from Mara.Library import CommonLib, ChemLib

"""
-----------------------------------------------------------------------------
Chemistry section. Atoms, amino acids etc.
-----------------------------------------------------------------------------
"""

libChemistry = ChemLib('libChemistry: Chemical library')

libChemistry.nuclei = {
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
    'X': Nucleus('X',0,0)
}

# FIX THIS!!! Should have routines for loading atom coordinates from a file.
libChemistry.aminoacids = {} 


