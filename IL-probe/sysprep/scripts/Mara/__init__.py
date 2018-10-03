#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Main constructor for Mara's Accruement of Relevant Algorithms 
"""
__version__ = '0.9.1'
__all__ = ['aux', 'Atom', 'Biology', 'Interface', 'IO', 'Library', \
        'Linguistics', 'Logic', 'Math', 'mpi', 'NMR', 'Physics', 'svd']

import sys
if sys.version < '2.3':
    from Mara.aux import OldLogger
    loki = OldLogger()
else:
    import logging
    logging.basicConfig(level=logging.WARNING, format=\
            '%(asctime)s %(levelname)-8s [%(name)s.%(module)s] : %(message)s',
            datefmt='%H:%M:%S')
    loki = logging.getLogger('Mara')

