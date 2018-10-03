#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.Biology (v. 2.1):
    Biologically involved general algorithms.
    
    Requirements: Python 2.2->

    Author: Martti Louhivuori (martti.louhivuori@helsinki.fi)

    Date: 13.2.2006

    ---

    Copyright (C) 2006  Martti Louhivuori

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    The full license text can be found in the file LICENSE.
"""
from Mara import loki
from Mara.Library import omniTranslator, libLingua
import os, sys

def parse_sequence(args):
    """
    Parse an amino acid sequence.

    Arguments:
        args     -- a list of sequence items or a name of a file containing
                    them, e.g. 'GLU PRO GLU CYS' or 'EPEC GLK C EEK'
    Returns:
        sequence -- a list of 3-letter amino acid symbols
    """
    loki.debug('parse_sequence < %s' % repr(args))
    if isinstance(args, str) and os.path.isfile(args):
        fname = args
    elif len(args) == 1 and isinstance(args[0], str):
        if os.path.isfile(args[0]):
            fname = args[0]
        else:
            if args[0].count(' '):
                args = args[0].split()
            else:
                args = args[0]
            fname = None
    else:
        fname = None
    if fname:
        f = open(fname)
        seq = f.read()
        f.close()
        loki.info("Read sequence from file '%s'." % fname)
        args = seq.strip().split()
        loki.debug('args=%s' % repr(args))
#        sequence = []
#        for aa in seq.strip().split():
#            try:
#                sequence.append(omniTranslator(aa.capitalize(), \
#                        '3-letter-aminoacids'))
#            except KeyError:
#                loki.warn("Discarding unknown aminoacid '%s'." % repr(aa))
#    else:
        # check whether all the sequence items are separated from each other
    args = [x.capitalize() for x in args]
    separated = True
    for a in args:
        if not (a in libLingua.dictionaryAmino1 or \
                a in libLingua.dictionaryAmino3):
            separated = False
    loki.debug('separated=%s' % repr(separated))
    sequence = []
    if separated:
        # append each item after converting it to a 3-letter symbol
        for a in args:
            try:
                sequence.append(omniTranslator(a.capitalize(), \
                        '3-letter-aminoacids'))
            except KeyError:
                loki.warn("Discarding unknown aminoacid '%s'." % repr(a))
    else:
        # jam all symbols together (hope they are all 1-letter symbols)
        aa = ''
        for a in args:
            aa += str(a)
        aa = aa.replace(' ', '')
        loki.debug('aa=%s' % repr(aa))
        # append each item after converting it to a 3-letter symbol
        for a in list(aa):
            try:
                sequence.append(omniTranslator(a, '3-letter-aminoacids'))
            except KeyError:
                loki.warn("Discarding unknown aminoacid '%s'." % repr(a))
    loki.debug('parse_sequence > %s' % repr(sequence))
    return sequence

def write_sequence(seq, output=None, format='long', use_translator=True):
    """
    Write amino acid sequence in the desired format either to a stream
    given or to STDOUT.
    
    TODO: - add some flexibility to the output formats
    """
    if output is None:
        output = sys.stdout
    if format == 'long':
        l = []
        while seq:
            s = seq.pop(0)
            if use_translator:
                s = omniTranslator(s.capitalize(), '3-letter-aminoacids')
            l.append(s.upper())
            if len(l) == 15:
                output.write(' '.join(l) + '\n')
                l = []
        output.write(' '.join(l) + '\n')
    elif format == 'full':
        l = []
        while seq:
            s = seq.pop(0)
            if use_translator:
                s = omniTranslator(s.capitalize(), 'english')
            l.append(s)
        output.write(' '.join(l) + '\n')
    else:
        l = []
        w = ''
        while seq:
            s = seq.pop(0)
            if use_translator:
                s = omniTranslator(s.capitalize(), '1-letter-aminoacids')
            w += s.upper()
            if len(w) == 10:
                l.append(w)
                w = ''
            if len(l) == 6:
                output.write(' '.join(l) + '\n')
                l = []
        if w:
            l.append(w)
        output.write(' '.join(l) + '\n')


