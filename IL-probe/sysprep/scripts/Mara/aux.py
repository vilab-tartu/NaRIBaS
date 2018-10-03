#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.aux (v. 1.2):
    Auxiliary functions for all kinds of purposes.
    
    Requirements: Python 2.2->

    Author: Martti Louhivuori (martti.louhivuori@helsinki.fi), 30.12.2004

    Date: 20.10.2010

    ---

    Copyright (C) 2004-2010  Martti Louhivuori

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    The full license text can be found in the file LICENSE.
"""
from Mara import loki
import re

def convert_seconds(seconds):
    """
    Converts seconds to days, hours, minutes and seconds as appropriate.
    """
    days = 0
    hours = 0
    minutes = 0
    while seconds >= 86400:
        seconds -= 86400
        days += 1
    while seconds >= 3600:
        seconds -= 3600
        hours += 1
    while seconds >= 60:
        seconds -= 60
        minutes += 1
    if days:
        s = '%d days'%days
    else:
        s = ''
    if hours:
        s += ' %d h'%hours
    if minutes:
        s += ' %d min'%minutes
    if seconds:
        s += ' %.1f s'%seconds
    return s.strip()

class BasicIterator:
    """
    A iterator built from scratch. (obsolete ?)
    """
    def __init__(self, x):
        self.x = x
        self.max = len(x)-1
        self.i = -1
        self.stop = False
    def __iter__(self):
        return self
    def next(self):
        if not self.stop and self.i < self.max:
            self.i += 1
            return self.x[self.i]
        else:
            self.stop = True
            raise StopIteration

def quicksort(list, start=0, end=None):
    """
    Sort a list in ascending alphabetical order using the quicksort method.
    Originally written by Magnus Lie Hetland.
    """
    if end == None:
        end = len(list)-1
    if start < end:                             # If there are two or more elements...
        pivot = list[end]                       # Partition around the last value
        bottom = start-1                        # Start outside the area to be partitioned
        top = end                               # Ditto
        done = 0
        while not done:                         # Until all elements are partitioned...
            while not done:                     # Until we find an out of place element...
                bottom = bottom+1               # ... move the bottom up.
                if bottom == top:               # If we hit the top...
                    done = 1                    # ... we are done.
                    break
                if list[bottom] > pivot:        # Is the bottom out of place?
                    list[top] = list[bottom]    # Then put it at the top...
                    break                       # ... and start searching from the top.
            while not done:                     # Until we find an out of place element...
                top = top-1                     # ... move the top down.
                if top == bottom:               # If we hit the bottom...
                    done = 1                    # ... we are done.
                    break
                if list[top] < pivot:           # Is the top out of place?
                    list[bottom] = list[top]    # Then put it at the bottom...
                    break                       # ...and start searching from the bottom.
        list[top] = pivot                       # Put the pivot in its place.
        quicksort(list, start, top-1)           # ... and sort both halves.
        quicksort(list, top+1, end)


def vararg_callback(option, opt_str, value, parser):
    """
    A callback function for optparse that handles argument lists of
    variable length. NOTE: Won't handle negative numbers correctly!
    Author: Guido van Rossum, Python Library Reference
    fixed: option.dest -> opt_str[2:] (ML)
    improved: swap '-' -> '_' in option name (ML)
    """
    assert value is None
    done = 0
    value = []
    rargs = parser.rargs
    while rargs:
        arg = rargs[0]
        # Stop if we hit an arg like "--foo", "-a", "-fx", "--file=f",
        # etc.  Note that this also stops on "-3" or "-3.0", so if
        # your option takes numeric values, you will need to handle
        # this.
        if ((arg[:2] == "--" and len(arg) > 2) or \
                (arg[:1] == "-" and len(arg) > 1 and arg[1] != "-")):
            break
        elif arg == '-':
            del rargs[0]
            break
        else:
            value.append(arg)
            del rargs[0]
    setattr(parser.values, opt_str[2:].replace('-', '_'), value)

def intersection(x, y):
    """
    Calculate the intersection of two lists.

    NOTE: non-unique elements will be lost
    """
    i = []
    for e in x:
        if e in y:
            i.append(e)
    return i

def union(x, y):
    """
    Calculate the union of two lists.

    NOTE: non-unique elements will be lost
    """
    u = []
    for e in x:
        u.append(e)
    for e in y:
        if e not in u:
            u.append(e)
    return u

def complement(x, y):
    """
    Calculate the complement of two lists, i.e. those elements in the first 
    list that are not present in the second list.
    """
    sx = []
    for e in x:
        if e not in y:
            sx.append(e)
    return sx
sublist = complement

def symmetric_difference(x, y):
    """
    Calculate the symmetric difference of two lists, i.e. XOR.
    """
    return complement(union(x, y), intersection(x, y))
xor = symmetric_difference


def pruneTags(x, y):
    """
    Return the intersection of two 'tag' lists.
    """
    if type(x) == tuple and len(x) == 2: tx, vx = x
    else: tx, vx = x.getTags(), None
    if type(y) == tuple and len(y) == 2: ty, vy = y
    else: ty, vy = y.getTags(), None
    tags = intersection(tx, ty)
    kx, ky = sublist(tx, tags), sublist(ty, tags)
    if vx == None:
        for k in kx: x.remove(k)
    else:
        for k in kx: 
            del vx[tx.index(k)]
            tx.remove(k)
        x = (tx, vx)
    if vy == None:
        for k in ky: y.remove(k)
    else:
        for k in ky: 
            del vy[ty.index(k)]
            ty.remove(k)
        y = (ty, vy)
    return x, y

def parse_option_hash(hash, translator=None, default_keys=None, \
        default_value=None, default_function=None, translator_target='mara'):
    """
    Parse and fix/append an option hash.

    Arguments:
        hash             -- an unformatted option hash or a default value
        translator       -- a translator to format hash keys
                            (default: as is)
        default_keys     -- a list of keys to include
                            (default: keys of the original hash)
        default_value    -- a default value to use in the absense of a
                            proper input hash
        default_function -- a default function to apply to values
                            (default: as is)
    Returns:
        out              -- a formatted option hash
    """
    loki.debug('hash(in)=%s' % repr(hash))
    # a do-nothing translator will be used if none given
    if not translator:
        translator = lambda x,y: x
    # a do-nothing function will be used if none given
    if not default_function:
        default_function = lambda x: x
    # translate each key & apply default function to the value
    if isinstance(hash, dict):
        out = {}
        for key in hash.keys():
            if translator_target is not None:
                k = translator(key, translator_target)
            else:
                k = translator(key)
            out[k] = default_function(hash[key])
    # or create a hash from the default keys & the given value
    elif hash in [int, float, str]:
        out = {}.fromkeys(default_keys, default_function(hash))
    # or create a hash from the default keys & the default value
    else:
        loki.info('Incomprehensible hash -> using %s' % repr(default_value))
        out = {}.fromkeys(default_keys, default_value)
    loki.debug('hash(out)=%s' % repr(out))
    return out

class ArgumentError(Exception):
    """
    A custom error. (obsolete?)
    """
    def __init__(self, txt, name='?'):
        self.txt = txt
        self.name = name
    def __str__(self):
        return '[%s] %s'%(self.name, self.txt)

class SpongeStream:
    """
    A stream that simply absorbs anything thrown at it.
    """
    def __init__(self, filter=True):
        self.__filter__ = bool(filter)
        self.__data__ = ''
    def __iter__(self):
        return BasicIterator(self.__data__)
    def write(self, x):
        self.__data__ += str(x)
    def read(self):
        return self.__data__
    def readlines(self):
        return filter(None, \
                [x.strip() for x in self.__data__.strip().split('\n')])

class nameGenerator:
    """
    Generate (file)names either as a prefix to the argument of call or as 
    a consequent list of numbered names fulfilling the python string given 
    to the creator.

    Inputs (to creator):
    output -- name prefix / a python string that accepts an integer as an 
              argument
    suffix -- swap the ending (after the last dot) of the alternative name 
              to this or just add it to the end of the primary name 
              (default: 'alt')
    regex  -- a regular expression that must be found for the custom naming
              scheme to come into play, i.e. it should point to the python
              placeholder that will receive the integer argument
              (default: '\%[.0-9]*d')
    single -- return only the primary name
              (default: False)
    swap   -- return only the alternative name
              (default: False)

    Arguments (to a call):
    input -- original name to be prefixed or nothing at all

    Output(s) (from a call):
    name   -- prefixed / custom filename
    (alt)  -- similar filename w/ an alternative ending (only if not single)
    """
    def __init__(self, output, suffix='alt', regex='\%[.0-9]*d', \
            single=False, swap=False):
        self.__output__ = str(output)
        self.__suffix__ = str(suffix)
        if re.search(str(regex), self.__output__):
            self.__custom__ = True
            self.__index__ = 0
        else:
            self.__custom__ = False
            self.__index__ = None
        self.setSingle(single)
        self.setSwap(swap)
    def __call__(self, input=None):
        if self.__custom__:
            self.__index__ += 1
            name = self.__output__ % self.__index__
        elif input:
            name = self.__output__ + str(input)
        else:
            raise TypeError, 'Only in custom mode are no arguments needed.'
        if self.__single__ and not self.__swap__:
            return name
        if name.count('.'):
            alt = '.'.join(name.split('.')[:-1]) + '.' + self.__suffix__
        else:
            alt = name + '.' + self.__suffix__
        if self.__swap__:
            return alt
        else:
            return name, alt
    def rewind(self):
        """
        Reset the internal counter.
        """
        if self.__index__: self.__index__ = 0
    def isCustom(self):
        return self.__custom__
    def setSingle(self, x=True):
        self.__single__ = bool(x)
    def setSwap(self, x=True):
        self.__swap__ = bool(x)

class pseudoGenerator:
    """
    A pseudo nameGenerator that smells and tastes like the real one, but 
    returns only what it is explicitly told to at creation.

    Inputs (to creator):
    first  -- the first filename list to iterate through
              (default: [])
    second -- the second filename list to iterate through
              (default: [])
    single -- return only the primary name
              (default: False)
    swap   -- return only the alternative name
              (default: False)

    Arguments (to a call):
    x -- a dummy placeholder

    Outputs (from a call):
    first, second -- corresponding filenames or None if out of bounds
    """
    def __init__(self, first=[], second=[], single=False, swap=False):
        self.__first__ = first
        self.__second__ = second
        self.__index__ = 0
        self.setSingle(single)
        self.setSwap(swap)
    def __call__(self, x):
        self.__index__ += 1
        if self.__index__ < len(self.__first__):
            first = self.__first__[self.__index__]
        else:
            first = None
        if self.__single__:
            return first
        if self.__index__ < len(self.__second__):
            second = self.__second__[self.__index__]
        else:
            second = None
        if self.__swap__:
            return second
        else:
            return first, second
    def rewind(self):
        """
        Reset the internal counter.
        """
        self.__index__ = 0
    def setSingle(self, x=True):
        self.__single__ = bool(x)
    def setSwap(self, x=True):
        self.__swap__ = bool(x)

class OldLogger:
    """
    A do-it-yourself logger. (obsolete?)
    """
    def __init__(self, verbose=False, debug=False):
        self.verbose = bool(verbose)
        self.debug = bool(debug)
    def debug(self, txt):
        if self.debug: print 'DEBUG:',txt
    def info(self, txt):
        if self.verbose or self.debug: print 'INFO:',txt
    def warn(self, txt):
        print 'WARNING:',txt
    def error(self, txt):
        print 'ERROR:',txt
    def critical(self, txt):
        print 'CRITICAL:',txt

def compare_none_last(x, y):
    """
    Compare x and y and return -1, 0 or 1 on x<y, x=y and x>y. None is 
    treated to be more than anything, i.e. unlike the default behaviour in
    built-in sort().
    """
    if x is None:
        if y is None: return -1
        else: return 1
    elif y is None: return -1
    elif x < y: return -1
    elif x > y: return 1
    else: return 0

def sort_none_last(x):
    """
    Sort a list in an ascending order leaving Nones and (None, x) tuples to
    the end.
    """
    x.sort(compare_none_last)
    hold = []
    for i in x:
        if type(i) is tuple and i[0] is None:
            hold.append(i)
    for i in hold:
        x.remove(i)
    for i in hold:
        x.append(i)
    return x

def inverse_dict(x):
    """
    Invert a dictionary, i.e. from key:value to value:key.
    """
    new = {}
    for key, value in x.items():
        new[value] = key
    return new

def sort_dict_values(x, reverse=False):
    """
    Sort a dictionary by its values.
    """
    if reverse:
        d = sorted(x.items(), lambda a,b: cmp(b[1], a[1]))
    else:
        d = sorted(x.items(), lambda a,b: cmp(a[1], b[1]))
    return d

def escape_parentheses(x):
    """
    Escape parentheses by adding a \ in front of them.
    """
    for p in ['(', ')', '[', ']', '{', '}']:
        glue = '\\' + p
        x = glue.join(x.split(p))
    return x


class AtomBin:
    def __init__(self):
        self.template = '%6d %10s %6d %6s %6s %6d %10.3f %10.4   ; %s'
        self.formats = {
                'nr':'%6d', 'type':'%10s', 'resnr':'%6d', 
                'residue':'%6s', 'atom':'%6s', 'cgnr':'%6d', 
                'charge':'%10.3f', 'mass':'%10.4f', 'typeB':'%10s', 
                'chargeB':'%10.3f', 'massB':'%10.4f'
                }
        self.atoms = []
        self.residues = {}
    def add(self, atom):
        if atom is str:
            atom = self.parse_atom(atom)
        self.residues.setdefault(atom['resnr'], [])
        self.residues[atom['resnr']].append(len(self.atoms))
        self.atoms.append(atom)
    def get_residue(self, resid):
        atoms = []
        for i in self.residues[resid]:
            atoms.append(self.atoms[i])
        return atoms
    def parse_atom(line):
        atom = {}
        line = line.strip()
        if line[0] == ';':
            return None
        if line.count(';'):
            line, comment = line.split(';')
            line = line.strip()
            atom['comment'] = comment.strip()
        parts = line.split()
        logging.debug('parts=' + repr(parts))
        if len(parts) < 6:
            logging.warn('Not enough fields in atom definition. Ignoring.')
            return None
        atom['nr'] = int(parts.pop(0))
        atom['type'] = parts.pop(0)
        atom['resnr'] = int(parts.pop(0))
        atom['residue'] = parts.pop(0)
        atom['atom'] = parts.pop(0)
        atom['cgnr'] = int(parts.pop(0))
        if len(parts): atom['charge'] = float(parts.pop(0))
        if len(parts): atom['mass'] = float(parts.pop(0))
        if len(parts): atom['typeB'] = parts.pop(0)
        if len(parts): atom['chargeB'] = float(parts.pop(0))
        if len(parts): atom['massB'] = float(parts.pop(0))
        return atom
    def format_atom(atom):
        elements = []
        for key in ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr']:
            if atom.has_key(key):
                elements.append(formats[key] % atom[key])
            else:
                raise KeyError, "required element '%s' is missing" % key
        for key in ['charge', 'mass', 'typeB', 'chargeB', 'massB']:
            if atom.has_key(key):
                elements.append(formats[key] % atom[key])
            else:
                break
        if atom.has_key('comment'):
            elements.append('; %s' % atom['comment'])
        return ' '.join(elements)
