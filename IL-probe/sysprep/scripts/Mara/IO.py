#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.IO (v. 0.9.2):
    Input / output interfaces for various data formats as well as a 
    low-memory file crawler.

    BlockIO
    EOF
    FastaIO
    FileCrawler
    GeneralFixed
    GnuplotLogIO
    GromosIO
    HoleIO
    ITPIO
    IndexIO
    PDB
    PalesIO
    PalesOutputIO
    XVGIO

    Requirements: Python 2.4->
                  BioPython (for PDB)

    TODO: - add virtual avatars and program version detection
          - fix rawPDB()
          - fix PDB.write() to include header
          - work around TER mishandling in Bio.PDB.PDBParser

    Author: Martti Louhivuori (m.j.louhivuori@rug.nl), 29.09.2005

    Date: 20.10.2010

    ---

    Copyright (C) 2005-2010  Martti Louhivuori

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    The full license text can be found in the file LICENSE.
"""

import re, os, sys, string
from Mara.Atom import *
#from Mara.Library import libLingua, omniTranslator
from Mara.Library import atomTranslator, omniTranslator
from Mara.Library.Lingua import libLingua
from Mara.aux import ArgumentError, SpongeStream, escape_parentheses, \
        complement, intersection
from Mara.NMR import DCouplings, Coupling, Couplings
from Mara import loki
from tempfile import TemporaryFile
#from Mara.Biology import aminolib
try:
    from Bio.PDB import PDBParser
    from Bio.PDB.PDBIO import PDBIO
    from Bio.PDB.StructureBuilder import StructureBuilder
    from Bio import PDB as bPDB
    __biopython_installed__ = True
except ImportError:
    __biopython_installed__ = False

### DEFAULT PARAMETERS TO USE ###
__default_bfactor__ = 0.0
__default_occupancy__ = 1.0
__default_altloc__ = ' '
__default_hetatm__ = ' '
__default_segid__ = ''
__default_icode__ = ' '
### DEFAULT PARAMETERS TO USE ###


class EOF(Exception):
    def __init__(self): pass

class FileCrawler:
    """
    Crawl through a file reading back and forth without loading 
    anything to memory.
    """
    def __init__(self, filename):
        """
        Create a FileCrawler instance.

        arguments:
          filename -- path to the file to crawl
        """
        self.__fp__ = open(filename)
        self.tell = self.__fp__.tell
        self.seek = self.__fp__.seek
    def readline(self):
        """
        Read a single line from the current position onward.
        """
        return self.__fp__.readline()
    def readchar(self):
        """
        Read a single character from the current position.
        """
        return self.__fp__.read(1)
    def read(self):
        """
        Read rest of the file.
        """
        return self.__fp__.read()
    def nextline(self):
        """
        Go to the next line.
        """
        c = self.readline()
    def prevline(self):
        """
        Go to the previous line.
        """
        found = 0
        try:
            self.prev()
        except EOF:
            self.rewind()
            found = 2
        while found < 2:
            c = self.readchar()
            if c == '\n':
                found += 1
            if found < 2:
                try:
                    self.prev(2)
                except EOF:
                    self.rewind()
                    break
    def next(self, i=1):
        """
        Go to the next character.
        """
        self.seek(i, 1)
    def prev(self, i=1):
        """
        Go to the previous character.
        """
        try:
            self.seek(-i, 1)
        except IOError:
            raise EOF
    def wind(self):
        """
        Go to the end of the file.
        """
        self.seek(0,2)
    def rewind(self):
        """
        Go to the beginning of the file.
        """
        self.seek(0)
    def get_line(self, i):
        """
        Read a single line.

        arguments:
          i -- the number of the line to read
        """
        self.rewind()
        j = 0
        while j < i:
            self.nextline()
        return self.readline()

class PDB:
    """
    Read / write coordinates in PDB format.
    """
    
    def __init__(self, filename=None, verbose=False, debug=False):
        """
        Create a PDB instance.

        arguments:
          filename -- PDB file to read (default: None)
          verbose  -- run in verbose mode (default: False) [defunc]
          debug    -- run in debug mode (default: False) [defunc]
        """
        self.head = None
        self.name = None
        self.author = []
        self.deposition = None
        self.release = None
        self.journal = None
        self.keywords = []
        self.resolution = None
        self.method = None
        self.id = None
        self.ensemble = None
        self.remarks = []
        self.models = None
        self.__parser__ = PDBParser()
        self.__writer__ = PDBIO()
        self.__filename__ = ''
        if filename != None and os.path.isfile(filename):
            self.read(filename)
    def __str__(self):
        return self.__filename__
    def __repr__(self):
        return 'PDB(%s)' % str(self)

    def read_remarks(self, filename):
        remarks = []
        ro = re.compile('^REMARK')
        for line in open(filename).readlines():
            if ro.match(line):
                remarks.append(line)
        return remarks

    def read(self, filename, discard=[], join_models=True):
        """
        Read a PDB file.

        arguments:
          filename    -- PDB file to read
          discard     -- names of atoms that should be ignored
                         (default: []) [defunc]
          join_models -- join all models into one ensemble (default: True)
        """
        if not os.path.isfile(filename):
            raise ArgumentError, "File '%s' not found." % filename
        else:
            self.__filename__ = filename
            loki.debug("Reading file '%s'." % str(self))
        self.id = os.path.split(filename)[1].split('.')[0]
        # catch output & errors into loki
        so = sys.stdout
        se = sys.stderr
        ostream = SpongeStream()
        estream = SpongeStream()
        sys.stdout = ostream
        sys.stderr = estream
        structure = self.__parser__.get_structure(self.id, filename)
        sys.stdout = so
        sys.stderr = se
        for line in ostream.readlines():
            loki.debug('Bio.PDB -> %s' % line)
        for line in estream.readlines():
            loki.warn('Bio.PDB -> %s' % line)
        # header
        h = structure.header
        try: self.head = h['head'].strip()
        except KeyError: pass
        try: self.name = h['name'].strip()
        except KeyError: pass
        try: self.deposition = h['deposition_date'].strip()
        except KeyError: pass
        try: self.release = h['release_date'].strip()
        except KeyError: pass
        try: self.resolution = h['resolution']
        except KeyError: pass
        try: self.method = h['structure_method'].strip()
        except KeyError: pass
        try: self.author = [x.strip() for x in h['author'].split(',')]
        except KeyError: pass
        try:
            if h['journal_reference']:
                self.journal = h['journal_reference'].strip()
            else:
                self.journal = h['journal'].strip()
        except KeyError: pass
        if h.has_key('keywords'):
            self.keywords = [x.strip() for x in h['keywords'].split(',')]
        else:
            self.keywords = []
        self.remarks = self.read_remarks(filename)
        # molecules
        models = structure.get_list()
        if not join_models: self.models = []
        polypeptides = []
        molecules = []
        for model in models:
            residues = []
            for chain in model.get_list():
                for residue in chain.get_list():
                    atoms = []
                    for atom in residue.get_list():
                        name = atomTranslator.base(atom.get_name())
                        if atomTranslator.language(name) == 'atom':
                            nucleus = allNuclei[name]
                        else:
                            nucleus = Nucleus('CG', 0, 72)
                        if name == None:
                            loki.warn("Unknown atom '%s'."%str(atom.get_name()))
                        else:
                            atoms.append(Atom( \
                                    nucleus, \
                                    atom.get_name(), \
                                    position=list(atom.get_coord())))
                            atoms[-1].bfactor = atom.get_bfactor()
                            atoms[-1].occupancy = atom.get_occupancy()
                            atoms[-1].altloc = atom.get_altloc()
                            atoms[-1].fullname = atom.get_fullname()
                    resname = residue.get_resname().capitalize()
                    if resname in allAminoAcids:
                        residues.append(Residue(resname, \
                                residue.get_id()[1], atoms=atoms, \
                                chain=chain.get_id()))
                        residues[-1].segid = residue.get_segid()
                    else:
                        if resname == 'Hoh':
                            name = 'water'
                            nicks = ['hoh']
                        else:
                            name = resname
                            nicks = []
                        molecules.append(Molecule(name, atoms, nicks))
                        molecules[-1].segid = residue.get_segid()
            if residues:
                polypeptides.append(Polypeptide('', residues))
            if not join_models and (polypeptides or molecules):
                self.models.append(Ensemble(model.get_id() + 1, \
                        polypeptides + molecules))
                polypeptides = []
                molecules = []
        if join_models:
            self.models = [Ensemble(self.id, polypeptides + molecules)]
        try:
            self.ensemble = self.models[0]
        except TypeError, IndexError:
            loki.warn("Not a single structure/model was found in '%s'." \
                    % filename)
        # trailer ?
        loki.debug("File '%s' read." % str(self))

    def __build_structure__(self):
        # catch output & errors into loki
#        so = sys.stdout
        se = sys.stderr
#        ostream = SpongeStream()
        estream = SpongeStream()
#        sys.stdout = ostream
        sys.stderr = estream
        builder = StructureBuilder()
        builder.init_structure(str(self.models[0]))
        for i in range(len(self.models)):
            builder.init_model(i+1)
            queue = []
            id = 0
            for molecule in self.models[i]:
                _chain_ids = list('ABCDEFGHIJKLMNOPQRSTUVW')
                if isinstance(molecule, Polypeptide):
                    chain_ids = molecule.getChainIDs()
                    if not chain_ids:
                        chain_ids = _chain_ids
                    try:
                        current_chain = chain_ids.pop(0)
                    except IndexError:
                        current_chain = 'X'
                    builder.init_chain(current_chain)
                    for residue in molecule:
                        if residue.chain and residue.chain != current_chain:
                            try:
                                current_chain = chain_ids.pop(0)
                            except IndexError:
                                current_chain = 'X'
                            builder.init_chain(current_chain)
                        try:
                            segid = residue.segid
                        except AttributeError:
                            segid = __default_segid__
                        # fix this to include hetres and placement
                        builder.init_seg(segid)
                        builder.init_residue(str(residue).upper(), \
                                __default_hetatm__, \
                                int(residue), __default_icode__)
                        for atom in residue:
                            try: bfactor = atom.bfactor
                            except AttributeError: bfactor = __default_bfactor__
                            try: occupancy = atom.occupancy
                            except AttributeError: occupancy = __default_occupancy__
                            try: altloc = atom.altloc
                            except AttributeError: altloc = __default_altloc__
                            try: fullname = atom.fullname
                            except AttributeError: fullname = \
                                    string.center(str(atom), 4)
                            builder.init_atom(str(atom), atom.position, \
                                    bfactor, occupancy, ' ', fullname)
                        id += 1
                else:
                    queue.append(molecule)
            id += 10
            if queue:
                builder.init_chain('Z')
                for molecule in queue:
                    try: segid = molecule.segid
                    except AttributeError: segid = __default_segid__
                    if len(str(molecule)) == 3:
                        resname = str(molecule).upper()
                    else:
                        found = False
                        for nick in molecule.nicknames:
                            if len(nick) == 3:
                                resname = nick.upper()
                                found = True
                                break
                        if not found: resname = str(molecule)[:3].upper()
                    if str(molecule) == 'water': hetero = 'W'
                    else: hetero = __default_hetatm__
                    # fix this to include placement
                    builder.init_seg(segid)
                    builder.init_residue(resname, hetero, id, \
                            __default_icode__)
                    for atom in molecule:
                        try: bfactor = atom.bfactor
                        except AttributeError: bfactor = __default_bfactor__
                        try: occupancy = atom.occupancy
                        except AttributeError: occupancy = __default_occupancy__
                        try: altloc = atom.altloc
                        except AttributeError: altloc = __default_altloc__
                        try: fullname = atom.fullname
                        except AttributeError: fullname = \
                                string.center(str(atom), 4)
                        builder.init_atom(str(atom), atom.position, \
                                bfactor, occupancy, altloc, fullname)
                    id += 1
        # log messages
#        sys.stdout = so
        sys.stderr = se
#        for line in ostream.readlines():
#            loki.debug('Bio.PDB -> %s' % line)
        for line in estream.readlines():
            loki.warn('Bio.PDB -> %s' % line)
        return builder.get_structure()

    def write_remarks(self, fp):
        ro = re.compile('^(MODEL|ATOM)')
        # rewind file position before reading
        fp.seek(0)
        lines = fp.readlines()
        # rewind again to overwrite previous data
        fp.seek(0)
        inserted = False
        for line in lines:
            if inserted or not ro.match(line):
                fp.write(line)
            else:
                for remark in self.remarks:
                    fp.write(remark)
                inserted = True
                fp.write(line)

    def write(self, filename):
        """
        Write a PDB file based on the ensemble stored in self.ensemble.

        arguments:
          filename -- PDB file to be written
        """
        # use a temporary file in case we need to edit the output
        fp = TemporaryFile()
        self.__writer__.set_structure(self.__build_structure__())
        self.__writer__.save(fp)
        if self.remarks:
            self.write_remarks(fp)
        # rewind file position before reading
        fp.seek(0)
        txt = fp.read()
        fp.close()
        if hasattr(filename, 'write'):
            filename.write(txt)
        else:
            fp = open(filename, 'w')
            fp.write(txt)
            fp.close()



class rawPDB:
    """
    DEFUNC! Should read and write PDBs with no outside help.
    """

    def __init__(self, filename=None, verbose=False, debug=False):
        loki.warn('Biopython not installed. Using a defunctional version ' + \
                'of PDB parsing. == DIVE INTO THE CODE AND FIXME!')
        self.header = []
        self.title = []
        self.author = []
        self.verbose = bool(verbose)
        self.debug = bool(debug)
        self.atoms = []
        self.residues = []
        self.chains = []
        self.models = []
        self.ensemble = None
        if filename != None:
            self.read(filename)

    def read(self, filename, discard=[], ignore_models=True):
        # TODO: change translations to a more generic form
        re_space = '[ \t]+'
        re_deci = '[+-]?\d*\.\d+'
        re_atom = re_space.join(['', '\d+','([0-9A-Z]+)', \
                '([A-Z]{3})[ \tA-Z]+(\d+)'] \
                + ['(%s)'%re_deci]*5) + '[ \t0-9A-Z]+\n'
        re_header = '^HEADER'
        re_title = '^TITLE'
        for i in range(len(discard)):
            discard[i] = omniTranslator(str(discard[i]), 'atom')
        file = open(filename)
        m = re.match
        s = re.search
        models = []
        model = None
        chains = []
        chain = None
        residues = []
        reskeys = []
        hetmols = []
        hetkeys = []
        for line in file.readlines():
            if m('^HEADER', line):
                classification = line[10:50].strip()
                depDate = line[50:59].strip()
                idCode = line[62:66].strip()
                self.header.append((classification, depDate, idCode))
            elif m('^TITLE', line):
                self.title.append(line[10:70].strip())
            elif m('^AUTHOR', line):
                self.author.append(line[10:70].strip())
            elif m('^MODEL', line):
                model = int(line[6:].strip())
            elif m('^ATOM', line):
                try:
                    nucleus, residue, id, x, y, z, occ, temp = \
                            s(re_atom, line[6:]).groups()
                except AttributeError:
                    print line[6:]
                name = omniTranslator.base(nucleus, 'atom')
                atom = Atom(allNuclei[name], 
                        name=nucleus,
                        position=(x, y, z))
#                        name=omniTranslator(nucleus, 'atom'), 
                atom.occupancy = occ
                atom.tempFactor = temp
                if str(atom) not in discard:
                    if id in reskeys:
                        residues[reskeys.index(id)].append(atom)
                    else:
                        residues.append(Residue(allAminoAcids[residue], id, i\
                                [atom]))
                        reskeys.append(id)
            elif m('^HETATM', line):
                nucleus, name, id, x, y, z, occ, temp = \
                        s(re_atom, line[6:]).groups()
                name = omniTranslator.base(nucleus, 'atom')
                atom = Atom(allNuclei[name], 
                        name=nucleus,
                        position=(x, y, z))
#                        name=omniTranslator(nucleus, 'atom'), 
                atom.occupancy = occ
                atom.tempFactor = temp
                if str(atom) not in discard:
                    if id in hetkeys:
                        hetmols[hetkeys.index(id)].append(atom)
                    else:
                        hetmols.append(Molecule(name, [atom]))
                        hetkeys.append(id)
            elif m('^TER', line) or m('^END', line):
                if chain == None:
                    name = str(len(chains)+1)
                if name not in discard:
                    chains.append(Polypeptide(name, residues))
                residues = []
                reskeys = []
            elif m('^ENDMDL', line) and not ignore_models:
                models.append(Ensemble(model, chains+hetmols))
                chains = []
                hetmols = []
                hetkeys = []
            else:
                if self.verbose:
                    print 'Unknown entry: %s'%repr(line)
        if len(self.models) == 0 and (len(chains) or len(hetmols)):
            models.append(Ensemble('default', chains+hetmols))
        self.models = models
        if len(self.models):
            self.ensemble = self.models[0]
#        else:
#            raise ArgumentError, 'Empty PDB? No structures loaded.'

    def getModel(self, name):
        for model in self.models:
            if str(model) == name: return model
        return None

    def use_model(self, i):
        if type(i) != int:
            for j in range(len(self.models)):
                if str(self.models[j]) == i:
                    m = j
                    break
        else:
            m = i
        self.chains = []
        self.hetatoms = []
        for x in models[m]:
            if classify(x) == 'Chain':
                self.chains.append(x)
            elif classify(x) == 'Molecule':
                self.hetatoms.append(x)
        self.residues = []
        for x in self.chains:
            self.residues.extend(x[:])
        self.atoms = []
        for x in self.residues:
            self.atoms.extend(x[:])

    def view(self):
        cin, cout, cerr = os.popen3('vmd %s'%name, os.P_NOWAIT)

class GeneralFixed:
    def __init__(self, width=78, template='', key_alias='', \
            format_alias='', keys=[], format='', comments=[]):
        loki.debug('GeneralFixed() < ' + ' | '.join([repr(x) for x in \
                [width, template, key_alias, format_alias, keys, format, \
                comments]]))
        self.width = int(width)
        self.template = str(template)
        self.key_alias = str(key_alias)
        self.format_alias = str(format_alias)
        self.keys = keys
        self.format = str(format)
        self.comments = comments
    def __extract_comment__(self, alias, stream, strip=True, shift=False, \
            glue=' '):
        alias = escape_parentheses(alias)
        comment = re.findall('%s(.+)\n' % alias, stream)
        snippets = re.findall('(%s.+\n)' % alias, stream)
        if strip:
            comment = glue.join([x.strip() for x in comment])
        elif shift:
            comment = glue.join([x[1:] for x in comment])
        return comment, snippets

#    def __format_regex__(self, format):
#        decimals = re.findall('(%([0-9]*)d)', format)
#        for dec in decimals:
#            if dec[1][0] == '0':
#                n = int(dec[1][1:])
#                r = '([0-9]{%d,}|[+-][0-9]{%d,})' % (n, n - 1)
#            elif dec[1]:
#                n = int(dec[1])
#                r = ['[0-9]{%d,}' % n, '[+-][0-9]{%d,}' % (n - 1)]
#                for i in range(1,n):
#                    r.append('[ ]{%d}[0-9]{%d,}' % (i, n-i))
#                    if i < n - 1:
#                        r.append('[ ]{%d}[+-][0-9]{%d,}' % (i, n-i))
#                r = '(' + '|'.join(r) + ')'
#            else:
#                r = '([+-]?[0-9]+)'
#            regex = format.replace(dec[0], r)
#        strings = re.findall('(%([0-9]*)s)', format)
#        for string in strings:
#            if string[1]:
#                r = '([^\t\n\a\b\f\r\v]{%d,})' % int(string[1])
#            else:
#                r = '(\w+)'
#            regex = regex.replace(string[0], r)
#        floats = re.findall('(%([0-9]*)\.?([0-9]*)f)', format)
#        for fl in floats:

    def __format_mapper__(self, x):
        mapping = {int:'diu', float:'eEfFgG', str:'oxXcrs'}
        for k,v in mapping.items():
            if x in mapping[k]:
                return k
        raise ValueError, "unknown format code '%s'" % x
    def read(self, filename, key_alias=None, format_alias=None, keys=None, \
            format=None, comments=[], wild_comments=False):
        loki.debug('GeneralFixed.read() < ' + ' | '.join([repr(x) for x in \
                [filename, key_alias, format_alias, keys, format, comments, \
                wild_comments]]))
        keys = keys or self.keys
        format = format or self.format
        key_alias = key_alias or self.key_alias
        format_alias = format_alias or self.format_alias
        comments = comments or self.comments
        f = open(filename)
        stream = f.read()
        f.close()
        if not keys:
            if not (self.key_alias and stream.count(self.key_alias)):
                raise ValueError, \
                        'either keys or a valid keyalias must be defined'
            else:
                keys, snips = self.__extract_comment__(self.key_alias, stream)
                keys = keys.split()
                for s in snips:
                    stream = stream.replace(s, '')
        if not format:
            if not (self.format_alias and stream.count(self.format_alias)):
                raise ValueError, \
                        'either format or a valid formatalias must be defined'
            else:
                format, snips = self.__extract_comment__(self.format_alias, \
                        stream, strip=False, shift=True)
                for s in snips:
                    stream = stream.replace(s, '')
        loki.debug('keys=%s', keys)
        loki.debug('format=%s', format)
        places = re.findall('%[0-9.]+([diouxXeEfFgGcrs])', format)
        loki.debug('places=%s', places)
        funks = map(self.__format_mapper__, places)
        loki.debug('funks=%s', funks)
        data = {'comments':{}, 'table':[]}
        for c in comments:
            comment, snips = self.__extract_comment__(c, stream)
            loki.debug("comment " + repr(c) + ": " + repr(comment))
            data['comments'][c] = comment
            for s in snips:
                stream = stream.replace(s, '')
        for line in stream.split('\n'):
            nonitem = False
            items = line.split()
            if len(items) == len(funks):
                try:
                    items = [f(x) for f,x in zip(funks, items)]
                except ValueError:
                    nonitem = True
            else:
                nonitem = True
            if nonitem:
                if wild_comments:
                    key = re.search('(^[ \t]*[a-zA-Z0-9]+)[ \t]+[a-zA-Z0-9]', \
                            line)
                    if key:
                        s = line.remove(key).strip()
                        data['comments'].setdefault(key, '')
                        if data['comments'][key]:
                            s = ' ' + s
                        data['comments'][key] += s
                    else:
                        loki.debug("ignoring '%s'" % line)
                else:
                    loki.debug("ignoring '%s'" % line)
            else:
                data['table'].append(dict(zip(keys, items)))
        return data
    def write(self, data, filename, keys=None, format=None):
        keys = keys or self.keys
        format = format or self.format
        if not keys or not format:
            raise ValueError, 'keys and format must be defined'
        if isinstance(filename, file):
            f = filename
        else:
            f = open(filename, 'w')
        template = self.template
        for key, value in data['comments'].items():
            if self.template:
                if re.search('\${%s, [0-9]+}' % key, self.template):
                    n = int(re.search('\$\{%s, ([0-9]+)\}' % key, \
                            self.template).groups()[0])
                    room = n - len(key) - 1
                    tag = '${%s, %d}' % (key, n)
                elif self.template.count('${%s}' % key):
                    room = self.width - len(key) - 1
                    tag = '${%s}' % key
                else:
                    raise ValueError, "comment '%s' not found in template" \
                            % str(key)
            else:
                room = self.width - len(key) - 1
            snips = []
            while value:
                snips.append(value[:room])
                value = value.replace(value[:room], '')
            if template:
                all = ''
                for s in snips:
                    all += key + ' ' + s + '\n'
                loki.debug("replacing '%s' w/ '%s' in PalesIO template", \
                        repr(tag), repr(all))
                template = template.replace(tag + '\n', all)
            else:
                for s in snips:
                    f.write(key + ' ' + s + '\n')
        if self.key_alias:
            line = self.key_alias + ' ' + ' '.join(keys) + '\n'
            if template:
                template = template.replace('${keys}\n', line)
            else:
                f.write(line)
        if self.format_alias:
            line = self.format_alias + ' ' + format + '\n'
            if template:
                template = template.replace('${format}\n', line)
            else:
                f.write(line)
        if template:
            f.write(template)
        for item in data['table']:
            values = []
            for key in keys:
                values.append(item[key])
            f.write(format % tuple(values) + '\n')
        f.close()


class PalesIO:
    def __init__(self, filename=None, \
            keys=['RESID_I', 'RESNAME_I', 'ATOMNAME_I', 'RESID_J', \
            'RESNAME_J', 'ATOMNAME_J', 'D', 'DD', 'W'], \
            format='%5d  %6s  %6s  %5d  %6s  %6s  %9.3f  %9.3f  %.2f', \
            key_alias='VARS', format_alias='FORMAT', \
            template='${DATA SEQUENCE, 69}\n\n${keys}\n${format}\n\n', \
            filter=None):
        self.keys = keys
        self.format = str(format)
        self.key_alias = str(key_alias)
        self.format_alias = str(format_alias)
        self.template = str(template)
        self.parser = GeneralFixed(width=72, key_alias=self.key_alias, \
                format_alias=self.format_alias, keys=self.keys, \
                format=self.format, comments=['DATA SEQUENCE'], \
                template=self.template)
        self.__coupling_map__ = {}
        for name, dc in DCouplings.items():
            self.__coupling_map__[dc.atoms] = name
        if not filter:
            filter = lambda x: True
        self.__filter__ = filter
        self.__filename__ = ''
        if filename:
            self.couplings, self.seq = self.read(filename)
        else:
            self.couplings, self.seq = Couplings(), []
    def __str__(self):
        return self.__filename__
    def __repr(self):
        return "PalesIO('" + str(self) + "', keys=" + repr(self.keys) + \
                ", format=" + repr(self.format) + ", key_alias=" + \
                repr(self.key_alias) + ", format_alias=" + \
                repr(self.format_alias) + ', template=' + \
                repr(self.template) + ", filter=" + repr(self.filter) + ")"

    def read(self, filename, filter=None):
        loki.debug('PalesIO.read() < %s | %s' % (repr(filename), \
                repr(filter)))
        self.__filename__ = filename
        loki.debug("Reading file '%s'." % str(self))
        filter = filter or self.__filter__
        data = self.parser.read(filename)
        if data['comments'].has_key('DATA SEQUENCE'):
            seq = data['comments']['DATA SEQUENCE']
            self.seq = [x for x in seq.replace(' ', '')]
        self.couplings = Couplings()
        switch_count = 0
        for item in data['table']:
            atom1, atom2 = item['ATOMNAME_I'], item['ATOMNAME_J']
            if atom1 == 'H':
                atom1 = 'HN'
                switch_count += 1
            if atom2 == 'H':
                atom2 = 'HN'
                switch_count += 1
            if atom1 not in libLingua.dictionaryAtom:
                raise ValueError, "unknown atom '%s'" % atom1
            if atom2 not in libLingua.dictionaryAtom:
                raise ValueError, "unknown atom '%s'" % atom2
            if self.__coupling_map__.has_key((atom1, atom2)):
                id1, id2 = item['RESID_I'], item['RESID_J']
                res1, res2 = item['RESNAME_I'], item['RESNAME_J']
            elif self.__coupling_map__.has_key((atom2, atom1)):
                id2, id1 = item['RESID_I'], item['RESID_J']
                res2, res1 = item['RESNAME_I'], item['RESNAME_J']
                atom2, atom1 = atom1, atom2
            else:
                raise ValueError, "unknown atom combination '%s, %s'" % \
                        (repr(item['ATOMNAME_I']), repr(item['ATOMNAME_J']))
            name = self.__coupling_map__[(atom1, atom2)]
            if id1 < id2:
                shift = (0, id2 - id1)
            else:
                shift = (id1 - id2, 0)
            if shift != DCouplings[name].shift:
                loki.warn("something fishy in residue ids '%s, %s' for " + \
                        "atoms '%s, %s'", str(id1), str(id2), atom1, atom2)
            dc = Coupling(item['D'], DCouplings[name], id2, \
                    residues=(res1.capitalize(), res2.capitalize()), \
                    error=item['DD'])
            if filter(dc):
                self.couplings.append(dc)
                if item['W']:
                    self.couplings[-1].weight = item['W']
        if switch_count:
            loki.info('switched %d H to HN' % switch_count)
        loki.debug("File '%s' read." % str(self))
        return self.couplings, self.seq

    def write(self, filename, couplings=None, seq=None, errors={}, \
            default_error=0.5, weights={}, default_weight=1.0):
        loki.debug('PalesIO.write() < ' + ' | '.join([repr(x) for x in \
                [filename, couplings, seq, errors, default_error, weights, \
                default_weight]]))
        couplings = couplings or self.couplings
        seq = seq or self.seq
        for c in couplings:
            if not hasattr(c, 'weight'):
                if weights.has_key(str(c)):
                    c.weight = float(weights[str(c)])
                else:
                    c.weight = float(default_weight)
            if not c.error:
                if errors.has_key(str(c)):
                    c.error = float(weights[str(c)])
                else:
                    c.error = float(default_error)
        sseq = ''
        for i in range(len(seq)):
            sseq += omniTranslator(seq[i], '1-letter-aminoacids')
            if not (len(sseq) + 1) % 11: sseq += ' '
        data = {'comments':{'DATA SEQUENCE':sseq}, 'table':[]}
        for c in couplings:
            res1, res2 = c.getResidues()
            atom1, atom2 = c.type().atoms
            if atom1 == 'HN': atom1 = 'H'
            if atom2 == 'HN': atom2 = 'H'
            shift1, shift2 = c.type().shift
            item = {'RESID_I': int(c) + shift1, \
                    'RESNAME_I': res1.upper(), \
                    'ATOMNAME_I': atom1.upper(), \
                    'RESID_J': int(c) + shift2, \
                    'RESNAME_J': res2.upper(), \
                    'ATOMNAME_J': atom2.upper(), \
                    'D': float(c), \
                    'DD': c.error, \
                    'W': c.weight}
            data['table'].append(item)
        loki.debug("Writing file '%s'." % str(filename))
        self.parser.write(data, filename)


__pales_output_template__ = """REMARK Molecular Alignment Simulation.

REMARK Simulation parameters.

${DATA PALES_MODE}

${DATA PALES LC_TYPE}
${DATA PALES LC_CONCENTRATION}
${DATA PALES ORIENT_SPHERE}
${DATA PALES ORIENT_PSI}
${DATA PALES GRID_SPACING}
${DATA PALES MODEL_RADIUS}
${DATA PALES LC_ORDER}
${DATA PALES ATOM_RADIUS}
${DATA PALES SEL_SIMPLE_FLAG}
${DATA PALES SURF_FLAG}


REMARK Order matrix.

${DATA SAUPE_MATRIX}
${DATA SAUPE}

${DATA IRREDUCIBLE_REP}
${DATA IRREDUCIBLE}
${DATA IRREDUCIBLE GENERAL_MAGNITUDE}


REMARK Mapping of coordinates.

${DATA MAPPING_COOR}
${DATA MAPPING}
${DATA MAPPING INV}


REMARK Eigensystem & Euler angles for clockwise rotation about z, y', z''.

${DATA EIGENVALUES (Sxx_d,Syy_d,Szz_d)}
${DATA EIGENVECTORS          (x_coor      y_coor      z_coor)}
${DATA EIGENVECTORS X_AXIS}
${DATA EIGENVECTORS Y_AXIS}
${DATA EIGENVECTORS Z_AXIS}

${DATA Q_EULER_SOLUTIONS}
${DATA Q_EULER_ANGLES  1}
${DATA Q_EULER_ANGLES  2}
${DATA Q_EULER_ANGLES  3}
${DATA Q_EULER_ANGLES  4}


REMARK Euler angles (psi/theta/phi) for rotation about x, y, z.

${DATA EULER_SOLUTIONS}
${DATA EULER_ANGLES}
${DATA EULER_ANGLES}

${DATA Da}
${DATA Dr}


${DATA Aa}
${DATA Ar}


${DATA Da_HN}
${DATA Rhombicity}


REMARK Dipolar couplings.

${DATA N}
${DATA RMS}
${DATA Chi2}
${DATA CORR R}
${DATA Q SAUPE}
${DATA REGRESSION OFFSET}
${DATA REGRESSION SLOPE}
${DATA REGRESSION BAX SLOPE}

${keys}
${format}

"""

class PalesOutputIO:
    def __init__(self, filename=None, \
            keys=['RESID_I', 'RESNAME_I', 'ATOMNAME_I', 'RESID_J', \
            'RESNAME_J', 'ATOMNAME_J', 'DI', 'D_OBS', 'D', 'D_DIFF', \
            'DD', 'W'], format=\
            '%4d %4s %4s %4d %4s %4s %9.2f %9.3f %9.3f %9.3f %.2f %.2f', \
            key_alias='VARS', format_alias='FORMAT', template=None, \
            filter=None, comments=['REMARK', 'DATA PALES_MODE', \
            'DATA PALES LC_TYPE', 'DATA PALES LC_CONCENTRATION', \
            'DATA PALES ORIENT_SPHERE', 'DATA PALES ORIENT_PSI', \
            'DATA PALES GRID_SPACING', 'DATA PALES MODEL_RADIUS', \
            'DATA PALES LC_ORDER', 'DATA PALES ATOM_RADIUS', \
            'DATA PALES SEL_SIMPLE_FLAG', 'DATA PALES SURF_FLAG', \
            'DATA SAUPE_MATRIX', 'DATA SAUPE', 'DATA IRREDUCIBLE_REP', \
            'DATA IRREDUCIBLE', 'DATA IRREDUCIBLE GENERAL_MAGNITUDE', \
            'DATA MAPPING_COOR', 'DATA MAPPING', 'DATA MAPPING INV', \
            'DATA EIGENVALUES (Sxx_d,Syy_d,Szz_d)', \
            'DATA EIGENVECTORS          (x_coor      y_coor      z_coor)', \
            'DATA EIGENVECTORS X_AXIS', 'DATA EIGENVECTORS Y_AXIS', \
            'DATA EIGENVECTORS Z_AXIS', 'DATA Q_EULER_SOLUTIONS', \
            'DATA Q_EULER_ANGLES  1', 'DATA Q_EULER_ANGLES  2', \
            'DATA Q_EULER_ANGLES  3', 'DATA Q_EULER_ANGLES  4', \
            'DATA EULER_SOLUTIONS', 'DATA EULER_ANGLES', \
            'DATA Da', 'DATA Dr', 'DATA Aa', 'DATA Ar', 'DATA Da_HN', \
            'DATA Rhombicity', 'DATA N', 'DATA RMS', 'DATA Chi2', \
            'DATA CORR R', 'DATA Q SAUPE', 'DATA REGRESSION OFFSET', \
            'DATA REGRESSION SLOPE', 'DATA REGRESSION BAX SLOPE'], \
            seq=None):
        self.keys = keys
        self.format = str(format)
        self.key_alias = str(key_alias)
        self.format_alias = str(format_alias)
        if template:
            self.template = str(template)
        else:
            self.template = __pales_output_template__
        self.comments = comments
        self.seq = seq or []
        self.parser = GeneralFixed(width=72, key_alias=self.key_alias, \
                format_alias=self.format_alias, keys=self.keys, \
                format=self.format, comments=self.comments, \
                template=self.template)
        self.__coupling_map__ = {}
        for name, dc in DCouplings.items():
            self.__coupling_map__[dc.atoms] = name
        if not filter:
            filter = lambda x: True
        self.__filter__ = filter
        if filename:
            self.couplings, self.meta = self.read(filename)
        else:
            self.couplings, self.meta = Couplings(), {}
    def __str__(self):
        return self.__filename__
    def __repr(self):
        return "PalesOutputIO('" + str(self) + "', keys=" + repr(self.keys) + \
                ', format=' + repr(self.format) + ', key_alias=' + \
                repr(self.key_alias) + ', format_alias=' + \
                repr(self.format_alias) + ', template=' + \
                repr(self.template) + ', filter=' + repr(self.filter) + \
                ', comments=' + repr(self.comments) + ')'

    def read(self, filename, filter=None, store=True):
        loki.debug('PalesOutputIO.read() < %s | %s | %s' % (repr(filename), \
                repr(filter), repr(store)))
        self.__filename__ = filename
        loki.debug("Reading file '%s'." % str(self))
        filter = filter or self.__filter__
        data = self.parser.read(filename)
# TODO: remove below
        fp = open(filename)
        txt = fp.read()
        pales2 = 'MAPPING' in txt
# TODO: remove above
        meta = {}
        for comment in data['comments'].keys():
            if comment[:5] == 'DATA ':
                meta[comment[5:]] = data['comments'][comment]
# TODO: remove this
        if pales2:
            meta['MAPPING'], meta['MAPPING INV'] = \
                    meta['MAPPING'].split('INV')
            meta['MAPPING'] = meta['MAPPING'].strip()
            meta['MAPPING INV'] = meta['MAPPING INV'].strip()
            meta['Da'], meta['Da_HN'] = meta['Da'].split('_HN')
            meta['Da'] = meta['Da'].strip()
            meta['Da_HN'] = meta['Da_HN'].strip()
        couplings = Couplings()
        switch_count = 0
        for item in data['table']:
            atom1, atom2 = item['ATOMNAME_I'], item['ATOMNAME_J']
            if atom1 == 'H':
                atom1 = 'HN'
                switch_count += 1
            if atom2 == 'H':
                atom2 = 'HN'
                switch_count += 1
            if atom1 not in libLingua.dictionaryAtom:
                raise ValueError, "unknown atom '%s'" % atom1
            if atom2 not in libLingua.dictionaryAtom:
                raise ValueError, "unknown atom '%s'" % atom2
            if self.__coupling_map__.has_key((atom1, atom2)):
                id1, id2 = item['RESID_I'], item['RESID_J']
                res1, res2 = item['RESNAME_I'], item['RESNAME_J']
            elif self.__coupling_map__.has_key((atom2, atom1)):
                id2, id1 = item['RESID_I'], item['RESID_J']
                res2, res1 = item['RESNAME_I'], item['RESNAME_J']
                atom2, atom1 = atom1, atom2
            else:
                raise ValueError, "unknown atom combination '%s, %s'" % \
                        (repr(item['ATOMNAME_I']), repr(item['ATOMNAME_J']))
            name = self.__coupling_map__[(atom1, atom2)]
            if (id1 - id2, 0) != DCouplings[name].shift:
                loki.warn("something fishy in residue ids '%s, %s' for " + \
                        "atoms '%s, %s'", str(id1), str(id2), atom1, atom2)
            dc = Coupling(item['D'], DCouplings[name], id2, \
                    residues=(res1.capitalize(), res2.capitalize()), \
                    error=item['DD'])
            if filter(dc):
                couplings.append(dc)
                if item['W']:
                    couplings[-1].weight = item['W']
                if item['DI']:
                    couplings[-1].pales_scale = item['DI']
                loki.debug('dc=%s' % repr(dc))
        if switch_count:
            loki.info('switched %d H to HN' % switch_count)
        if store:
            self.meta = meta
            self.couplings = couplings
        loki.debug("File '%s' read." % str(self))
        return couplings, meta

    def write(self, filename, couplings=None, meta={}, observed=None, \
            errors={}, default_error=0.5, weights={}, default_weight=1.0, \
            scales={}, default_scale=1.0):
        loki.debug('PalesOutputIO.write() < ' + ' | '.join([repr(x) for x in \
                [filename, couplings, meta, observed, errors, default_error, \
                weights, default_weight, scales, default_scale]]))
        couplings = couplings or self.couplings
        if not observed:
            observed = [0.0] * len(couplings)
        elif observed in [int, float, str]:
            observed = [float(observed)] * len(couplings)
        elif len(observed) != len(couplings):
            loki.warn('Coupling sets do not match. ' + \
                    'Discarding observed couplings.')
            observed = [0.0] * len(couplings)
        else:
            observed = [float(x) for x in observed]
        meta = meta or self.meta
        for c in couplings:
            if not hasattr(c, 'weight'):
                if weights.has_key(str(c)):
                    c.weight = float(weights[str(c)])
                else:
                    c.weight = float(default_weight)
            if not c.error:
                if errors.has_key(str(c)):
                    c.error = float(weights[str(c)])
                else:
                    c.error = float(default_error)
            if not hasattr(c, 'pales_scale'):
                if scales.has_key(str(c)):
                    c.pales_scale = float(scales[str(c)])
                else:
                    c.pales_scale = float(default_scale)
        comments = {}
        for m,v in meta.items():
            comments['DATA ' + m] = v
        data = {'comments':comments, 'table':[]}
        for c,o in zip(couplings, observed):
            res1, res2 = c.getResidues()
            atom1, atom2 = c.type().atoms
            if atom1 == 'HN': atom1 = 'H'
            if atom2 == 'HN': atom2 = 'H'
            shift1, shift2 = c.type().shift
            item = {'RESID_I': int(c) + shift1, \
                    'RESNAME_I': res1.upper(), \
                    'ATOMNAME_I': atom1.upper(), \
                    'RESID_J': int(c) + shift2, \
                    'RESNAME_J': res2.upper(), \
                    'ATOMNAME_J': atom2.upper(), \
                    'DI': c.pales_scale, \
                    'D_OBS': o, \
                    'D': float(c), \
                    'D_DIFF': o - float(c), \
                    'DD': c.error, \
                    'W': c.weight}
            data['table'].append(item)
        loki.debug("Writing file '%s'." % str(filename))
        self.parser.write(data, filename)


class TextFile:
    def __init__(self, filename=None, ignore_comments=True, delimiter='#', \
            meta_reco='(\$([\w\d]+)=([\w\d]+|".+"|'+"'.+'))", mode='r'):
        self.ignore_comments = ignore_comments
        if isinstance(delimiter, str):
            self.delimiter = delimiter
        else:
            loki.warn("using default comment delimiter '#'")
            self.delimiter = '#'
        self.__meta_reco__ = re.compile(meta_reco)
        if filename:
            self.open(filename, mode)
        else:
            self.set_mode(mode)
            self.__fp__ = None
        self.flush()
    def __iter__(self):
        return self
    def next(self):
        try:
            line = self.__fp__.next()
        except StopIteration:
            raise StopIteration
        if self.ignore_comments:
            line = self.remove_comments(line)
        return line
    def remove_comments(self, line):
        if not line or not isinstance(line, str):
            return ''
        ending = ''
        if line[-1] == '\n':
            ending = '\n'
            line = line[:-1]
        parts = line.split(self.delimiter, 1)
        if len(parts) > 1:
            line = parts[0].rstrip()
            meta, comment = self.pick_meta(parts[1])
            if comment:
                self.comments.append(comment)
            if meta:
                self.meta[meta[0]] = meta[1]
        return line + ending
    def open(self, filename, mode='r'):
        self.set_mode(mode)
        self.__fp__ = open(filename, self.__mode__)
    def close(self):
        self.__fp__.close()
    def set_mode(self, mode='r'):
        if mode in ['r','w','a', 'r+', 'rb', 'wb', 'ab', 'r+b']:
            self.__mode__ = mode
        else:
            loki.warn("Ignoring unknown file operation mode '%s'." % mode)
    def flush(self):
        self.comments = []
        self.meta = {}
        if hasattr(self.__fp__, 'seek'):
            self.__fp__.seek(0)
    def readline(self):
        try:
            return self.next()
        except StopIteration:
            return None
    def readlines(self):
        self.flush()
        lines = []
        for line in self:
            lines.append(line)
        return lines
    def read(self):
        self.flush()
        txt = ''
        for line in self:
            txt += line
        return txt
    def pick_meta(self, comment):
        ro = self.__meta_reco__.search(comment)
        meta = None
        if ro:
            comment = comment.replace(ro.group(1), '').strip()
            meta = (ro.group(2), ro.group(3).strip('"' + "'"))
        return meta, comment

class GromosIO:
    def __init__(self, filename=None):
        if filename is not None:
            self.__noread__ = False
            self.file = FileCrawler(filename)
            self.read_meta()
        else:
            self.__noread__ = True
        self.reset()
        self.__max__ = 99999 # max integer to fit in the columns
    def read_meta(self):
        if self.__noread__:
            raise IOError, 'no file open for reading'
        self.file.wind()
        self.file.prevline()
        self.box = None
        while self.box is None:
            mark = self.file.tell()
            line = self.file.readline()
            if line.strip():
                self.box = [float(i) for i in line.split()]
            else:
                self.file.seek(meta)
                self.file.prevline()
        self.file.rewind()
        self.title = self.file.readline().rstrip()
        self.count = int(self.file.readline().strip())
        self.__reset_point__ = self.file.tell()
        return self.title, self.count, self.box
    def read_atom(self):
        if self.__noread__:
            raise IOError, 'no file open for reading'
        line = self.file.readline()
        if len(line.split()) in [3,0]:
            raise EOF
        record = {}
        record['resno'] = int(line[:5])
        record['resname'] = line[5:10]
        record['atomname'] = line[10:15]
        record['atomno'] = int(line[15:20])
        xv = line[20:].split()
        record['position'] = [float(xv[0]), float(xv[1]), float(xv[2])]
        if len(xv) > 3:
            record['velocity'] = [float(xv[3]), float(xv[4]), float(xv[5])]
        if self.__reada__ == self.__max__ and record['atomno'] == 0:
            self.__wrapa__ += 1
        if self.__readr__ == self.__max__ and record['resno'] == 0:
            self.__wrapr__ += 1
        self.__reada__ = record['atomno']
        self.__readr__ = record['resno']
        record['atomno'] += self.__wrapa__ * (self.__max__ + 1)
        record['resno'] += self.__wrapr__ * (self.__max__ + 1)
        return record
    def read_residue(self):
        if self.__noread__:
            raise IOError, 'no file open for reading'
        try:
            atoms = [self.read_atom()]
        except EOF:
            raise EOF
        id = atoms[0]['resno']
        loki.debug("reading residue '%d'" % id)
        while id == atoms[-1]['resno']:
            try:
                atoms.append(self.read_atom())
            except EOF:
                return atoms
        # retrace one step
        del atoms[-1]
        self.file.prevline()
        if self.__readr__ == 0:
            self.__wrapr__ -= 1
        if self.__reada__ == 0:
            self.__wrapa__ -= 1
        self.__reada__ = atoms[-1]['atomno']
        self.__readr__ = id
        return atoms
    def read_atoms(self):
        if self.__noread__:
            raise IOError, 'no file open for reading'
        atoms = []
        while True:
            try:
                atoms.append(self.read_atom())
            except EOF:
                break
        return atoms
    def read_residues(self):
        if self.__noread__:
            raise IOError, 'no file open for reading'
        residues = []
        while True:
            try:
                residues.append(self.read_residue())
            except EOF:
                break
        return residues
    def reset(self):
        self.__ai__ = 0
        self.__ri__ = 0
        self.__reada__ = None
        self.__readr__ = None
        self.__preva__ = None
        self.__prevr__ = None
        self.__wrapa__ = 0
        self.__wrapr__ = 0
        if not self.__noread__:
            self.file.seek(self.__reset_point__)
    def declare(self, residue):
        if residue != self.__prevr__:
            self.__ri__ += 1
            self.__prevr__ = residue
        self.__ai__ += 1
        return self.__ai__ % (self.__max__ + 1), \
                self.__ri__ % (self.__max__ + 1)
    def set_residue_index(self, i):
        self.__ri__ = i
    def set_atom_index(self, i):
        self.__ai__ = i
    def get_residue_index(self):
        return self.__ri__
    def get_atom_index(self):
        return self.__ai__
    def write_atom(self, atom):
        ai, ri = self.declare(atom['resno'])
        if atom.has_key('velocity'):
            print '%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f' % (
                    ri, atom['resname'], atom['atomname'], ai, 
                    atom['position'][0], atom['position'][1], 
                    atom['position'][2], 
                    atom['velocity'][0], atom['velocity'][1], 
                    atom['velocity'][2])
        else:
            print '%5d%-5s%5s%5d%8.3f%8.3f%8.3f' % (
                    ri, atom['resname'], atom['atomname'], ai,
                    atom['position'][0], atom['position'][1],
                    atom['position'][2])
    def write_box(self, box=None):
        box = box or self.box
        print '   ' + '   '.join(['%.5f' % x for x in box])
    def write_count(self, count=None):
        count = count or self.count
        print '  %d' % count
    def write(self, atoms=None, title=None, count=None, box=None):
        title = title or self.title
        print title
        self.write_count(len(atoms))
        atoms = atoms or self.atoms
        for atom in atoms:
            self.write_atom(atom)
        self.write_box(box)

class BlockIO:
    def __init__(self, filename=None, bra='[', ket=']'):
        if filename is not None:
            self.__noread__ = False
            self.file = FileCrawler(filename)
        else:
            self.__noread__ = True
        self.__marks__ = {}
        self.__found_all__ = False
        self.__order__ = []
        self.__bra__ = str(bra)
        self.__ket__ = str(ket)
        self.comment = ['#', '%', ';']
    def __readline__(self):
        line = self.file.readline()
        if line == '':
            raise EOF
        if line[0] in self.comment:
            return ''
        else:
            return line.strip()
    def __check__(self, line):
        if line and line[0] == self.__bra__ and line[-1] == self.__ket__:
            name = line[len(self.__bra__):-1*len(self.__ket__)].strip()
            if self.__marks__.has_key(name):
                if self.__marks__[name] != self.file.tell():
                    loki.warn("Duplicate index group '%s' found." % name)
            else:
                self.__marks__[name] = self.file.tell()
            return name
        return False
    def __find__(self, name):
        if self.__noread__:
            raise IOError, 'no file open for reading'
        if self.__marks__.has_key(name):
            return self.__marks__[name]
        while True:
            line = self.__readline__()
            if self.__check__(line) and self.__marks__.has_key(name):
                break
        return self.__marks__[name]
    def __find_all__(self):
        if self.__noread__:
            raise IOError, 'no file open for reading'
        self.file.rewind()
        while True:
            try:
                line = self.__readline__()
            except EOF:
                break
            self.__check__(line)
        self.__found_all__ = True
    def __keep_order__(self):
        if not self.__found_all__:
            self.__find_all__()
        marks = []
        for k,v in zip(self.__marks__.keys(), self.__marks__.values()):
            marks.append((v,k))
        marks.sort()
        self.__order__ = []
        for mark in marks:
            self.__order__.append(mark[1])
    def atomise(self, line):
        return [line]
    def unify(self, items):
        return items
    def is_comment(self, line):
        return len(line) > 0 and line[0] in self.comment
    def group(self, name):
        if self.__noread__:
            raise IOError, 'no file open for reading'
        try:
            mark = self.__find__(name)
        except EOF:
            raise ValueError, "Unknown group name '%s'." % name
        loki.debug('group %s @ mark %s' % (name, repr(mark)))
        self.file.seek(mark)
        items = []
        while True:
            try:
                line = self.__readline__()
            except EOF:
                break
            if line and line[0] != self.__bra__:
                if line[0] not in self.comment:
                    items.extend(self.atomise(line))
            else:
                break
        return items
    def get_order(self):
        if not self.__noread__:
            self.__keep_order__()
            return self.__order__
        else:
            return None
    def read_comments(self):
        comments = []
        while True:
            try:
                line = self.__readline__()
            except EOF:
                break
            if self.is_comment(line):
                comments.append(line)
            elif line:
                self.file.prevline()
                break
        return comments
    def read_block(self):
        name = None
        block = []
        while True:
            try:
                line = self.__readline__()
            except EOF:
                break
            x = self.__check__(line)
            if name is None:
                if x:
                    name = x
                elif self.is_comment(line):
                    raise ValueError, 'Block name is not defined.'
            else:
                if x:
                    self.file.prevline()
                    break
                else: 
                    block.extend(self.atomise(line))
        if name is None:
            raise EOF
        return (name, block)
    def read(self):
        if self.__noread__:
            raise IOError, 'no file open for reading'
        self.file.rewind()
        name = None
        blocks = {}
        while True:
            try:
                line = self.__readline__()
            except EOF:
                break
            block = self.__check__(line)
            if block:
                name = block
                blocks[name] = []
            elif line:
                if name is not None:
                    blocks[name].extend(self.atomise(line))
                else:
                    raise ValueError, 'Block name is not defined.'
        self.__found_all__ = True
        return blocks
    def write(self, blocks, filename=None):
        if filename:
            # TODO: re-direct STDOUT etc.
            print 'WARNING: not yet implemented'
        if not self.__noread__:
            self.__keep_order__()
        extra = complement(blocks.keys(), self.__order__)
        extra.sort()
        keys = intersection(self.__order__, blocks.keys())
        keys.extend(extra)
        for key in keys:
            print self.__bra__, key, self.__ket__
            for line in self.unify(blocks[key]):
                print line
        # TODO: restore STDOUT etc.
        #if filename:

class IndexIO(BlockIO):
    def atomise(self, line):
        return [int(x) for x in line.split()]
    def unify(self, items):
        lines = []
        buffer = []
        for id in items:
            buffer.append('%4d' % id)
            if len(buffer) == 15:
                lines.append(' '.join(buffer))
                buffer = []
        if len(buffer):
            lines.append(' '.join(buffer))
        return lines

class ITPIO(BlockIO):
    def __init__(self, filename=None, bra='[', ket=']'):
        if filename is not None:
            self.__noread__ = False
            self.file = FileCrawler(filename)
        else:
            self.__noread__ = True
        self.__marks__ = {}
        self.__found_all__ = False
        self.__order__ = []
        self.__bra__ = str(bra)
        self.__ket__ = str(ket)
        self.comment = ['#', '%', ';']
        self.header = []
        self.template = '%6d %10s %6d %6s %6s %6d %10.3f %10.4   ; %s'
    def __readline__(self):
        line = self.file.readline()
        if line == '':
            raise EOF
        else:
            return line.strip()
    def read(self):
        if self.__noread__:
            raise IOError, 'no file open for reading'
        self.file.rewind()
        name = None
        blocks = {}
        while True:
            try:
                line = self.__readline__()
            except EOF:
                break
            block = self.__check__(line)
            if block:
                name = block
                blocks.setdefault(name, [])
#                blocks[name] = []
            elif line:
                if name is not None:
                    blocks[name].extend(self.atomise(line))
                elif self.is_comment(line):
                    self.header.append(line)
                else:
                    raise ValueError, 'Block name is not defined.'
        self.__found_all__ = True
        return blocks
    def keep_order(self):
        self.__keep_order__()
    def write(self, blocks, filename=None):
        if filename:
            # TODO: re-direct STDOUT etc.
            print 'WARNING: not yet implemented'
        if not len(self.__order__):
            self.__keep_order__()
        # echo header
        for line in self.header:
            print line
        print
        # echo blocks
        extra = complement(blocks.keys(), self.__order__)
        extra.sort()
        keys = intersection(self.__order__, blocks.keys())
        keys.extend(extra)
        for key in keys:
            print self.__bra__, key, self.__ket__
            for line in self.unify(blocks[key]):
                print line
            if key != keys[-1]:
                print
    def replacement(self, old, new):
        if old in self.__order__:
            self.__order__[self.__order__.index(old)] = new
        else:
            print "WARNING: unknown target '%s' for replacement!" % old
    def new_block(self, name, after=None, before=None):
        if after is not None and after in self.__order__:
            i = self.__order__.index(after)
            end = self.__order__[i+1:]
            self.__order__ = self.__order__[:i+1]
            self.__order__.append(name)
            self.__order__.extend(end)
        elif before is not None and before in self.__order__:
            i = self.__order__.index(before)
            end = self.__order__[i:]
            self.__order__ = self.__order__[:i]
            self.__order__.append(name)
            self.__order__.extend(end)
        else:
            self.__order__.append(name)
    def get_sequence(self):
        atoms = self.group('atoms')
        sequence = []
        prev = None
        for atom in atoms:
            if atom['resid'] != prev:
                sequence.append(atom['resname'])
                prev = atom['resid']
        return sequence

class FastaIO(BlockIO):
    def __init__(self, filename=None, bra='>', ket=''):
        if filename is not None:
            self.__noread__ = False
            self.file = FileCrawler(filename)
        else:
            self.__noread__ = True
        self.__marks__ = {}
        self.__found_all__ = False
        self.__order__ = []
        self.__bra__ = str(bra)
        self.__ket__ = str(ket)
        self.comment = ['#', '%', ';']
    def __check__(self, line):
        if line and line[0] == self.__bra__:
            name = line[len(self.__bra__):].strip()
            if self.__marks__.has_key(name):
                if self.__marks__[name] != self.file.tell():
                    loki.warn("Duplicate index group '%s' found." % name)
            else:
                self.__marks__[name] = self.file.tell()
            return name
        return False
    def atomise(self, line):
        return [x for x in line.strip()]
    def unify(self, items):
        lines = []
        buffer = []
        for id in items:
            buffer.append(str(id))
            if len(buffer) == 80:
                lines.append(' '.join(buffer))
                buffer = []
        if len(buffer):
            lines.append(' '.join(buffer))
        return lines

class XVGIO:
    def __init__(self, filename=None):
        self.__noread__ = True
        self.__reset_point__ = 0
        self.__ro__ = re.compile('@ s\d+ legend "(.*)"')
        self.__found_all__ = False
        self.comments = []
        self.legends = []
        self.commands = []
        self.keys = []
        self.sets = []
        if filename is not None:
            self.open(filename)
    def __read_meta__(self):
        self.file.rewind()
        while True:
            try:
                line = self.__readline__()
            except EOF:
                break
            if not line: pass
            elif line[0] == '#':
                self.comments.append(line)
            elif line[0] == '@':
                if self.__ro__.match(line):
                    self.legends.append(self.__ro__.search(line).group(1))
                else:
                    self.commands.append(line)
            else:
                self.__size__ = len(line.split()) - 1
                self.__init_sets__()
                self.file.prevline()
                self.__reset_point__ = self.file.tell()
                self.__blocks__ = [self.__reset_point__]
                break
    def __init_sets__(self):
        self.sets = []
        while len(self.sets) < self.__size__:
            self.sets.append([])
    def __readline__(self):
        line = self.file.readline()
        if line == '':
            raise EOF
        return line.strip()
    def open(self, filename):
        if not self.__noread__:
            print "Reading new file '%s'." % filename
        self.file = FileCrawler(filename)
        self.__noread__ = False
        self.__read_meta__()
    def rewind(self):
        self.file.seek(self.__reset_point__)
    def read(self):
        if self.__noread__:
            raise IOError, 'no file open for reading'
        self.rewind()
        self.keys = []
        self.blocks = []
        self.__init_sets__()
        while True:
            try:
                keys, sets = self.read_block()
            except EOF:
                break
            if len(self.keys) and keys != self.keys:
                loki.warn('key mismatch between blocks.')
                loki.debug('self.keys='+repr(self.keys))
            self.blocks.append((keys, sets))
            if len(self.blocks) == 1:
                self.keys.extend(keys)
                for i in range(self.__size__):
                    self.sets[i].extend(sets[i])
        self.__found_all__ = True
        return self.keys, self.sets
    def read_block(self):
        keys = []
        sets = []
        for i in range(self.__size__):
            sets.append([])
        while True:
            try:
                line = self.__readline__()
            except EOF:
                if len(keys):
                    break
                raise EOF
            if not line:
                break
            loki.debug('line=' + repr(line))
            parts = line.split()
            key = parts.pop(0)
            loki.debug('key=' + repr(key))
            loki.debug('parts=' + repr(parts))
            loki.debug('__size__=' + repr(self.__size__))
            if len(parts) != self.__size__:
                if key == '&':
                    x = self.file.tell()
                    if x not in self.__blocks__:
                        self.__blocks__.append(x)
                    break
                else:
                    print 'WARNING: wrong number of columns... halt.'
                    sys.exit()
            keys.append(key)
            for i in range(len(parts)):
                sets[i].append(float(parts[i]))
        loki.debug('keys=' + repr(keys))
        loki.debug('sets=' + repr(sets))
        return keys, sets
    def write(self, filename=None):
        if filename is not None:
            old_stdout = sys.stdout
            sys.stdout = open(filename, 'w')
        self.echo_meta()
        i = 0
        for (keys, sets) in self.blocks:
            if i:
                print '&'
            self.echo_sets(keys, sets)
            i += 1
        if filename is not None:
            sys.stdout.close()
            sys.stdout = old_stdout
    def echo_meta(self):
        if len(self.comments):
            for line in self.comments:
                if line[0] != '#': 
                    line = '# ' + line
                print line
        if len(self.commands):
            for line in self.commands:
                print line
        if len(self.legends):
            # legend box turned on? FIXME: check-it properly
            if not len(self.commands):
                print '@ view 0.15, 0.15, 0.75, 0.85'
                print '@ legend on'
                print '@ legend box on'
                print '@ legend loctype view'
                print '@ legend 0.78, 0.8'
                print '@ legend length 2'
            for i in range(len(self.legends)):
                print '@ s%d legend "%s"' % (i, self.legends[i])
    def echo_sets(self, keys=None, sets=None):
        keys = keys or self.keys
        sets = sets or self.sets
        loki.debug('keys=' + repr(keys))
        loki.debug('sets=' + repr(sets))
        for i in range(len(keys)):
            data = []
            for j in range(len(sets)):
                data.append('%.6f' % (sets[j][i]))
            print '%s  %s' % (keys[i], '  '.join(data))

class GnuplotLogIO:
    def __init__(self, filename=None):
        self.__noread__ = True
        self.__eof__ = False
        self.__ro_file__ = re.compile(r'data read from "(\S+?)"')
        self.__ro_func__ = re.compile(
                r'function used for fitting: ([\w\d]+?)\(x\)')
        self.__ro_misc__ = re.compile(
                r'(\S[^:\n]*\S)\s+:\s+?([-+]?(?:\d+(?:\.\d*)?|\.\d+)' + \
                r'(?:[eE][-+]?\d+)?)')
        self.__misc_map__ = {
                'WSSR':'WSSR', 'delta(WSSR)':'delta', 'lambda':'lambda', 
                'delta(WSSR)/WSSR':'delta/WSSR', 
                'limit for stopping':'limit', 
                'final sum of squares of residuals':'sum-squared', 
                'rel. change during last iteration':'change', 
                'degrees of freedom (ndf)':'ndf', 
                'rms of residuals      (stdfit) = sqrt(WSSR/ndf)':'stdfit', 
                'variance of residuals (reduced chisquare) = WSSR/ndf':'chisquare'
                }
        self.__ro_snip__ = re.compile(
                r'Asymptotic Standard Error[=\n\s]+([\S\s]*?)' + \
                r'correlation matrix')
        self.__ro_para__ = re.compile(
                r'([\w\d]+)\s+=\s+?([-+]?(?:\d+(?:\.\d*)?|\.\d+)' + \
                r'(?:[eE][-+]?\d+)?)\s+\+\/\-\s' + \
                r'([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)')
        self.__ro_mark__ = re.compile(r'\*{79}')
        if filename is not None:
            self.open(filename)
    def __readline__(self):
        line = self.file.readline()
        if line == '':
            self.__eof__ = True
            raise EOF
        return line
    def rewind(self):
        self.file.seek(0)
        self.__eof__ = False
    def next_block(self):
        count = 0
        while count < 2:
            line = self.__readline__()
            if self.__ro_mark__.match(line):
                count += 1
        self.file.prevline()
    def __read__(self):
        if self.__eof__:
            raise EOF
        count = 0
        self.__buffer__ = ''
        while count < 2:
            try:
                line = self.__readline__()
            except EOF:
                return self.__buffer__
            else:
                if self.__ro_mark__.match(line):
                    count += 1
                elif count:
                    self.__buffer__ += line
        self.file.prevline()
        return self.__buffer__
    def __analyse__(self, txt=None):
        txt = txt or self.__buffer__
        if not txt:
            raise ValueError, 'No text to analyse.'
        try:
            snippet = self.__ro_snip__.search(txt).group(1)
            block = {
                    'filename': self.__ro_file__.search(txt).group(1),
                    'function': self.__ro_func__.search(txt).group(1),
                    'parameters': self.__ro_para__.findall(snippet)
                    }
            misc = self.__ro_misc__.findall(txt)
            for item in misc:
                if item[0] in self.__misc_map__.keys():
                    block[self.__misc_map__[item[0]]] = item[1]
                else:
                    print "Unknown item '%s'." % item[0]
        except AttributeError:
            print 'Unknown syntax. Is the text really from a Gnuplot log?'
            return {}
        return block
    def read_block(self):
        self.__read__()
        return self.__analyse__()
    def open(self, filename):
        if not self.__noread__:
            print "Reading new file '%s'." % filename
        self.file = FileCrawler(filename)
        self.__noread__ = False

class HoleIO:
    def __init__(self, filename=None):
        self.__noread__ = True
        self.__reset_point__ = 0
        self.__end_point__ = None
        self.__ro__ = re.compile('Minimum radius found:\s+(\d*\.\d+)\s+angstroms.')
        self.keys = []
        self.data = {'radius': [], 'distance': [], 'thickness': []}
        if filename is not None:
            self.open(filename)
    def __read_meta__(self):
        self.file.rewind()
        while True:
            before = self.file.tell()
            try:
                line = self.__readline__()
            except EOF:
                break
            if not line: pass
            elif line == 'cenxyz.cvec      radius  cen_line_D sum{s/(area point sourc':
                self.__reset_point__ = self.file.tell()
            elif self.__ro__.match(line):
                self.minimum = self.__ro__.search(line).group(1)
                end_point = before
            else: pass # FIXME: catch metadata
        self.__end_point__ = end_point
    def __readline__(self):
        if self.file.tell() == self.__end_point__:
            raise EOF
        line = self.file.readline()
        if line == '':
            raise EOF
        return line.strip()
    def open(self, filename):
        if not self.__noread__:
            print "Reading new file '%s'." % filename
        self.file = FileCrawler(filename)
        self.__noread__ = False
        self.__read_meta__()
    def rewind(self):
        self.file.seek(self.__reset_point__)
    def read(self, key=None):
        if self.__noread__:
            raise IOError, 'no file open for reading'
        self.rewind()
        self.keys = []
        self.data = {'radius': [], 'distance': [], 'thickness': []}
        while True:
            try:
                line = self.__readline__()
            except EOF:
                break
            if not line:
                break
            parts = line.split()
            self.keys.append(parts[0])
            self.data['radius'].append(float(parts[1]))
            self.data['distance'].append(float(parts[2]))
            self.data['thickness'].append(float(parts[3]))
        if key is not None:
            return self.keys, self.data[key]
        return self.keys, self.data


# choose correct components and clear namespace
if not __biopython_installed__:
    PDB = rawPDB
if __name__ == '__main__':
    if __biopython_installed__:
        print 'Using BioPython for parsing PDBs.'
    else:
        print 'Using internal PDB parsers.'


