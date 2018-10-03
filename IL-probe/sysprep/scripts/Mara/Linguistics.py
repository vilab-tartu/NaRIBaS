#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.Linguistics (v. 0.8):
    Translation routines that provide an easy way to switch from one naming
    convention to another provided they are both known.
    
    Requirements: Python 2.2->

    TODO: - add documentation

    Author: Martti Louhivuori (martti.louhivuori@helsinki.fi), 6.10.2005

    Date: 14.2.2006

    ---

    Copyright (C) 2006  Martti Louhivuori

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    The full license text can be found in the file LICENSE.
"""
# TODO: - add synonyms to replications

from Mara import loki
from Mara.aux import quicksort, ArgumentError
from copy import deepcopy

class Translator:
    """
    Translates words from one language to another to its best ability, i.e.
    uses knowledge of third languages as well as any direct connection between
    the two. 
    In other words: A general purpose straight-forward word substitutor. ;)
    """

    def __init__(self, dictionaries, codeces, default=None, unknown=None):
        """
        Creates a Translator instance.

        arguments:
          dictionaries -- a list of Dictionaries
          codeces      -- a list of Codeces
          default      -- default language (default: None)
          unknown      -- how to denote an unknown word (default: None)
        """
        if False in [isinstance(d, Dictionary) for d in dictionaries]:
            raise ArgumentError, 'Unknown dictionary format encountered.'
        if False in [isinstance(c, Codex) for c in codeces]:
            raise ArgumentError, 'Unknown codex format encountered.'
        self.__dictionaries__ = {}
        for d in dictionaries:
            self.__dictionaries__[str(d)] = d
        self.__codeces__ = codeces
        self.__init_paths__()
        self.__init_links__()
        self.default = default
        self.unknown = unknown
    
    def __contains__(self, x):
        for d in self.__dictionaries__.values():
            if x in d:
                return True
        return False

    def __str__(self):
        s = 'Dictionaries: '
        tab = '\n'+' '*len(s)
        s += tab.join([str(d) for d in self.__dictionaries__])
        s2 = 'Codeces: '
        tab = '\n'+' '*len(s2)
        s2 += tab.join([str(c) for c in self.__codeces__])
        return s + '\n' + s2

    def __repr__(self):
        return 'Translator(%s, %s, %s, %s)'%(repr(self.__dictionaries__), \
                repr(self.__codeces__), repr(self.default), \
                repr(self.unknown))

    def __match_attributes__(self, x, y):
        matches, mismatches = 0, 0
        for attr in x:
            if attr in y:
                matches += 1
            else:
                mismatches += 1
        mismatches += len(y)-matches
        return matches, mismatches

    def __call__(self, word, target=None, attributes=None, find_best=True):
        """
        Translate a word into a target language matching desired attributes.

        arguments:
          word       -- word of interest
          target     -- target language (default: None == self.default)
          attributes -- list of desired attributes (default: [])
          find_best  -- whether to actually search for a best match
                        (default: True)

        TODO: find_best redundant, remove print statements etc.
        """
        if target == None:
            target = self.default
        if attributes is None:
            attributes = []
        elif type(attributes) is str:
            attributes = [attributes]
        loki.debug('attributes=%s', attributes)
        loki.debug('word=%s' % repr(word))
        if not self.__dictionaries__.has_key(target):
            raise ArgumentError, "Unknown target language '%s'." % target
        if not self.__links__.has_key(target):
            raise ArgumentError, "No known links to target language '%s'." \
                    % target
        if word in self.__dictionaries__[target]:
            return word
        else:
            name = None
            original = None
            for dictionary in self.__dictionaries__.values():
                if word in dictionary:
                    name = str(dictionary)
                    original = dictionary[word]
                    for a in original.attributes:
                        if a not in attributes:
                            attributes.append(a)
                    break
            loki.debug('original=%s' % repr(original))
            loki.debug('attributes=%s', attributes)
            loki.debug('name=%s' % repr(name))
            if name != None:
                path = self.__paths__[name][target]
                loki.debug('path=%s' % repr(path))
                if path != None:
                    # start from the base form
                    if original.base:
                        word = original.base
                    loki.debug('word=%s' % repr(word))
                    # traverse through the links to find the target word
                    for p in path:
                        word = str(self.__links__[name][p][word])
                        loki.debug('word=%s @ %s' % (repr(word), repr(p)))
                        name = p
                    # find the best match for the original form
                    if find_best:
                        new = self.__dictionaries__[target][word]
                        loki.debug('new=%s' % repr(new))
                        matches, mismatches = self.__match_attributes__(
                                attributes, new.attributes)
                        loki.debug('matches=%s', matches)
                        loki.debug('mismatches=%s', mismatches)
                        if mismatches:
                            alt = [word] + new.synonyms + new.conjugates
                            loki.debug('alt=%s' % repr(alt))
                            record = []
                            best = None
                            for a in alt:
                                m, mm = self.__match_attributes__( \
                                        attributes, \
                                        self.__dictionaries__\
                                                [target][a].attributes)
                                record.append((m,mm))
                                if best == None or m > record[best][0] or \
                                        (m == record[best][0] and \
                                        mm < record[best][1]):
                                    best = len(record)-1
                            loki.debug('record=%s', record)
                            loki.debug('best=%s', best)
                            if record[best][0] > matches or \
                                    (record[best][0] == matches and \
                                    record[best][1] < mismatches):
                                return alt[best]
                    return word
            return self.unknown

    def __init_paths__(self):
        pairs = []
        for c in self.__codeces__:
            pairs.append(c.get_languages())
        languages = []
        for p in pairs:
            if p[0] not in languages:
                languages.append(p[0])
            if p[1] not in languages:
                languages.append(p[1])
        # init distance table
        row = {}
        for l in languages:
            row[l] = None
        matrix = {}
        for l in languages:
            matrix[l] = deepcopy(row)
            matrix[l][l] = []
        # neighbours
        neighbours = {}
        for pair in pairs:
            matrix[pair[0]][pair[1]] = [pair[1]]
            matrix[pair[1]][pair[0]] = [pair[0]]
            neighbours[pair[0]] = neighbours.get(pair[0],[]) + [pair[1]]
            neighbours[pair[1]] = neighbours.get(pair[1],[]) + [pair[0]]
        # rest
        unconnected = 0
        for key1 in matrix.keys():
            for key2 in matrix[key1].keys():
                if matrix[key1][key2] == None:
                    unconnected += 1
        level = 2
        while level < len(languages) and unconnected:
            for key1 in matrix.keys():
                for key2 in matrix[key1].keys():
                    if matrix[key1][key2] != None and \
                            len(matrix[key1][key2]) == level-1:
                        for n in neighbours[key1]:
                            if matrix[key2][n] == None:
                                matrix[key2][n] = matrix[key2][key1] + [n]
                                matrix[n][key2] = [key1] + matrix[key1][key2]
                                unconnected -= 2
                        for n in neighbours[key2]:
                            if matrix[key1][n] == None:
                                matrix[key1][n] = matrix[key1][key2] + [n]
                                matrix[n][key1] = [key2] + matrix[key2][key1]
                                unconnected -= 2
            level += 1
        self.__paths__ = matrix

    def __init_links__(self):
        self.__links__ = {}
        for i in range(len(self.__codeces__)):
            x,y = self.__codeces__[i].get_languages()
            q = self.__links__.get(x,{})
            q[y] = self.__codeces__[i]
            self.__links__[x] = q
            q = self.__links__.get(y,{})
            q[x] = self.__codeces__[i]
            self.__links__[y] = q

    def base(self, word, language=None):
        """
        Return the base form of a word.

        arguments:
          word     -- word of interest
          language -- the language of the word if known (default: None)
        """
        if not language:
            for dictionary in self.__dictionaries__.values():
                if word in dictionary:
                    return dictionary[word].base or word
        elif language in self.__dictionaries__:
            w = self.__dictionaries__[language][word]
            if hasattr(w, 'base'):
                return w.base or word
            else:
                return None
        else:
            return None

    def recognise_language(self, language):
        return language in self.__dictionaries__.keys()
    def recognise_word(self, word):
        return word in self
    def language(self, word):
        """
        Return the language of a word.

        arguments:
          word -- word of interest
        """
        for l,d in self.__dictionaries__.items():
            if word in d:
                return l
        return None

class CodexError(Exception):
    def __init__(self): pass
    def __repr__(self):
        print 'The source file for codex must contain either a pickled or a'
        print "shelved alphabet. If shelved the corresponding 'key' must be"
        print 'provided.'

class AmbiguousSynonym:
    def __init__(self, words):
        if type(words) == str:
            self.__words__ = [words]
        else:
            self.__words__ = list(words)
    def __str__(self):
        return self.__words__[0]
    def __repr__(self):
        return 'AmbiguousSynonym(%s)'%repr(self.__words__)
    def __len__(self):
        return len(self.__words__)
    def __contains__(self, word):
        return self.__words__.__contains__(word)

class Codex:
    """
    Connections between words in different languages. May be used as a
    simplistic code of translation. NOTE: real languages are so complex 
    that a comprehensive one-to-one mapping between two languages is 
    next to impossible.
    """

    def __init__(self, languages, source=None, key=None, unknown=None):
        if source == None:
            self.__alphabet__ = {}
        elif type(source) == dict:
            self.__alphabet__ = source
        elif type(source) != str:
            if key != None:
                try:
                    dict = shelve.open(source)
                    self.__alphabet__ = dict[key]
                except KeyError:
                    raise CodexError
            else:
                self.__alphabet__ = pickle.load(source)
        else:
            raise CodexError
        self.__languages__ = languages

    def __setitem__(self, key, value):
        self.__alphabet__[key] = AmbiguousSynonym(value)

    def __getitem__(self, key):
        return self.__alphabet__[key]

    def __str__(self):
        return repr(self.__languages__)

    def __repr__(self):
        return 'Codex(%s)'%repr(self.__alphabet__)
    
    def __len__(self):
        return len(self.__alphabet__)
    
    def __contains__(self, key):
        return self.__alphabet__.has_key(key)

    def associate(self, x, y):
        if type(x) == str:
            self[x] = y
        else:
            for z in x:
                self[z] = y
        if type(y) == str:
            self[y] = x
        else:
            for z in y:
                self[z] = x

    def store(self, file, key=None):
        if key != None:
            f = open(file, 'w')
            pickle.dump(self.__alphabet__, f)
            f.close()
        else:
            f = shelve.open(file)
            f[key] = self.__alphabet__
            f.close()

    def get_languages(self):
        return tuple(self.__languages__)


class Word:
    def __init__(self, word, description='', base=None, attributes=[], \
            synonyms=[], conjugates=[]):
        self.__word__ = str(word)
        self.description = str(description)
        self.base = base
        if type(synonyms) in [tuple, list]:
            self.synonyms = list(synonyms)
        else:
            self.synonyms = [str(synonyms)]
        if type(conjugates) in [tuple, list]:
            self.conjugates = list(conjugates)
        else:
            self.conjugates = [str(conjugates)]
        if type(attributes) is str:
            self.attributes = [attributes]
        else:
            self.attributes = list(attributes)
    def __str__(self):
        return self.__word__
    def __repr__(self):
        return "Word('%s', '%s', base='%s', attributes='%s')"%(
                self.__word__, self.description, self.base,
                repr(self.attributes))
    def __lt__(self, other):
        return self.__word__ < other
    def __le__(self, other):
        return self.__word__ <= other
    def __eq__(self, other):
        return self.__word__ == other
    def __ne__(self, other):
        return self.__word__ != other
    def __gt__(self, other):
        return self.__word__ > other
    def __ge__(self, other):
        return self.__word__ >= other
    def replica(self, word):
        return Word(word, self.description, base=self.base, \
                attributes=self.attributes, synonyms=self.synonyms, \
                conjugates=self.conjugates)


class Dictionary:
    def __init__(self, language, words=[]):
        self.__language__ = str(language)
        self.__index__ = {}
        for w in words:
            self.__index__[str(w)] = w
    def __str__(self):
        return self.__language__
    def __repr__(self):
        return "Dictionary('%s', %s)"%(self.__language__, \
                repr(self.__index__.values()))
    def __contains__(self, x):
        return self.__index__.has_key(x)
    def __setitem__(self, key, value):
        self.__index__[key] = value
    def __getitem__(self, key):
        try:
            return self.__index__[key]
        except KeyError:
            loki.warn("Unknown word '%s'." % str(key))
            return None
    def sort(self, what='words', direction='ascending'):
        if what == 'words':
            sorted = self.__index__.keys()
            quicksort(sorted)
        else:
            descs = []
            flag = '<<<DICTIONARY-TEMP-FLAG>>>'
            for k in self.__index__.keys():
                descs.append(self.__index__[k].description + flag + \
                        str(self.__index__[k]))
            quicksort(descs)
            sorted = []
            for desc in descs:
                d,w = desc.split(flag)
                sorted.append(w)
        if direction != 'ascending':
            sorted.reverse()
        return sorted
    def set_synonyms(self, words):
        for i in range(len(words)):
            for j in range(i+1,len(words)):
                self.__index__[words[i]].synonyms.append(words[j])
                self.__index__[words[j]].synonyms.append(words[i])
    def parse_conjugates(self):
        for k in self.__index__.keys():
            w = self.__index__[k]
            if w.base:
                self.__index__[str(w.base)].conjugates.append(str(w))
    def replicate(self, old, new):
        word = self[old].replica(new)
        self[new] = word

