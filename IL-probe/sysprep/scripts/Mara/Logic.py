#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.Logic (v. 0.7):
    Simple wrapper functions for chaining together complex logical 
    expressions.

    Functions with self-evident names include:
      LessThan
      GreaterThan
      Equal
      And
      Or
      
    Requirements: Python 2.2->
    
    Author: Martti Louhivuori (martti.louhivuori@helsinki.fi), 23.12.2006

    Date: 23.12.2006

    ---

    Copyright (C) 2006  Martti Louhivuori

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    The full license text can be found in the file LICENSE.
"""

class LessThan:
    def __init__(self, x):
        self.__x__ = x
    def __call__(self, y):
        return y < self.__x__
    def __repr__(self):
        return '< ' + str(self.__x__)

class GreaterThan:
    def __init__(self, x):
        self.__x__ = x
    def __call__(self, y):
        return y > self.__x__
    def __repr__(self):
        return '> ' + str(self.__x__)

class Equal:
    def __init__(self, x):
        self.__x__ = x
    def __call__(self, y):
        return y == self.__x__
    def __repr__(self):
        return '== ' + str(self.__x__)

class And:
    def __init__(self, rules):
        self.__rules__ = rules
    def __call__(self, y):
        for rule in self.__rules__:
            if not rule(y):
                return False
        return True
    def __repr__(self):
        return ' && '.join([repr(x) for x in self.__rules__])

class Or:
    def __init__(self, rules):
        self.__rules__ = rules
    def __call__(self, y):
        for rule in self.__rules__:
            if rule(y):
                return True
        return False
    def __repr__(self):
        return ' || '.join([repr(x) for x in self.__rules__])

