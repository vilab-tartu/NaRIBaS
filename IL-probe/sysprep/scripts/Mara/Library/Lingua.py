#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.Library.Lingua (v. 0.8):
    A collection of dictionaries & translation codeces for natural
    languages as well as symbol sets.
    NOTE: Mara.Library provides general purpose translators with no
          custom tweaking needed.

    Requirements: Python 2.2->

    Author: Martti Louhivuori (martti.louhivuori@helsinki.fi), 8.2.2006

    Date: 16.3.2009

    ---

    Copyright (C) 2006  Martti Louhivuori

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    The full license text can be found in the file LICENSE.
"""
# TODO: - add synonyms to replications (maybe in Linguistics.py)

from Mara.Linguistics import * # Dictionary, Word, Codex
from Mara.Library import CommonLib

"""
-----------------------------------------------------------------------------
Linguistic section. Dictionaries, codices etc.
-----------------------------------------------------------------------------
"""

libLingua = CommonLib('libLingua: Linguistic library')

#############################################################################
# DICTIONARY DATABASE                                                       #
#############################################################################

"""
The one letter names of aminoacids.
"""
libLingua.dictionaryAmino1 = Dictionary('1-letter-aminoacids', 
        [
        Word('A', 'Alanine'),
        Word('R', 'Arginine'),
        Word('N', 'Asparganine'),
        Word('D', 'Aspartate'),
        Word('C', 'Cysteine'),
        Word('E', 'Glutamate'),
        Word('Q', 'Glutamine'),
        Word('G', 'Glycine'),
        Word('H', 'Histidine'),
        Word('I', 'Isoleucine'),
        Word('L', 'Leucine'),
        Word('K', 'Lysine'),
        Word('M', 'Methionine'),
        Word('F', 'Phenylalanine'),
        Word('P', 'Proline'),
        Word('S', 'Serine'),
        Word('T', 'Threonine'),
        Word('W', 'Tryptophan'),
        Word('Y', 'Tyrosine'),
        Word('V', 'Valine'),
        Word('R+', 'Arginine (charged)', 'R', 'charged'),
        Word('D-', 'Aspartate (charged)', 'D', 'charged'),
        Word('E-', 'Glutamate (charged)', 'E', 'charged'),
        Word('H+', 'Histidine (charged)', 'H', 'charged'),
        Word('K+', 'Lysine (charged)', 'K', 'charged'),
        Word('Ms', 'Selenomethionine', 'M', 'selenium')
        ])
libLingua.dictionaryAmino1.parse_conjugates()

"""
The three letter names of aminoacids.
"""
# FIX HISTIDINES ETC. !!!
libLingua.dictionaryAmino3 = Dictionary('3-letter-aminoacids', 
        [
        Word('Ala', 'Alanine'),
        Word('Arg', 'Arginine'),
        Word('Asn', 'Asparganine'),
        Word('Asp', 'Aspartate'),
        Word('Cys', 'Cysteine'),
        Word('Glu', 'Glutamate'),
        Word('Gln', 'Glutamine'),
        Word('Gly', 'Glycine'),
        Word('His', 'Histidine'),
        Word('Ile', 'Isoleucine'),
        Word('Leu', 'Leucine'),
        Word('Lys', 'Lysine'),
        Word('Met', 'Methionine'),
        Word('Phe', 'Phenylalanine'),
        Word('Pro', 'Proline'),
        Word('Ser', 'Serine'),
        Word('Thr', 'Threonine'),
        Word('Trp', 'Tryptophan'),
        Word('Tyr', 'Tyrosine'),
        Word('Val', 'Valine'),
        Word('Arg+', 'Arginine (charged)', 'Arg', 'charged'),
        Word('Asp-', 'Aspartate (charged)', 'Asp', 'charged'),
        Word('Cyss', 'Cysteine (free sulphur)', 'Cys', 'free sulphur'),
        Word('Glu-', 'Glutamate (charged)', 'Glu', 'charged'),
        Word('His+', 'Histidine (charged)', 'His', 'charged'),
        Word('Hist', 'Histidine (alt.)', 'His', 'alt'),
        Word('Hist+', 'Histidine (alt., charged)', 'His', ['alt','charged']),
        Word('Lys+', 'Lysine (charged)', 'Lys', 'charged'),
        Word('Mse', 'Selenomethionine', 'Met', 'selenium')
        ])
libLingua.dictionaryAmino3.parse_conjugates()

"""
Notations that are in general use among the NMR community. TODO: add more.
"""
libLingua.dictionaryNMR = Dictionary('general NMR',
        [Word('HACA', 'D- / J-coupling btw alpha carbon and alpha proton.'),
        Word('CAHA', 'D- / J-coupling btw alpha carbon and alpha proton.', 
            'HACA', 'proton last'),
        Word('CACB', 'D- / J-coupling btw alpha carbon and beta carbon.'),
        Word('CBCA', 'D- / J-coupling btw alpha carbon and beta carbon.', 
            'CACB', 'backwards'),
        Word('HACO', 'D- / J-coupling btw alpha proton and carboxyl carbon.'),
        Word('COHA', 'D- / J-coupling btw alpha proton and carboxyl carbon.', 
            'HACO', 'proton last'),
        Word('CACO', 'D- / J-coupling btw alpha carbon and carboxyl carbon.'),
        Word('COCA', 'D- / J-coupling btw alpha carbon and carboxyl carbon.', 
            'CACO', 'backwards'),
        Word('HAN', 'D- / J-coupling btw nitrogen and alpha proton.'),
        Word('NHA', 'D- / J-coupling btw nitrogen and alpha proton.', 
            'HAN', 'proton last'),
        Word('NCA', 'D- / J-coupling btw nitrogen and alpha carbon.'),
        Word('CAN', 'D- / J-coupling btw nitrogen and alpha carbon.',
            'NCA', 'backwards'),
        Word('HN', 'D- / J-coupling btw nitrogen and nitrogen proton.'),
        Word('NH', 'D- / J-coupling btw nitrogen and nitrogen proton.', 
            'HN', 'proton last'),
        Word('CON', 'D- / J-coupling btw nitrogen and carboxyl carbon.'),
        Word('NCO', 'D- / J-coupling btw nitrogen and carboxyl carbon.', 
            'CON', 'backwards'),
        Word('NCO[-1]', 'D- / J-coupling btw nitrogen and carboxyl carbon.', 
            'CON', ['backwards','index']),
        Word('HCO', 'D- / J-coupling btw carboxyl carbon and nitrogen proton.'),
        Word('COH', 'D- / J-coupling btw carboxyl carbon and nitrogen proton.', 
            'HCO', 'proton last'),
        Word('HCO[-1]', 'D- / J-coupling btw carboxyl carbon and nitrogen proton.', 
            'HCO', 'index'),
        Word('HNCO', 'D- / J-coupling btw carboxyl carbon and nitrogen proton.', 
            'HCO', 'explicit'),
        Word('HNCO[-1]', 'D- / J-coupling btw carboxyl carbon and nitrogen proton.', 
            'HCO', ['explicit','index'])])
libLingua.dictionaryNMR.parse_conjugates()

"""
English dictionary. Huh huh. TODO: use some open dictionary as base, write
maintenance routines, move actual data out of here and only read it here
"""
libLingua.dictionaryEnglish = Dictionary('english',
        [
        Word('alanine','A simple nonessential crystalline amino acid. Chemical structure: CH3-CH(NH2)-COOH'),
        Word('arginine','A crystalline basic amino acid derived from guanidine. Chemical structure: HN=C(NH2)-NH-(CH2)3-CH(NH2)-COOH'),
        Word('asparagine','A nonessential amino acid that is an amide of aspartic acid. Chemical structure: H2N-CO-CH2-CH(NH2)-COOH'),
        Word('aspartate','Syn. aspartic acid. Aminoacid / a salt or ester of aspartic acid. Chemical structure: HOOC-CH2-CH(NH2)-COOH', 
            synonyms='aspartic acid'),
        Word('cysteine','A crystalline sulfur-containing amino acid readily oxidizable to cystine. Chemical structure: HS-CH2-CH(NH2)-COOH'),
        Word('glutamine','A crystalline amino acid C5H10N2O3 that is found both free and in proteins in plants and animals and that yields glutamic acid and ammonia on hydrolysis. Chemical structure: H2N-CO-(CH2)2-CH(NH2)-COOH'),
        Word('glutamate','Syn. glutamic acid. Aminoacid / a salt or ester of glutamic acid. Chemical structure: HOOC-(CH2)2-CH(NH2)-COOH', 
            synonyms='glutamic acid'),
        Word('glycine','A sweet crystalline amino acid obtained especially by hydrolysis of proteins. Chemical structure: NH2-CH2-COOH.'),
        Word('histidine','A crystalline essential amino acid formed by the hydrolysis of most proteins. Chemical structure: N*H-CH=N-CH=C*-CH2-CH(NH2)-COOH (* = connected)'),
        Word('isoleucine','A crystalline essential amino acid isomeric with leucine. Chemical structure: CH3-CH2-CH(CH3)-CH(NH2)-COOH'),
        Word('leucine','A white crystalline essential amino acid obtained by the hydrolysis of most dietary proteins. Chemical structure: (CH3)2-CH-CH2-CH(NH2)-COOH'),
        Word('lysine','A crystalline essential amino acid obtained from the hydrolysis of various proteins. Chemical structure: H2N-(CH2)4-CH(NH2)-COOH'),
        Word('methionine','A crystalline sulfur-containing essential amino acid. Chemical structure: CH3-S-(CH2)2-CH(NH2)-COOH'),
        Word('phenylalanine','An essential amino acid that is converted in the normal body to tyrosine. Chemical structure: Ph-CH2-CH(NH2)-COOH'),
        Word('proline','An amino acid that can be synthesized by animals from glutamate. Chemical structure: N*H-(CH2)3-C*H-COOH (* = connected)'),
        Word('serine','A crystalline nonessential amino acid that occurs especially as a structural part of many proteins. Chemical structure: HO-CH2-CH(NH2)-COOH'),
        Word('threonine','A colorless crystalline essential amino acid. Chemical structure: CH3-CH(OH)-CH(NH2)-COOH'),
        Word('tryptophan','A crystalline essential amino acid that is widely distributed in proteins. Chemical structure: Ph*-NH-CH=C*-CH2-CH(NH2)-COOH (* = connected)'),
        Word('tyrosine','A phenolic amino acid that is a precursor of several important substances (as epinephrine and melanin). Chemical structure: HO-p-Ph-CH2-CH(NH2)-COOH'),
        Word('valine','A crystalline essential amino acid that is one of the building blocks of plant and animal proteins. Chemical structure: (CH3)2-CH-CH(NH2)-COOH'),
        Word('selenomethionine','Selenomethionine is a methionine containing selenium.')
        ])

"""
The standard Mara notations. TODO: add more.
"""
libLingua.dictionaryMara = Dictionary('mara',
    [Word('DHaCa', 'D- / J-coupling btw alpha carbon and alpha proton.'),
    Word('DCaHa', 'D- / J-coupling btw alpha carbon and alpha proton.',
        'DHaCa', 'proton last'),
    Word('DCaCb', 'D- / J-coupling btw alpha carbon and beta carbon.'),
    Word('DHaCo', 'D- / J-coupling btw alpha proton and carboxyl carbon.'),
    Word('DCoHa', 'D- / J-coupling btw alpha proton and carboxyl carbon.',
        'DHaCo', 'proton last'),
    Word('DCaCo', 'D- / J-coupling btw alpha carbon and carboxyl carbon.'),
    Word('DHaN', 'D- / J-coupling btw nitrogen and alpha proton.'),
    Word('DNHa', 'D- / J-coupling btw nitrogen and alpha proton.',
        'DHaN', 'proton last'),
    Word('DNCa', 'D- / J-coupling btw nitrogen and alpha carbon.'),
    Word('DHnN', 'D- / J-coupling btw nitrogen proton and nitrogen.'),
    Word('DNHn', 'D- / J-coupling btw nitrogen proton and nitrogen.',
        'DHnN', 'proton last'),
    Word('DCoN', 'D- / J-coupling btw carboxyl carbon and nitrogen.'),
    Word('DHnCo', 'D- / J-coupling btw carboxyl carbon and nitrogen proton.'),
    Word('DCoHn', 'D- / J-coupling btw carboxyl carbon and nitrogen proton.',
        'DHnCo', 'proton last')])
libLingua.dictionaryMara.parse_conjugates()

"""
Atom names. Default
"""
libLingua.dictionaryAtom = Dictionary('atom',
        [Word('H', 'hydrogen'),
        Word('He', 'helium'),
        Word('B', 'boron'),
        Word('C', 'carbon'),
        Word('N', 'nitrogen'),
        Word('O', 'oxygen'),
        Word('F', 'phosphor'),
        Word('P', 'led'),
        Word('S', 'sulphur'),
        Word('SE', 'selenium'),
        Word('CA', 'alpha carbon', 'C', 'alpha'),
        Word('CB', 'beta carbon', 'C', 'beta'),
        Word('CG', 'gamma carbon', 'C', 'gamma'),
        Word('CG1', 'gamma carbon (1st)', 'C', ['gamma', '1st']),
        Word('CG2', 'gamma carbon (2nd)', 'C', ['gamma', '2nd']),
        Word('CD', 'delta carbon', 'C', 'delta'),
        Word('CD1', 'delta carbon (1st)', 'C', ['delta', '1st']),
        Word('CD2', 'delta carbon (2nd)', 'C', ['delta', '2nd']),
        Word('CE', 'epsilon carbon', 'C', 'epsilon'),
        Word('CE1', 'epsilon carbon (1st)', 'C', ['epsilon', '1st']),
        Word('CE2', 'epsilon carbon (2nd)', 'C', ['epsilon', '2nd']),
        Word('CE3', 'epsilon carbon (3rd)', 'C', ['epsilon', '3rd']),
        Word('CZ', 'zeta carbon', 'C', 'zeta'),
        Word('CZ2', 'zeta carbon (2nd)', 'C', ['zeta', '2nd']),
        Word('CZ3', 'zeta carbon (3rd)', 'C', ['zeta', '3rd']),
        Word('CH1', 'eta carbon (1st)', 'C', ['eta', '1st']),
        Word('CH2', 'eta carbon (2nd)', 'C', ['eta', '2nd']),
        Word('CH3', 'eta carbon (3rd)', 'C', ['eta', '3rd']),
        Word('C1', 'carbon (1st)', 'C', '1st'),
        Word('C2', 'carbon (2nd)', 'C', '2nd'),
        Word('HA', 'alpha proton', 'H', 'alpha'),
        Word('HB', 'beta proton', 'H', 'beta'),
        Word('HG', 'gamma proton', 'H', 'gamma'),
        Word('HG1', 'gamma proton (1st)', 'H', ['gamma', '1st']),
        Word('HG2', 'gamma proton (2nd)', 'H', ['gamma', '2nd']),
        Word('HD', 'delta proton', 'H', 'delta'),
        Word('HD1', 'delta proton (1st)', 'H', ['delta', '1st']),
        Word('HD2', 'delta proton (2nd)', 'H', ['delta', '2nd']),
        Word('HE', 'epsilon proton', 'H', 'epsilon'),
        Word('HE1', 'epsilon proton (1st)', 'H', ['epsilon', '1st']),
        Word('HE2', 'epsilon proton (2nd)', 'H', ['epsilon', '2nd']),
        Word('HZ', 'zeta proton', 'H', 'zeta'),
        Word('HH', 'eta proton', 'H', 'eta'),
        Word('HH1', 'eta proton (1st)', 'H', ['eta', '1st']),
        Word('HH2', 'eta proton (2nd)', 'H', ['eta', '2nd']),
        Word('HN', 'proton (amide)', 'H', 'amide'),
        Word('HT1', 'terminal proton (1st)', 'H', ['terminus', '1st']),
        Word('HT2', 'terminal proton (2nd)', 'H', ['terminus', '2nd']),
        Word('HT3', 'terminal proton (3rd)', 'H', ['terminus', '3rd']),
        Word('ND', 'delta nitrogen', 'N', 'delta'),
        Word('ND1', 'delta nitrogen (1st)', 'N', ['delta', '1st']),
        Word('ND2', 'delta nitrogen (2nd)', 'N', ['delta', '2nd']),
        Word('NE', 'epsilon nitrogen', 'N', 'epsilon'),
        Word('NE1', 'epsilon nitrogen (1st)', 'N', ['epsilon', '1st']),
        Word('NE2', 'epsilon nitrogen (2nd)', 'N', ['epsilon', '2nd']),
        Word('NZ', 'zeta nitrogen', 'N', 'zeta'),
        Word('NH', 'eta nitrogen', 'N', 'eta'),
        Word('NH1', 'eta nitrogen (1st)', 'N', ['eta', '1st']),
        Word('NH2', 'eta nitrogen (2nd)', 'N', ['eta', '2nd']),
        Word('OG', 'gamma oxygen', 'O', 'gamma'),
        Word('OG1', 'gamma oxygen (1st)', 'O', ['gamma', '1st']),
        Word('OG2', 'gamma oxygen (2nd)', 'O', ['gamma', '2nd']),
        Word('OD', 'delta oxygen', 'O', 'delta'),
        Word('OD1', 'delta oxygen (1st)', 'O', ['delta', '1st']),
        Word('OD2', 'delta oxygen (2nd)', 'O', ['delta', '2nd']),
        Word('OE', 'epsilon oxygen', 'O', 'epsilon'),
        Word('OE1', 'epsilon oxygen (1st)', 'O', ['epsilon', '1st']),
        Word('OE2', 'epsilon oxygen (2nd)', 'O', ['epsilon', '2nd']),
        Word('OH', 'eta oxygen / hydroxyl', 'O', 'epsilon'),
        Word('OT1', 'terminal oxygen (1st)', 'O', ['terminus', '1st']),
        Word('OT2', 'terminal oxygen (2nd)', 'O', ['terminus', '2nd']),
        Word('OXT', 'oxygen (?)', 'O', 'xt'),
        Word('O1P', 'oxygen (?)', 'O', '1p'),
        Word('O2P', 'oxygen (?)', 'O', '2p'),
        Word('O3P', 'oxygen (?)', 'O', '3p'),
        Word('O4P', 'oxygen (?)', 'O', '4p'),
        Word('O1', 'oxygen (1st)', 'O', '1st'),
        Word('O2', 'oxygen (2nd)', 'O', '2nd'),
        Word('SD', 'delta sulphur', 'S', 'delta'),
        Word('SG', 'gamma sulphur', 'S', 'gamma'),
        Word('1H', 'hydrogen (1st)', 'H', '1st'),
        Word('2H', 'hydrogen (2nd)', 'H', '2nd'),
        Word('3H', 'hydrogen (3rd)', 'H', '3rd'),
        Word('1HA', 'alpha hydrogen (1st)', 'H', ['alpha', '1st']),
        Word('1HB', 'beta hydrogen (1st)', 'H', ['beta', '1st']),
        Word('1HD', 'delta hydrogen (1st)', 'H', ['delta', '1st']),
        Word('1HD1', 'delta hydrogen (1st of 1st)', 'H', ['delta', '1st', '1st branch']),
        Word('1HD2', 'delta hydrogen (1st of 2nd)', 'H', ['delta', '1st', '2nd branch']),
        Word('1HE', 'epsilon hydrogen (1st)', 'H', ['epsilon', '1st']),
        Word('1HE1', 'epsilon hydrogen (1st of 1st)', 'H', ['epsilon', '1st', '1st branch']),
        Word('1HE2', 'epsilon hydrogen (1st of 2nd)', 'H', ['epsilon', '1st', '2nd branch']),
        Word('1HG', 'gamma hydrogen (1st)', 'H', ['gamma', '1st']),
        Word('1HG1', 'gamma hydrogen (1st of 1st)', 'H', ['gamma', '1st', '1st branch']),
        Word('1HG2', 'gamma hydrogen (1st of 2nd)', 'H', ['gamma', '1st', '2nd branch']),
        Word('1HH1', 'eta hydrogen (1st of 1st)', 'H', ['eta', '1st', '1st branch']),
        Word('1HH2', 'eta hydrogen (1st of 2nd)', 'H', ['eta', '1st', '2nd branch']),
        Word('1HZ', 'zeta hydrogen (1st)', 'H', ['alpha', '1st']),
        Word('2HA', 'alpha hydrogen (2nd)', 'H', ['alpha', '2nd']),
        Word('2HB', 'beta hydrogen (2nd)', 'H', ['beta', '2nd']),
        Word('2HD', 'delta hydrogen (2nd)', 'H', ['delta', '2nd']),
        Word('2HD1', 'delta hydrogen (2nd of 1st)', 'H', ['delta', '2nd', '1st branch']),
        Word('2HD2', 'delta hydrogen (2nd of 2nd)', 'H', ['delta', '2nd', '2nd branch']),
        Word('2HE', 'epsilon hydrogen (2nd)', 'H', ['epsilon', '2nd']),
        Word('2HE1', 'epsilon hydrogen (2nd of 1st)', 'H', ['epsilon', '2nd', '1st branch']),
        Word('2HE2', 'epsilon hydrogen (2nd of 2nd)', 'H', ['epsilon', '2nd', '2nd branch']),
        Word('2HG', 'gamma hydrogen (2nd)', 'H', ['gamma', '2nd']),
        Word('2HG1', 'gamma hydrogen (2nd of 1st)', 'H', ['gamma', '2nd', '1st branch']),
        Word('2HG2', 'gamma hydrogen (2nd of 2nd)', 'H', ['gamma', '2nd', '2nd branch']),
        Word('2HH1', 'eta hydrogen (2nd of 1st)', 'H', ['eta', '2nd', '1st branch']),
        Word('2HH2', 'eta hydrogen (2nd of 2nd)', 'H', ['eta', '2nd', '2nd branch']),
        Word('2HZ', 'zeta hydrogen (2nd)', 'H', ['zeta', '2nd']),
        Word('3HB', 'beta hydrogen (3rd)', 'H', ['beta', '3rd']),
        Word('3HD', 'delta hydrogen (3rd)', 'H', ['delta', '3rd']),
        Word('3HD1', 'delta hydrogen (3rd of 1st)', 'H', ['delta', '3rd', '1st branch']),
        Word('3HD2', 'delta hydrogen (3rd of 2nd)', 'H', ['delta', '3rd', '2nd branch']),
        Word('3HE', 'epsilon hydrogen (3rd)', 'H', ['epsilon', '3rd']),
        Word('3HG1', 'gamma hydrogen (3rd of 1st)', 'H', ['gamma', '3rd', '1st branch']),
        Word('3HG2', 'gamma hydrogen (3rd of 2nd)', 'H', ['gamma', '3rd', '2nd branch']),
        Word('3HZ', 'zeta hydrogen (3rd)', 'H', ['zeta', '3rd']),
        Word('Q', 'pseudo atom'),
        Word('MC', 'pseudo atom', 'Q', 'pseudo')
        ])
libLingua.dictionaryAtom.replicate('1HD1', 'HD11')
libLingua.dictionaryAtom.replicate('2HD1', 'HD12')
libLingua.dictionaryAtom.replicate('3HD1', 'HD13')
libLingua.dictionaryAtom.replicate('1HD2', 'HD21')
libLingua.dictionaryAtom.replicate('2HD2', 'HD22')
libLingua.dictionaryAtom.replicate('3HD2', 'HD23')
#libLingua.dictionaryAtom.replicate('1HE1', 'HE11')
#libLingua.dictionaryAtom.replicate('2HE1', 'HE12')
libLingua.dictionaryAtom.replicate('1HE2', 'HE21')
libLingua.dictionaryAtom.replicate('2HE2', 'HE22')
libLingua.dictionaryAtom.replicate('1HG1', 'HG11')
libLingua.dictionaryAtom.replicate('2HG1', 'HG12')
libLingua.dictionaryAtom.replicate('3HG1', 'HG13')
libLingua.dictionaryAtom.replicate('1HG2', 'HG21')
libLingua.dictionaryAtom.replicate('2HG2', 'HG22')
libLingua.dictionaryAtom.replicate('3HG2', 'HG23')
libLingua.dictionaryAtom.replicate('1HH1', 'HH11')
libLingua.dictionaryAtom.replicate('2HH1', 'HH12')
#libLingua.dictionaryAtom.replicate('3HH1', 'HH13')
libLingua.dictionaryAtom.replicate('1HH2', 'HH21')
libLingua.dictionaryAtom.replicate('2HH2', 'HH22')
#libLingua.dictionaryAtom.replicate('3HH2', 'HH23')
libLingua.dictionaryAtom.replicate('1HA', 'HA1')
libLingua.dictionaryAtom.replicate('2HA', 'HA2')
libLingua.dictionaryAtom.replicate('1HB', 'HB1')
libLingua.dictionaryAtom.replicate('2HB', 'HB2')
libLingua.dictionaryAtom.replicate('3HB', 'HB3')
libLingua.dictionaryAtom.replicate('3HD', 'HD3')
#libLingua.dictionaryAtom.replicate('1HE', 'HE1')
#libLingua.dictionaryAtom.replicate('2HE', 'HE2')
libLingua.dictionaryAtom.replicate('3HE', 'HE3')
libLingua.dictionaryAtom.replicate('1HZ', 'HZ1')
libLingua.dictionaryAtom.replicate('2HZ', 'HZ2')
libLingua.dictionaryAtom.replicate('3HZ', 'HZ3')

libLingua.dictionaryAtom.parse_conjugates()

libLingua.dictionaryMartini = Dictionary('martini',
        [Word('BB', 'backbone'),
        Word('S1', 'side-chain (1st)'), 
        Word('S2', 'side-chain (2nd)'), 
        Word('S3', 'side-chain (3rd)'), 
        Word('S4', 'side-chain (4th)'), 
        Word('BBc', 'backbone (coil)', 'BB', 'coil'),
        Word('BBh', 'backbone (helix)', 'BB', 'helix'),
        Word('BBs', 'backbone (bend)', 'BB', 'bend'),
        Word('BBt', 'backbone (turn)', 'BB', 'turn'),
        Word('BBb', 'backbone (beta)', 'BB', 'beta'),
        Word('BBg', 'backbone (3-helix)', 'BB', '3-helix'),
        Word('BBi', 'backbone (5-helix)', 'BB', '5-helix'),
        Word('BBe', 'backbone (extended)', 'BB', 'extended'),
        Word('S1c', 'side-chain (1st, coil)', 'S1', 'coil'),
        Word('S1h', 'side-chain (1st, helix)', 'S1', 'helix'),
        Word('S1s', 'side-chain (1st, bend)', 'S1', 'bend'),
        Word('S1t', 'side-chain (1st, turn)', 'S1', 'turn'),
        Word('S1b', 'side-chain (1st, beta)', 'S1', 'beta'),
        Word('S1g', 'side-chain (1st, 3-helix)', 'S1', '3-helix'),
        Word('S1i', 'side-chain (1st, 5-helix)', 'S1', '5-helix'),
        Word('S1e', 'side-chain (1st, extended)', 'S1', 'extended'),
        Word('S2c', 'side-chain (2nd, coil)', 'S2', 'coil'),
        Word('S2h', 'side-chain (2nd, helix)', 'S2', 'helix'),
        Word('S2s', 'side-chain (2nd, bend)', 'S2', 'bend'),
        Word('S2t', 'side-chain (2nd, turn)', 'S2', 'turn'),
        Word('S2b', 'side-chain (2nd, beta)', 'S2', 'beta'),
        Word('S2g', 'side-chain (2nd, 3-helix)', 'S2', '3-helix'),
        Word('S2i', 'side-chain (2nd, 5-helix)', 'S2', '5-helix'),
        Word('S2e', 'side-chain (2nd, extended)', 'S2', 'extended'),
        Word('S3c', 'side-chain (3rd, coil)', 'S3', 'coil'),
        Word('S3h', 'side-chain (3rd, helix)', 'S3', 'helix'),
        Word('S3s', 'side-chain (3rd, bend)', 'S3', 'bend'),
        Word('S3t', 'side-chain (3rd, turn)', 'S3', 'turn'),
        Word('S3b', 'side-chain (3rd, beta)', 'S3', 'beta'),
        Word('S3g', 'side-chain (3rd, 3-helix)', 'S3', '3-helix'),
        Word('S3i', 'side-chain (3rd, 5-helix)', 'S3', '5-helix'),
        Word('S3e', 'side-chain (3rd, extended)', 'S3', 'extended'),
        Word('S4c', 'side-chain (4th, coil)', 'S4', 'coil'),
        Word('S4h', 'side-chain (4th, helix)', 'S4', 'helix'),
        Word('S4s', 'side-chain (4th, bend)', 'S4', 'bend'),
        Word('S4t', 'side-chain (4th, turn)', 'S4', 'turn'),
        Word('S4b', 'side-chain (4th, beta)', 'S4', 'beta'),
        Word('S4g', 'side-chain (4th, 3-helix)', 'S4', '3-helix'),
        Word('S4i', 'side-chain (4th, 5-helix)', 'S4', '5-helix'),
        Word('S4e', 'side-chain (4th, extended)', 'S4', 'extended')
        ])
libLingua.dictionaryMartini.parse_conjugates()


libLingua.allDictionaries = [ \
        libLingua.dictionaryAmino1, \
        libLingua.dictionaryAmino3, \
        libLingua.dictionaryNMR, \
        libLingua.dictionaryEnglish, \
        libLingua.dictionaryMara, \
        libLingua.dictionaryAtom, \
        libLingua.dictionaryMartini]

#############################################################################
# CODEX DATABASE                                                            #
#############################################################################

"""
Conversions btw the standard Mara notation and those in general use among the
NMR community.
"""
libLingua.codexNMR = Codex(('mara', 'general NMR'))
libLingua.codexNMR.associate('DHaCa', 'HACA')
libLingua.codexNMR.associate('DCaCb', 'CACB')
libLingua.codexNMR.associate('DHaCo', 'HACO')
libLingua.codexNMR.associate('DCaCo', 'CACO')
libLingua.codexNMR.associate('DHaN', 'HAN')
libLingua.codexNMR.associate('DNCa', 'NCA')
libLingua.codexNMR.associate('DHnN', 'HN')
libLingua.codexNMR.associate('DCoN', 'CON')
libLingua.codexNMR.associate('DHnCo', 'HCO')

"""
Conversions btw one and three letter aminoacid symbols.
"""
libLingua.codexAminoSymbols = Codex(('1-letter-aminoacids','3-letter-aminoacids'))
libLingua.codexAminoSymbols.associate('G','Gly')
libLingua.codexAminoSymbols.associate('A','Ala')
libLingua.codexAminoSymbols.associate('V','Val')
libLingua.codexAminoSymbols.associate('L','Leu')
libLingua.codexAminoSymbols.associate('I','Ile')
libLingua.codexAminoSymbols.associate('P','Pro')
libLingua.codexAminoSymbols.associate('C','Cys')
libLingua.codexAminoSymbols.associate('M','Met')
libLingua.codexAminoSymbols.associate('H','His')
libLingua.codexAminoSymbols.associate('F','Phe')
libLingua.codexAminoSymbols.associate('Y','Tyr')
libLingua.codexAminoSymbols.associate('W','Trp')
libLingua.codexAminoSymbols.associate('N','Asn')
libLingua.codexAminoSymbols.associate('Q','Gln')
libLingua.codexAminoSymbols.associate('S','Ser')
libLingua.codexAminoSymbols.associate('T','Thr')
libLingua.codexAminoSymbols.associate('K','Lys')
libLingua.codexAminoSymbols.associate('R','Arg')
libLingua.codexAminoSymbols.associate('D','Asp')
libLingua.codexAminoSymbols.associate('E','Glu')

"""
Conversions btw english and three letter aminoacid symbols.
"""
libLingua.codexAmino = Codex(('english','3-letter-aminoacids'))
libLingua.codexAmino.associate('alanine','Ala')
libLingua.codexAmino.associate('arginine','Arg')
libLingua.codexAmino.associate('asparagine','Asn')
libLingua.codexAmino.associate('aspartate','Asp')
libLingua.codexAmino.associate('cysteine','Cys')
libLingua.codexAmino.associate('glutamine','Gln')
libLingua.codexAmino.associate('glutamate','Glu')
libLingua.codexAmino.associate('glycine','Gly')
libLingua.codexAmino.associate('histidine','His')
libLingua.codexAmino.associate('isoleucine','Ile')
libLingua.codexAmino.associate('leucine','Leu')
libLingua.codexAmino.associate('lysine','Lys')
libLingua.codexAmino.associate('methionine','Met')
libLingua.codexAmino.associate('phenylalanine','Phe')
libLingua.codexAmino.associate('proline','Pro')
libLingua.codexAmino.associate('serine','Ser')
libLingua.codexAmino.associate('threonine','Thr')
libLingua.codexAmino.associate('tryptophan','Trp')
libLingua.codexAmino.associate('tyrosine','Tyr')
libLingua.codexAmino.associate('valine','Val')
libLingua.codexAmino.associate('selenomethionine','Mse')

libLingua.allCodices = [libLingua.codexNMR, libLingua.codexAminoSymbols, \
        libLingua.codexAmino]

#omniTranslator = Translator(libLingua.allDictionaries, libLingua.allCodices)
#NMRTranslator = Translator(\
#        [libLingua.dictionaryMara, libLingua.dictionaryNMR], \
#        [libLingua.codexNMR], 'mara')

