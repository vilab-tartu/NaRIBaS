#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.Interface (v. 0.1):
    Front-ends to various external programs.

    Cyana -- automated structure calculations w/ NMR constraints
    Pales -- prediction of alignment from structure
    
    Requirements: Python 2.4->
                  Cyana
                  Pales

    TODO: - should Pales input parameters be checked?
          - more features to Cyana
    
    Author: Martti Louhivuori (martti.louhivuori@helsinki.fi)

    Date: 14.2.2006

    ---

    Copyright (C) 2006  Martti Louhivuori

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    The full license text can be found in the file LICENSE.
"""
import os, sys, re
from random import Random
from Mara.aux import ArgumentError
from Mara import loki

class Cyana:
    """
    Front-end for CYANA. Provides methods needed to run CYANA.

    Current implementation includes only methods needed to run CYANA in batch 
    mode using macros. TODO: more features, control of screen output, support
    for a session mode CYANA, i.e. a true API to CYANA/INCLAN.
    """
    def __init__(self, input, output=None, fixed_angles=None, rg_limits=None, 
            use_limits=True, seed=None, state=None, jump=None, tally=0, 
            tally_digits=3, name=None, verbose=False, debug=False):
        # store and/or process inputs
        self.input = str(input)
        self.rg_limits = rg_limits
        self.use_limits = bool(use_limits)
        self.tally = int(tally)
        self.tally_digits = int(tally_digits)
        if output != None:
            self.output = str(output)
        else:
            self.output = self.input + '-'
        self.fixed_angles = fixed_angles
        if fixed_angles != None:
            if type(fixed_angles) != list:
                if type(fixed_angles) == tuple:
                    self.fixed_angles = list(fixed_angles)
                else:
                    self.fixed_angles = [fixed_angles]
            for i in range(len(self.fixed_angles)):
                if type(self.fixed_angles[i]) in [list, tuple]:
                    lo, hi = str(self.fixed_angles[i][0]), \
                            str(self.fixed_angles[i][1])
                elif re.search('([0-9]+).+([0-9]+)', str(self.fixed_angles[i])):
                    lo, hi = re.search('([0-9]+).+([0-9]+)', \
                            str(self.fixed_angles[i])).groups()
                else:
                    lo = hi = str(self.fixed_angles[i])
                self.fixed_angles[i] = '%s..%s'%(lo,hi)
        if name is None:
            self.name = 'temp-%d'%os.getpid()
        else:
            self.name = str(name)
        if debug: loki.setLevel(logging.DEBUG)
        elif verbose: loki.setLevel(logging.INFO)
        # check that required input files exist
        if not os.path.isfile(self.input+'.seq'):
            raise ArgumentError, ('No sequence file found.', 'Cyana')
        if self.fixed_angles and not os.path.isfile(self.input + '.pdb'):
            raise ArgumentError, ('No model PDB file found for angle fixes.', 
                    'Cyana')
        # setup temporary files
        os.system('cp -f %s.seq %s.seq'%(self.input, self.name))
        loki.info('Readied the sequence file for further use.')
        if self.use_limits:
            for suffix in ['lol', 'upl', 'aco', 'xplor']:
                if os.path.isfile(self.input + '.' + suffix):
                    os.system('cp -f %s.%s %s.%s'%(self.input, suffix, 
                        self.name, suffix))
                    loki.info("Readied a '%s' file for further use."%suffix)
        if self.fixed_angles:
            os.system('cp -f %s.pdb %s.pdb'%(self.input, self.name))
        # initialise the random generator
        self.randgen = Random(seed)
        if state == None:
            if jump != None:
                if sys.version < '2.2':
                    loki.warn('Python older than 2.2. Problems may arise.')
                elif sys.version > '2.2' and sys.version < '2.3':
                    # this should provide enough room for each run to be in 
                    # unique state space
                    self.randgen.jumpahead(jump*self.amount*1.1)
                else:
                    self.randgen.jumpahead(jump)
        else:
            self.randgen.setstate(state)
        self.__last_macro_size__ = None

    def anneal(self, amount, snippets=None, max_temp=4.0, min_temp=0.9, 
        steps=5000):
        """
        Anneal.
        """
        amount = int(amount)
        if snippets is None:
            snippets = amount
        else:
            snippets = int(snippets)
        max_temp = float(max_temp)
        min_temp = float(min_temp)
        steps = int(steps)
        if snippets > amount:
            loki.warn('Snippet size larger than total amount. '
                    'Running in one go.')
            snippets = amount
# Current version of Cyana (2.1?) seems to be broken as tend is not
# supported by anneal. FIXME when corrected.
#        self.__calc_end__ = 'anneal thigh=%1.3f tend=%1.3f steps=%d'%( \
#                max_temp, min_temp, steps)
        self.__calc_end__ = 'anneal thigh=%1.3f steps=%d'%( \
                max_temp, steps)
        if len(str(amount)) > self.tally_digits:
            self.tally_digits = len(str(amount))
        loki.info('--- Anneal ---\nthigh=%1.3f tend=%1.3f steps=%d\n'%( \
                max_temp, min_temp, steps) + '-'*14)
        done = 0
        while done < amount:
            if amount-done > snippets:
                batch = snippets
            else:
                batch = amount - done
            self.run(batch)
            done += batch
            loki.info('%d structures ready'%done)
        loki.info('all done. Thank you for all the fish.')

    def checkRg(self, filename):
        """
        Check whether a given conformation is within the limits set for radius
        of gyration.
        """
        rg = calculate_rg(filename)
        loki.info('Rg(%s) = %f)'%(filename, rg))
        if (self.rg_limits[0] != None and rg < self.rg_limits[0]) or \
                (self.rg_limits[1] != None and rg > self.rg_limits[1]):
            return False
        return True

    def run(self, size):
        """
        Run CYANA using a macro created earlier.
        """
        loki.info('A new Cyana round started for %d structures.'%size)
        # if no Rg check, point to something always true 
        if self.rg_limits:
            check_rg = self.checkRg
        else:
            check_rg = len
        # run cyana as many time as needed until it completes succesfully
        run_cyana = "cyana < %s.cmd" % self.name
        todo = size
        while todo:
            # check if a new macro is needed
            if todo != self.__last_macro_size__:
                self.writeMacro(todo)
            self.writeScript()
            os.system(run_cyana)
            for i in range(1,todo+1):
                name = '%s-%03d.pdb'%(self.name, i)
                if os.path.exists(name) and check_rg(name):
                    todo -= 1
                    self.tally += 1
                    os.system('mv -f %s %s%s.pdb'%(name, self.output, 
                        str(self.tally).zfill(self.tally_digits)))
                else:
                    os.system('rm -f %s'%name)

    def writeMacro(self, size):
        """
        Write all the required CYANA-commands into a macro.
        """
        loki.info('A new Cyana macro required for %d structures.'%size)
        outf = open(self.name+'.cya', 'w')
        outf.write('cyanalib\n')
        outf.write('readdata %s\n'%self.name)
        # all angles in given amino acids as well as all omegas should be fixed
        if self.fixed_angles:
            outf.write('angle select *\n')
            outf.write('angle fix\n')
            outf.write('angle select -OMEGA\n')
            for fix in self.fixed_angles:
                outf.write('angle select -%s\n' % fix)
            outf.write('angle free\n')
            outf.write('read pdb %s\n' % self.name)
        outf.write('calc_all %d %s\n'% (size, self.__calc_end__))
        outf.write('write_all %s- pdb\n'%self.name)
# this could be used to find out the end seed
#        outf.write('system echo $seed > %s.seed'%self.name)
        outf.close()
        self.__last_macro_size__ = size

    def writeScript(self):
        """
        Write a top-level macro used to launch CYANA.
        """
        outf = open(self.name+'.cmd', 'w')
        seed = int(self.randgen.uniform(0,2147483648))
        outf.write('seed = %d\n' % seed)
        outf.write('%s\n' % self.name)
        outf.write('quit\n')
        outf.close()

    def shutdown(self):
        """
        Perform a clean exit by deleting all temporary files.
        """
        loki.info('Removing temporary files.')
        cmd = "rm -f %s*" % self.name
        os.system(cmd)


class Pales:
    def __init__(self, command='PALES', verbose=False, debug=False, **free):
#euler, noeuler, qDa, qRms, qStd, w, excl, incl, verb, noverb, yzInv, noyzInv, rem, fixedDI, nofixedDI): pass
        self.args, self.options = self.__parse_free__(free)
        self.command = str(command)
        self.verbose = bool(verbose)
        self.debug = bool(debug)
        if debug:  loki.setLevel(logging.DEBUG)
        elif verbose: loki.setLevel(logging.INFO)
    
    def __parse_free__(self, free, args={}, options=[]):
        for key, value in free.items():
            if value in [None, '']:
                options.append(key)
            else:
                args[key] = value
        return args, options

    def __join__(self, args, options):
        txt = ''
        for key, value in args.items():
            txt += ' -%s %s'%(str(key), str(value))
        if len(options):
            txt += ' -' + ' -'.join(options)
        return txt

    def run(self, txt=''):
        loki.info('arguments: %s'%txt)
        os.system('%s %s'%(self.command, txt))

    def steric(self, args={}, options=[], **free):
#stPales, pdb, inD, outD, pdbRot, rotID, bic, pf1, wv, dot, digPsi, dGrid, rM, lcS, surf, nosurf, rA, H, noH, res, r1, rN, atom, resName, chain, seg, id, inv, clear, set, exact, regex, glob, aExact, aRegex, aGlob, rExact, rRegex, rGlob, cExact, cRegex, cGlob, sExact, sRegex, sGlob, a1, s1, sN, aD1, aD2): pass
        txt = self.__join__(*self.__parse_free__(free, \
                dict(self.args.items() + args.items()), self.options + options))
        self.run('-stPales' + txt)
    def electrostatic(self, args={}, options=[], **free):
#elPales, pdb, inD, pot, pka, outD, outPka, pdbRot, rotID, bic, pf1, wv, dot, digPsi, dGrid, rM, lcS, surf, nosurf, rA, pH, nacl, chSurf, cut, pkaDef, nopkaDef, nocharmm, charmm, chTerm, nochTerm, kEne, kPot, kMono, kDipo, kQuad, temp, dIFace, noelInfo, elInfo, nodipole, dipole, nochSel, chSel, H, noH, res, r1, rN, atom, resName, chain, seg, id, inv, clear, set, exact, regex, glob, aExact, aRegex, aGlob, rExact, rRegex, rGlob, cExact, cRegex, cGlob, sExact, sRegex, sGlob, a1, s1, sN, aD1, aD2):pass
        txt = self.__join__(*self.__parse_free__(free, \
                dict(self.args.items() + args.items()), self.options + options))
        self.run('-elPales' + txt)
    def bestfit(self, args={}, options=[], **free):
#bestFit, single, jack, mcDc, mcStruc, pdb, inD, outD, pdbRot, rotID, map, outDa, outMap, outAng, outFail, err, noerr, vecNoise, angNoise, stdCone, nofixCone, fixCone, digCone, autoCone, noautoCone, lCone, incCone, nCone, kAdjustMc, autoMc, noautoMc, nAdjustMc, fracAdjustMc, kLAdjustMc, incAdjustMc, dFrac, nDfixed, nonDfixed, svd, saupe, fixed, nofixed, dadr, da, daMin, daMax, dr, drMin, drMax, psi, psiMax, psiMin, theta, thetaMax, thetaMin, phi, phiMax, phiMin, res, r1, rN, atom, resName, chain, seg, id, inv, clear, set, exact, regex, glob, aExact, aRegex, aGlob, rExact, rRegex, rGlob, cExact, cRegex, cGlob, sExact, sRegex, sGlob, a1, s1, sN, aD1, aD2):pass
        txt = self.__join__(*self.__parse_free__(free, \
                dict(self.args.items() + args.items()), self.options + options))
        self.run('-bestFit' + txt)
    def bestfit_unassigned(self, args={}, options=[], **free):
#bestFit, pdb, inD, exhaust, fixAa, fixRes, dcID, exhaustSel, nSel, dcOpt, nodcOpt, rand, map, outD, pdbRot, rotID, lDa, hDa, lR, hR, res, r1, rN, atom, resName, chain, seg, id, inv, clear, set, exact, regex, glob, aExact, aRegex, aGlob, rExact, rRegex, rGlob, cExact, cRegex, cGlob, sExact, sRegex, sGlob, a1, s1, sN, aD1, aD2):pass
        txt = self.__join__(*self.__parse_free__(free, \
                dict(self.args.items() + args.items()), self.options + options))
        self.run('-bestFit' + txt)
    def likelihood(self, args={}, options=[], **free):
#daMl, inD, outD, outSc, lDa, hDa, incDa, lR, hR, incR, s1, sN, aD1, aD2, nomapMl, mapMl):pass
        txt = self.__join__(*self.__parse_free__(free, \
                dict(self.args.items() + args.items()), self.options + options))
        self.run('-daMl' + txt)
    def freeformat(self, args={}, options=[], **free):
#stPalesFree, pdbF, outA, bic, pf1, wv, dot, digPsi, dGrid, rM, lcS, surf, nosurf, rA):pass
        txt = self.__join__(*self.__parse_free__(free, \
                dict(self.args.items() + args.items()), self.options + options))
        self.run('-stPalesFree' + txt)
    def histogram(self, args={}, options=[], **free):
#daHist, inD, outD, bin, nFit, sdevDzz, sdevDyy, daLo, daHi, rhLo, rhHi, s1, sN, aD1, aD2, noave, ave):pass
        txt = self.__join__(*self.__parse_free__(free, \
                dict(self.args.items() + args.items()), self.options + options))
        self.run('-daHist' + txt)
    def conversion(self, args={}, options=[], **free):
#conv, inD, outD, xplor, cns, dyana, tbl, pales, rosetta, tensID, scale, unscale, dcNoise, s1, sN, aD1, aD2, nodcHH, dcHH, euler, noeuler, qDa, qRms, qStd, w, excl, incl, verb, noverb, yzInv, noyzInv, rem, fixedDI, nofixedDI):pass
        txt = self.__join__(*self.__parse_free__(free, \
                dict(self.args.items() + args.items()), self.options + options))
        self.run('-conv' + txt)
    def noise(self, args={}, options=[], **free):
#struc, pdb, inD, pdbRot, idum, stdCone, nofixCone, fixCone, digCone, surf, nosurf, H, noH, res, r1, rN, atom, resName, chain, seg, id, inv, clear, set, exact, regex, glob, aExact, aRegex, aGlob, rExact, rRegex, rGlob, cExact, cRegex, cGlob, sExact, sRegex, sGlob, a1, s1, sN, aD1, aD2):pass
        txt = self.__join__(*self.__parse_free__(free, \
                dict(self.args.items() + args.items()), self.options + options))
        self.run('-struc' + txt)
    def analyseDC(self, args={}, options=[], **free):
#anDc, inD1, inD2, outD, s1, sN, aD1, aD2, add, sub, mul, div, c1, c2): pass
        txt = self.__join__(*self.__parse_free__(free, \
                dict(self.args.items() + args.items()), self.options + options))
        self.run('-anDc' + txt)
    def analysePDB(self, args={}, options=[], **free):
#anPdb, pdb, inS, pka, outP, outPka, pdbRot, rotID, H, noH, res, r1, rN, atom, resName, chain, seg, id, inv, clear, set, exact, regex, glob, aExact, aRegex, aGlob, rExact, rRegex, rGlob, cExact, cRegex, cGlob, sExact, sRegex, sGlob, pH, pkaDef, nopkaDef, nocharmm, charmm, chTerm, nochTerm, kMono, kDipo, kQuad, nodipole, dipole, nochSel, chSel, noelInfo, elInfo, user, ster, el):pass
        txt = self.__join__(*self.__parse_free__(free, \
                dict(self.args.items() + args.items()), self.options + options))
        self.run('-anPdb' + txt)
    def analyseA(self, args={}, options=[], **free):
#anA, inS1, inS2, outA, add, sub, mul, div, max, c1, c2, twist, rotID, hinge): pass
        txt = self.__join__(*self.__parse_free__(free, \
                dict(self.args.items() + args.items()), self.options + options))
        self.run('-anA' + txt)
    def multiple(self, function, batch, args={}, options=[], **free):
        keys = batch.keys()
        n = []
        for key in keys: n.append(len(batch[key]))
        if len(n) != n.count(n[0]):
            raise ArgumentError, 'Batch items need to be of equal length.'
        else: n = n[0]
        for i in range(n):
            for key in keys: args[key] = batch[key][i]
            function(args, options, **free)


