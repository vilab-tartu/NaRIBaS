#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
#---------------------------------------------------------------------------#
# Function: Calculate an average dataset over all datasets in one or more   #
#           xvg-files.                                                      #
# Usage: average-xvg.py [options] <file1.xvg> {<file2.xvg> ... }            #
# Help: run join-xvg.py --help                                              #
# Author: Martti Louhivuori (m.j.louhivuori@rug.nl)                         #
# Version: 1.0 (20.06.2008)                                                 #
#---------------------------------------------------------------------------#
from optparse import OptionParser
import logging, sys, re
from numpy import array, sqrt
from Mara.IO import XVGIO

if __name__ == '__main__':
    usage = 'usage: %prog [options] <file1.xvg> <file2.xvg> {<file3.xvg> ... }'
    desc = 'Calculate an average dataset over all datasets in one or ' + \
            'more xvg-files.'
    parser = OptionParser(usage=usage, description=desc)
    parser.add_option('-o', '--output', metavar='FILE', default=None,
            help='write output to FILE (default: STDOUT)')
    parser.add_option('-l', '--label', action='store_true', default=False, 
            help='create legend keys from filenames')
    parser.add_option('--verbose', action='store_true', default=False,
            help='display additional information while running')
    parser.add_option('--debug', action='store_true', default=False, 
            help='run in debug mode, i.e. maximum information')

    options, args = parser.parse_args()

    # set logger format etc.
    logging.basicConfig(level=logging.WARNING, format='%(levelname)s ' + \
            '%(message)s @ %(asctime)s %(module)s line %(lineno)s',
            datefmt='%H:%M:%S')
    # set logging thresholds
    if options.debug:
        logging.getLogger('').setLevel(logging.DEBUG)
    elif options.verbose:
        logging.getLogger('').setLevel(logging.WARNING)
    else:
        logging.getLogger('').setLevel(logging.CRITICAL)
    logging.debug('options: %s' % repr(options))
    logging.debug('args: %s' % repr(args))

    # redirect STDOUT?
    if options.output:
        try:
            sys.stdout = open(options.output, 'w')
        except IOError, errmsg:
            print '#', errmsg
            print '# Directing output to STDOUT instead.'

    # read xvg files
    legends = []
    commands = None
    keys = None
    sets = []
    for file in args:
        xvg = XVGIO(file)
        key, set = xvg.read()
        logging.debug('file=' + repr(file))
        logging.debug('key=' + repr(key))
        logging.debug('set=' + repr(set))
        if keys is None:
            keys = key
        sets.extend(set)
        if commands is None and len(xvg.commands):
            commands = xvg.commands
        if options.label:
            legends.append(file.strip('.xvg'))
        elif len(xvg.legends):
            legends.extend(xvg.legends)
    logging.debug('keys=' + repr(keys))
    logging.debug('sets=' + repr(sets))
    logging.debug('commands=' + repr(commands))
    logging.debug('legends=' + repr(legends))
    # calculate average
    matrix = array(sets)
    avg = sum(matrix) / len(matrix)
    dev = sqrt( sum((matrix - avg)**2) / float(len(matrix) * (len(matrix) -1)) )
    # prepare output
    xvg.comments = ['# Created by average-xvg.py from', 
            '# ' + '  %s' % ' '.join(args)]
    xvg.legends = ['average', 'std.err.']
    xvg.sets = [avg, dev]
    # output average dataset
    xvg.echo_meta()
    xvg.echo_sets()

    # the end.
    if options.output:
        sys.stdout.close()

