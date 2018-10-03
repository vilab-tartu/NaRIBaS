#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.mpi (v. 0.1):
    A multiprocessor interface. (Untested)

    Factory -- A batch runner that uses multiple processing nodes.
    Worker  -- Interface to a single processing node.
    
    Requirements: Python 2.4->
    
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
from os import spawnvp, spawnlp, P_NOWAIT, waitpid
from time import time

class Factory:
    """
    A batch runner that uses multiple processing nodes.
    """
    def __init__(self, nodes, jobs=[], workers=None, options=[]):
        if workers == None:
            workers = len(nodes)
        self.idle = []
        for i in range(workers):
            self.idle.append(Worker(nodes[i % len(nodes)], options))
        self.busy = {}
        self.queue = jobs

    def produce(self, amount='all'):
        times = []
        while (amount and self.queue) or self.busy:
            if self.idle and amount and self.queue:
                w = self.idle.pop(0)
                pid = w.engage(*self.queue.pop(0))
                self.busy[pid] = w
                if amount != 'all':
                    amount -= 1
            else:
                pid, status = waitpid(-1,0)
                if pid in self.busy.keys():
                    w = self.busy[pid]
                    del self.busy[pid]
                    t = w.release()
                    times.append(t)
                    self.idle.append(w)
        sum = 0
        for i in range(len(times)):
            sum += times[i]
        return (sum, sum/len(times))

class Worker:
    """
    Interface to a single processing node.
    """

    def __init__(self, node, options=[], command='mosrun'):
        self.node = node
        self.options = options
        self.command = command
        self.start = 0
        self.pid = 0
        self.busy = False

    def engage(self, job, options=[]):
        ops = [self.command, '-j%d' % self.node, job]
        ops.extend(self.options)
        ops.extend(options)
        pid = spawnvp(P_NOWAIT, self.command, ops)
        self.start = time()
        self.pid = pid
        return pid

    def release(self):
        t = time()-self.start
        self.pid = 0
        self.busy = False
        return t

