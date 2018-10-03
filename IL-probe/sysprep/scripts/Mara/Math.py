#!/usr/bin/python
# -*- coding: ISO-8859-1 -*-
"""
Mara.Math (v. 1.3):
    Common math routines.

    PearsonCorrelation          -- Calculate the Pearson correlation between 
                                   two vectors.
    CornilescuQ                 -- Calculate Cornilescu q-factor for observed 
                                   RDCs.
    CosinesArbitraryAxes        -- Calculate the cosines from a vector to 
                                   arbitrary axes.
    EulerRotationMatrix         -- Calculate rotation matrix using Euler angles.
    RotateAroundArbitraryLine   -- Rotate a point around an arbitrary line 
                                   segment.
    RotationMatrixArbitraryLine -- A rotation matrix for a rotation around the
                                   direction of an arbitrary line segment.
    InertiaMatrix               -- Matrix of inertia for the calculation of 
                                   an inertia tensor.
    makeOrderMatrix             -- Make an order matrix from a list of cosines.
    gdo                         -- Calculate the general degree of order from 
                                   diagonal Saupe elements.
    gdo5e                       -- Calculate the general degree of order from 
                                   Saupe elements.
    sampleAverage               -- Calculate the average value w/ error from a 
                                   list of values.
    ParetoOptimal               -- Determine a pareto optimal co-ordinate set.

    Requirements: Python 2.2->
                  numpy
    
    Author: Martti Louhivuori (m.j.louhivuori@rug.nl), 6.10.2005

    Date: 07.10.2010

    ---

    Copyright (C) 2006-2008  Martti Louhivuori

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    The full license text can be found in the file LICENSE.
"""

from numpy import sum, sqrt, array, dot, sin, cos, pi, zeros, float, \
        transpose, min, max

def PearsonCorrelation(x, y):
    """
    Calculate the Pearson correlation between two vectors.
    """
    if len(x) != len(y):
        raise TypeError, "Vectors need to be of same length."
    x, y = array(x), array(y)
    n = float(len(x))
    sx, sy, sx2, sy2 = sum(x), sum(y), sum(x**2), sum(y**2)
    return (sum(x*y)-sx*sy/n) / sqrt((sx2-sx**2/n)*(sy2-sy**2/n))

def CornilescuQ(reference, observed):
    """
    Calculate Cornilescu q-factor for observed RDCs.
    """
    diff = 0.0
    total = 0.0
    for x,y in zip(reference, observed):
        diff += (x - y)**2
        total += x**2
    return sqrt(diff/total)

def CosinesArbitraryAxes(vectors, axes=[[1,0,0],[0,1,0],[0,0,1]]):
    """
    Calculate the cosines from a vector to arbitrary axes.
    """
    axes = array(axes, float)
    vectors = array(vectors, float)
    cosines = transpose( transpose(dot(vectors, axes)) / \
            sqrt(sum(vectors**2, 1)) )
    return [tuple(x) for x in cosines]

def EulerRotationMatrix(alpha, beta, gamma, order='zyz', direction='ccw'):
    """
    Calculate rotation matrix using Euler angles.
    
    Assumes counter-clockwise rotations around z-y'-z''. Rotation direction 
    can be given as 'cw' for a clockwise rotation. Currently implemented 
    rotation orders: 'zyz'.
    
    TODO: add more rotation orders / generalise the generation of the matrix 
    """
    alpha, beta, gamma = float(alpha)/180.0*pi, float(beta)/180.0*pi, \
            float(gamma)/180.0*pi
    sa, sb, sg = sin((alpha, beta, gamma))
    ca, cb, cg = cos((alpha, beta, gamma))
    rotmat = zeros((3,3),float)
    if order == 'zyz':
        rotmat[0,0] = cg * cb * ca - sg * sa
        rotmat[0,1] = -cg * cb * sa - sg * ca
        rotmat[0,2] = cg * sb
        rotmat[1,0] = sg * cb * ca + cg * sa
        rotmat[1,1] = -sg * cb * sa + cg * ca
        rotmat[1,2] = sg * sb
        rotmat[2,0] = -sb * ca
        rotmat[2,1] = sb * sa
        rotmat[2,2] = cb
        asymmetric = [(0,1),(0,2),(1,0),(2,0)]
    if direction == 'cw':
        for a in asymmetric:
            rotmat[a] *= -1
    return rotmat

def RotateAroundArbitraryLine(p, theta, p1, p2):
    """
    Rotate a point p anticlockwise by angle theta around an arbitrary line 
    segment p1->p2 and return the rotated point assuming a right handed 
    coordinate system.
    """
    p = array(p, float)
    p1 = array(p1, float)

    matrix = RotationMatrixArbitraryLine(theta, p1, p2)

    return dot(matrix, p - p1) + p1

def RotationMatrixArbitraryLine(theta, p1, p2):
    """
    Return a rotation matrix for an anticlockwise rotation around the 
    direction of an arbitrary line segment p1->p2.

    Note: to rotate around an arbitrary line, one must also translate the 
    coordinates prior to a rotation (see RotateAroundArbitraryLine)
    """
    q = zeros((3,3),float)
    p1 = array(p1, float)
    p2 = array(p2, float)

    r = p2 - p1
    r /= sqrt(sum(r**2))
    costheta = cos(theta)
    sintheta = sin(theta)

    q[0,0] = costheta + (1 - costheta) * r[0]**2
    q[0,1] = (1 - costheta) * r[0] * r[1] - r[2] * sintheta
    q[0,2] = (1 - costheta) * r[0] * r[2] + r[1] * sintheta
    
    q[1,0] = (1 - costheta) * r[0] * r[1] + r[2] * sintheta
    q[1,1] = costheta + (1 - costheta) * r[1]**2
    q[1,2] = (1 - costheta) * r[1] * r[2] - r[0] * sintheta

    q[2,0] = (1 - costheta) * r[0] * r[2] - r[1] * sintheta
    q[2,1] = (1 - costheta) * r[1] * r[2] + r[0] * sintheta
    q[2,2] = costheta + (1 - costheta) * r[2]**2

    return q

def InertiaMatrix(x, y, z):
    """
    Calculate the inertia matrix for a single mass sans mass.
    """
    q = zeros((3,3), float)
    q[0,0] = y**2 + z**2
    q[0,1] = -x * y
    q[0,2] = -x * z
    q[1,0] = -y * x
    q[1,1] = x**2 + z**2
    q[1,2] = -y * z
    q[2,0] = -z * x
    q[2,1] = -z * y
    q[2,2] = x**2 + y**2
    return q 

def makeCosineMatrix(cosines):
    """
    Make a cosine matrix a la SVD-RDCs from a list of cosines.
    """
    matrix = zeros((len(cosines),5), float)
    for i in range(len(cosines)):
        x, y, z = cosines[i]
        matrix[i] = [y**2 - x**2, z**2 - x**2, 2 * x * y, 2 * x * z, 2 * y * z]
    return matrix
makeOrderMatrix = makeCosineMatrix # FIXME: remove?

def gdo(Sxx, Syy, Szz):
    """
    Calculate the general degree of order from diagonal Saupe elements.
    """
    return sqrt( 2/3. * (Sxx**2 + Syy**2 + Szz**2) )

def gdo5e(Szz, Syy, Sxy, Sxz, Syz):
    """
    Calculate the general degree of order from Saupe elements.
    """
    return sqrt( 2/3. * ((Szz + Syy)**2 + Syy**2 + Szz**2 - 2*Sxy*Sxz \
            - 2*Sxy*Syz - 2*Sxz*Syz) )

def sampleAverage(x):
    """
    Calculate the average value w/ error from a list of values.
    """
    avg = sum(x) / float(len(x))
    error = 0.0
    for e in x:
        error += (avg - e)**2
    error /= float(len(x) * (len(x) - 1))
    return avg, sqrt(error)

class ParetoOptimal:
    """
    Determine a pareto optimal co-ordinate set.
    """
    def __init__(self):
        self.__quads__ = [None] * 8
        self.__ids__ = [None] * 8
    def get_quadrant(self, point):
        id = 0
        if point[0] < 0.0:
            id += 4
        if point[1] < 0.0:
            id += 2
        if point[2] < 0.0:
            id += 1
        return id
    def accept(self, point):
        quad = self.get_quadrant(point)
        if self.__quads__[quad] is None:
            return quad
        elif self.__quads__[quad].ndim == 1:
            if min(abs(self.__quads__[quad]) - abs(array(point))) < 0.0:
                return quad
        else:
            los = min(abs(self.__quads__[quad]) - abs(array(point)), 1)
            if max(los) < 0.0:
                return quad
        return -1
    def add(self, point, id=None):
        quad = self.accept(point)
        if quad >= 0:
            loki.debug('accepted in quadrant %s' % repr(quad))
            if self.__quads__[quad] is None:
                self.__quads__[quad] = array(point)
                self.__ids__[quad] = [id]
            else:
                # append new point
                flat = self.__quads__[quad].tolist()
                if self.__quads__[quad].ndim == 1:
                    flat = [flat]
                flat.append(point)
                # remove obsolete points
                if self.__quads__[quad].ndim == 1:
                    diff = [max(abs(self.__quads__[quad]) - abs(array(point)))]
                else:
                    diff = max(abs(self.__quads__[quad]) - abs(array(point)), 1)
                kill = []
                for i in range(len(diff)):
                    if diff[i] < 0.0:
                        kill.append(i)
                        loki.debug('%s made obsolete by %s' % (repr(self.__quads__[quad][i]), repr(point)))
                kill = sorted(kill, reverse=True)
                loki.debug('kill=' + repr(kill))
                for i in kill:
                    del flat[i]
                    del self.__ids__[quad][i]
                # replace pareto set & append ID
                self.__quads__[quad] = array(flat)
                self.__ids__[quad].append(id)
            return True
        else:
            loki.debug('rejected')
            return False
    def get_set(self, quad=None):
        if quad is None:
            all = []
            ids = []
            for id,quad in zip(self.__ids__, self.__quads__):
                if quad is not None:
                    all.extend(quad.tolist())
                    ids.extend(id)
        elif self.__quads__[quad] is not None:
            all = self.__quads__[quad].tolist()
            ids = array(self.__ids__[quad]).tolist()
        return (ids,all)

if __name__ == '__main__':
    print '-'*7 + ' testing ' + '-'*7
    print 'PearsonCorrelation([1,2,3],[3,4,4]) = %f' % \
            PearsonCorrelation([1,2,3],[3,4,4])
    print 'CornilescuQ([4.5,0.4,-2.5],[4.6,0.4,-2.9]) = %f' % \
            CornilescuQ([4.5,0.4,-2.5],[4.6,0.4,-2.9])

