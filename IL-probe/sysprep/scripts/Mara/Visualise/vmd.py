"""
This module is a wrapper for Scientific.Visualization.VMD to make it work at
my machine... only classes Scene and SceneFile were modified.
"""

from Scientific.IO.TextFile import TextFile
from Scientific.Geometry import Transformation, Vector, VectorModule
import Numeric
import os, string, sys, tempfile
from Scientific.Visualization.Color import *
from Scientific.Visualization.VMD import VMDObject, Molecules, ShapeObject, Sphere, Cube, Cylinder, Cone, Line, Group, isGroup, Arrow, Material, DiffuseMaterial

"""
Produces the VMD files
"""
class SceneFile:

    def __init__(self, filename, mode = 'r', scale = 1., delete = 0):
	if mode == 'r':
	    raise TypeError, 'Not yet implemented.'
	self.file = TextFile(filename, 'w')
	self.memo = {}
	self.delete = delete
	self.scale = scale
	self.filename = filename
	self.writeString('proc mara_graphics {} {\n')
	self.writeString('mol new\n')
	self.writeString('mol rename top mara\n')
	
    def __del__(self):
	self.close()

    def writeString(self, data):
	self.file.write(data)

    def writeVector(self, v):
	self.writeString(" {%g %g %g}" % tuple(v))

    def close(self):
	if self.file is not None:
	    self.writeString('}\nmara_graphics\n')
	    self.writeString('display resetview\n')
	    if self.delete:
		self.writeString('file delete ' + self.filename)
	    self.file.close()
	    self.file = None

    def write(self, object):
	object.writeToFile(self)

"""
The scene
"""
class Scene:

    """VMD scene

    A VMD scene is a collection of graphics objects that can be
    written to a VMD script file or fed directly to VMD.

    Constructor: Scene(|objects|=None, **|options|)

    Arguments:

    |objects| -- a list of graphics objects or 'None' for an empty scene

    |options| -- options as keyword arguments. The only option available
                 is "scale", whose value must be a positive number which
                 specifies a scale factor applied to all coordinates of
                 geometrical objects *except* for molecule objects, which
                 cannot be scaled.
    """

    def __init__(self, objects=None, **options):
	if objects is None:
	    self.objects = []
	elif type(objects) == type([]):
	    self.objects = objects
	else:
	    self.objects = [objects]
	try:
	    self.scale = options['scale']
	except KeyError:
	    self.scale = 1.

    def __len__(self):
	return len(self.objects)

    def __getitem__(self, item):
	return self.object[item]

    def addObject(self, object):
        "Adds |object| to the list of graphics objects."
	self.objects.append(object)

    def writeToFile(self, filename, delete = 0):
        "Writes the scene to a VMD file with name |filename|."
	file = SceneFile(filename, 'w', self.scale, delete)
	for o in self.objects:
	    o.writeToFile(file)
	file.close()

    def view(self):
        "Start VMD for the scene."
	filename = tempfile.mktemp()
	self.writeToFile(filename, 1)
	os.system('vmd -e ' + filename + ' 1> /dev/null 2>&1')


