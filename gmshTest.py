from fipy import *

lines = [
    'Point(1) = { 0.0, 2.0, 0.0 , 1.0};\n',
    'Point(2) = { 2.0, 2.0, 0.0 , 1.0};\n',
    'Point(3) = { 2.0, 0.0, 0.0 , 1.0};\n',
    'Point(4) = { 0.0, 0.0, 0.0 , 1.0};\n',
    '\n',
    'Line(1)  = {1 ,2};\n',
    'Line(2)  = {2, 3};\n',
    'Line(3)  = {3, 4};\n',
    'Line(4)  = {4 ,1};\n',
    '\n',
    'Line Loop(1) = {1, 2, 3, 4};\n',
    '\n',
    'Plane Surface(1) = 1;\n',
    ]
    
import tempfile
(f, geomName) = tempfile.mkstemp('.geo')
file = open(geomName, 'w')
file.writelines(lines)
file.close()
import os
os.close(f)

import sys
if sys.platform == 'win32':
    meshName = 'tmp.msh'
else:
    (f, meshName) = tempfile.mkstemp('.msh')
os.system('gmsh ' + geomName + ' -2 -v 0 -format msh -o ' + meshName)

if sys.platform != 'win32':
    os.close(f)
os.remove(geomName)

from fipy.meshes.gmshImport import GmshImporter2D
mesh3 = GmshImporter2D(meshName)
os.remove(meshName)

var3 = VectorFaceVariable(mesh=mesh3, value=mesh3._getOrientedFaceNormals())
x = mesh3.getCellCenters()[...,0]
var4 = CellVariable(mesh=mesh3, value=x)

TSVViewer(vars=(var4.getFaceGrad(), var3)).plot("Gmsh.txt")