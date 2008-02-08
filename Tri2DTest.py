from fipy import *

mesh3 = Tri2D(nx=2,ny=2)

var3 = VectorFaceVariable(mesh=mesh3, value=mesh3._getOrientedFaceNormals())
x = mesh3.getCellCenters()[...,0]
var4 = CellVariable(mesh=mesh3, value=x)

TSVViewer(vars=(var4.getFaceGrad(), var3)).plot("Tri2D.txt")