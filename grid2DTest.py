from fipy import *

mesh3 = Grid2D(dx=(0.25, 0.75, 1.), dy=(0.1, 1.8, 0.1))

var3 = VectorFaceVariable(mesh=mesh3, value=mesh3._getOrientedFaceNormals())
x = mesh3.getCellCenters()[...,0]
var4 = CellVariable(mesh=mesh3, value=x)

TSVViewer(vars=(var4.getFaceGrad(), var3)).plot("Grid2D.txt")