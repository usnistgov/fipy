"""
Demonstration of my new VTKViewer.raw() method.
"""

from fipy import Grid3D, CellVariable, VTKViewer
import pyvista as pv

# FiPy setup
mesh = Grid3D(dx=1, dy=1, dz=1, nx=10, ny=10, nz=10)
var = CellVariable(mesh)
var.setValue(1, where=mesh.cellCenters[0] > 5)
var.setValue(2, where=mesh.cellCenters[0] <= 5)
viewer = VTKViewer(vars=var)

# pyvista plot
vtk = viewer.raw()
pv.wrap(vtk).plot()