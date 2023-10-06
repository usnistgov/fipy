"""
   >>> eq.solve(var,
   ...          boundaryConditions = BCs,
   ...          solver = solver)

   Using the Pysparse solvers, the answer is totally inaccurate. This is due to
   the 4th order term having a high matrix condition number. In this particular
   example, multigrid preconditioners such as those provided by Trilinos allow
   a more accurate solution.

   >>> print(var.allclose(mesh.cellCenters[0], atol = 10))
   1

"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy import input
from fipy import CellVariable, Grid1D, LinearLUSolver, NthOrderBoundaryCondition, DiffusionTerm, Viewer

Lx = 1.
nx = 100000
dx = Lx / nx

mesh = Grid1D(dx = dx, nx = nx)

var = CellVariable(mesh = mesh)

eq = DiffusionTerm((1.0, 1.0))

BCs = (NthOrderBoundaryCondition(mesh.facesLeft, 0., 0),
       NthOrderBoundaryCondition(mesh.facesRight, Lx, 0),
       NthOrderBoundaryCondition(mesh.facesLeft, 0., 2),
       NthOrderBoundaryCondition(mesh.facesRight, 0., 2))

solver = LinearLUSolver(iterations=10)

if __name__ == '__main__':
    eq.solve(var,
             boundaryConditions = BCs,
             solver = solver)

    viewer = Viewer(var)
    viewer.plot()

    print(var.allclose(mesh.cellCenters[0], atol = 10))
    input("finished")

