r"""

This example demonstrates the use of the :class:`~fipy.terms.vanLeerConvectionTerm.VanLeerConvectionTerm` as
defined by http://www.gre.ac.uk/~physica/phy2.12/theory/node173.htm

In this example a square wave is advected. The Van Leer discretization
should in theory do a good job of preserving the shape of the
wave. This may or may not be happening in this case. This example
needs further testing.

The test case is mainly to check that the periodic mesh is working
correctly. We advect the wave on different meshes one periodic and one
non-periodic but twice as long. The results are then compared. The
periodic wave wraps around the mesh.

>>> from builtins import range
>>> for step in range(steps):
...     eq1.solve(var=var1, dt=dt, solver=DefaultAsymmetricSolver(tolerance=1e-11, iterations=10000))
...     eq2.solve(var=var2, dt=dt, solver=DefaultAsymmetricSolver(tolerance=1e-11, iterations=10000))

>>> print(numerix.allclose(var1.globalValue[nx // 2:3 * nx // 4],
...                        var2.globalValue[:nx // 4], atol=1e-6))
1

Currently after 20 steps the wave has lost 23% of its height. Van Leer
should do better than this.

>>> print(var1.max() > 0.77)
1
"""
from __future__ import print_function
from __future__ import division
from __future__ import unicode_literals

from builtins import range
__docformat__ = 'restructuredtext'

from fipy import input
from fipy import CellVariable, Grid1D, PeriodicGrid1D, TransientTerm, VanLeerConvectionTerm, DefaultAsymmetricSolver, Viewer
from fipy.tools import numerix

L = 20.
nx = 40
dx = L / nx
cfl = 0.5
velocity = 1.0
dt = cfl * dx / velocity

steps = int(L /  4. / dt / velocity)

mesh = Grid1D(dx = dx, nx = nx)

periodicMesh = PeriodicGrid1D(dx=dx, nx=nx // 2)

startingArray = numerix.zeros(nx, 'd')
startingArray[2 * nx // 10: 3 * nx // 10] = 1.

var1 = CellVariable(
    name = "non-periodic",
    mesh = mesh,
    value = startingArray)

var2 = CellVariable(
    name = "periodic",
    mesh = periodicMesh,
    value = startingArray[:nx // 2])

eq1 = TransientTerm() - VanLeerConvectionTerm(coeff = (-velocity,))
eq2 = TransientTerm() - VanLeerConvectionTerm(coeff = (-velocity,))

if __name__ == '__main__':

    viewer1 = Viewer(vars=var1)
    viewer2 = Viewer(vars=var2)
    viewer1.plot()
    viewer2.plot()

    newVar2 = var2.copy()

    for step in range(steps):
        eq1.solve(var=var1, dt=dt, solver=DefaultAsymmetricSolver())
        eq2.solve(var=var2, dt=dt, solver=DefaultAsymmetricSolver())
        viewer1.plot()
        viewer2.plot()

    newVar2[:nx // 4] = var2[nx // 4:]
    newVar2[nx // 4:] = var2[:nx // 4]

    print('maximum absolute difference between periodic and non-periodic grids:', abs(var1[nx // 4:(3 * nx) // 4] - newVar2).max())

    input('finished')
