r"""
Test

>>> print((0.4 < max(q.globalValue[0]) < 0.5))
True

"""
from __future__ import division
from __future__ import unicode_literals

from builtins import range
__docformat__ = 'restructuredtext'

from fipy import input
from fipy import CellVariable, FaceVariable, Grid1D, TransientTerm, CentralDifferenceConvectionTerm
from fipy.tools import numerix

L = 2.
X0 = -1.
nx = 800
cfl = 0.1
K = 4.
rho = 1.

dx = L / nx
m = Grid1D(nx=nx, dx=dx) + X0
x, = m.cellCenters

q = CellVariable(mesh=m, rank=1, elementshape=(2,))

q[0,:] = numerix.exp(-50 * (x - 0.3)**2) * numerix.cos(20 * (x - 0.3))
q[0, x > 0.3] = 0.

Ax = FaceVariable(mesh=m, rank=3, value=[((0, K), (1 / rho, 0))], elementshape=(1, 2, 2))

eqn = TransientTerm() + CentralDifferenceConvectionTerm(Ax) == 0

if  __name__ == '__main__':
    from fipy import MatplotlibViewer as Viewer
    vi = Viewer((q[0], q[1]))
    vi.plot()

for step in range(500):
    eqn.solve(q, dt=cfl * dx)
    if step % 10 ==  0 and  __name__ == '__main__':
        vi.plot()

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    input('finished')

