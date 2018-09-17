#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "acoustics.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

r"""
Test

>>> print (0.4 < max(q.globalValue[0]) < 0.5)
True

"""

__docformat__ = 'restructuredtext'

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

    raw_input('finished')
