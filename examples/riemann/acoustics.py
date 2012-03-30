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
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  
 # ###################################################################
 ##

r"""

>>> from fipy import Grid1D, CellVariable, numerix, RoeConvectionTerm
>>> from fipy import Viewer, TransientTerm

>>> L = 2.
>>> X0 = -1.
>>> nx = 800
>>> cfl = 0.1
>>> K = 4.
>>> rho = 1.

>>> dx = L / nx

>>> m = Grid1D(nx=nx + 4, dx=dx) + (X0 - 2 * dx)
>>> x, = m.cellCenters
>>> X, = m.faceCenters

>>> q = CellVariable(mesh=m, rank=1, elementshape=(2,))

>>> q[0,:] = numerix.exp(-50 * (x - 0.3)**2) * numerix.cos(20 * (x - 0.3))
>>> q[0, x > 0.3] = 0.

>>> Ax = CellVariable(mesh=m, rank=3, value=[((0, K), (1 / rho, 0))], elementshape=(1, 2, 2))

>>> eqn = TransientTerm() + RoeConvectionTerm(Ax) == 0

>>> def setBCs(q):
...     q[0,0] = q[0,3]
...     q[0,1] = q[0,2]
...     q[1,0] = -q[1,3]
...     q[1,1] = -q[1,2]
...     q[:,-1] = 0.
...     q[:,-2] = 0.
 
>>> setBCs(q)
 
>>> if  __name__ == '__main__':
...     from fipy import MatplotlibViewer as Viewer
...     vi = Viewer((q[0], q[1]))
...     vi.plot()
...     timesteps = 10000
...     raw_input('press key')
... else:
...     timesteps = 100
 
>>> elapsedTime = 0.0
>>> dt = 0.9 * dx / eqn.maxeigenvalue(q)
>>> elapsedTime = 0.0
>>> timesteps = 100

>>> for step in range(timesteps):
...     eqn.solve(q, dt=dt)
...     setBCs(q)
...     elapsedTime += dt
...     if step % 20 ==  0 and  __name__ == '__main__':
...         vi.plot()

>>> import os 
>>> filepath = os.path.splitext(__file__)[0] + 'Clawpack100.gz'
>>> print q[:,2:-2].allclose(numerix.loadtxt(filepath, skiprows=6).swapaxes(0,1))
True

>>> timesteps = 600
>>> for step in range(timesteps):
...     eqn.solve(q, dt=dt)
...     setBCs(q)
...     elapsedTime += dt
...     if step % 20 ==  0 and  __name__ == '__main__':
...         vi.plot()

>>> filepath = os.path.splitext(__file__)[0] + 'Clawpack700.gz'
>>> print q[:,2:-2].allclose(numerix.loadtxt(filepath, skiprows=6).swapaxes(0,1))
True

"""

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input('finished')
