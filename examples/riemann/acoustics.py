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
Test

>>> print (0.4 < max(q.globalValue[0]) < 0.5)
True

"""

__docformat__ = 'restructuredtext'

from fipy import *

L = 2.
X0 = -1.
nx = 800
cfl = 0.1
K = 4.
rho = 1.

dx = L / nx

dt = Variable(1.)

m = Grid1D(nx=nx, dx=dx) + X0
x, = m.cellCenters
X, = m.faceCenters

q = CellVariable(mesh=m, rank=1, elementshape=(2,))

q[0,:] = numerix.exp(-50 * (x - 0.3)**2) * numerix.cos(20 * (x - 0.3))
q[0, x > 0.3] = 0.

q.faceGrad.constrain(0, [(X == -1) | (X == 1), X == 1])
q.constrain(0, [X == -10, X == -1])

Ax = CellVariable(mesh=m, rank=3, value=[((0, K), (1 / rho, 0))], elementshape=(1, 2, 2))

roeConvectionTerm = RoeConvectionTerm(Ax, dt=dt)
eqn = TransientTerm() + roeConvectionTerm == 0

if  __name__ == '__main__':
    from fipy import MatplotlibViewer as Viewer
    vi = Viewer((q[0], q[1]))
    vi.plot()
    raw_input('press key')

from profiler import Profiler
from profiler import calibrate_profiler


elapsedTime = 0.0
dt.setValue(0.1 * dx / roeConvectionTerm.maxeigenvalue(q))

##fudge = calibrate_profiler(10000)
##profile = Profiler('profile', fudge=fudge)

for step in range(10000):
    eqn.solve(q, dt=float(dt))
    print 'max(max(q[0]), max(q[1]))',max(max(q[0]), max(q[1]))
    elapsedTime += float(dt)
    if step % 100 ==  0 and  __name__ == '__main__':
        vi.plot()
        print 'step',step
        print 'elapsedTime',elapsedTime
##        raw_input('press key')
##profile.stop()

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input('finished')
