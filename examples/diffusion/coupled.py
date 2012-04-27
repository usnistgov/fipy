#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "coupled.py"
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
How to use coupled equations.


## Coupled


from fipy import Grid1D, CellVariable, TransientTerm, DiffusionTerm, Viewer

m = Grid1D(nx=100, Lx=1.)

v0 = CellVariable(mesh=m, hasOld=True, value=0.5)
v1 = CellVariable(mesh=m, hasOld=True, value=0.5)

v0.constrain(0, m.facesLeft)
v0.constrain(1, m.facesRight)

v1.constrain(1, m.facesLeft)
v1.constrain(0, m.facesRight)

eqn0 = TransientTerm(var=v0) == DiffusionTerm(-1, var=v1) + DiffusionTerm(0.01, var=v0)
eqn1 = TransientTerm(var=v1) == DiffusionTerm(1, var=v0) + DiffusionTerm(0.01, var=v1)

eqn = eqn0 & eqn1

vi = Viewer((v0, v1))

for t in range(1): 
    v0.updateOld()
    v1.updateOld()
    eqn.solve(dt=1.)
    vi.plot()

## Vector

v = CellVariable(mesh=m, hasOld=True, value=0.5, elementshape=(2,))

v.constrain(0, [m.facesLeft, m.facesRight])
v.constrain(1, [m.facesRight, m.facesLeft])

eqn = TransientTerm(((1, 0), (0, 1))) == DiffusionTerm([((0.01, -1), (1, 0.01))])

vi = Viewer((v[0], v[1]))

for t in range(1): 
    v.updateOld()
    eqn.solve(var=v, dt=1.)
    vi.plot()

## Uncoupled

m = Grid1D(nx=100, Lx=1.)

v0 = CellVariable(mesh=m, hasOld=True, value=0.5)
v1 = CellVariable(mesh=m, hasOld=True, value=0.5)

v0.constrain(0, m.facesLeft)
v0.constrain(1, m.facesRight)

v1.constrain(1, m.facesLeft)
v1.constrain(0, m.facesRight)

eq0 = TransientTerm() == -v1.faceGrad.divergence + DiffusionTerm(0.01)
eq1 = TransientTerm() == v0.faceGrad.divergence + DiffusionTerm(0.01)

vi = Viewer((v0, v1))

for t in range(10000): 
    v0.updateOld()
    v1.updateOld()
    eq0.solve(var=v0, dt=1e-5)
    eq1.solve(var=v1, dt=1e-5)
    vi.plot()

"""
