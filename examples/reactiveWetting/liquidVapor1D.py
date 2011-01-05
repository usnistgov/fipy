#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "liquidVapor1D.py"
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

>>> from fipy import *

>>> factor = 1

>>> Lx = 1.
>>> nx = 100
>>> dx = Lx / nx

>>> viscosity = 1e-3
>>> epsilon = 1e-2

>>> temperature = 1.854236326217e-01
>>> vbar = 8.149038895253e-01
>>> alpha = 7.0

>>> cfl = 0.1
>>> tolerance = 1e-3
>>> variance = 0.01
>>> dt = 1e-5

>>> liquidDensity = 1.0
>>> vaporDensity = 0.01
>>> averageDensity = (liquidDensity + vaporDensity) / 2.0

>>> mesh = Grid1D(nx=nx, dx=dx)

>>> density = CellVariable(mesh=mesh, hasOld=True, name='density')
>>> density[:] = averageDensity * (1  + variance * (2 * numerix.random.random(mesh.getNumberOfCells()) - 1))

>>> velocity = CellVariable(mesh=mesh, hasOld=True, name='velocity')
>>> velocity.constrain(0, mesh.getExteriorFaces())

>>> potential = CellVariable(mesh=mesh, hasOld=False, name='potential')
>>> ap = CellVariable(mesh=mesh, value=1e+20)
>>> faceDensity = density.getArithmeticFaceValue()

>>> ConvectionTerm = CentralDifferenceConvectionTerm

>>> potentialDerivative = -2. * temperature / density / (1 - vbar * density)**2
>>> potentialConst = temperature * (numerix.log(density / (1 - vbar * density)) + alpha + vbar * density / (1 - vbar * density)**2)

>>> apCoeff =  mesh._getFaceAreas() * mesh._getCellDistances() / ap.getArithmeticFaceValue()
>>> densityCoeff = apCoeff * faceDensity

>>> tensorAverage = (density * potential.getLeastSquaresGrad()).getArithmeticFaceValue()

>>> velocityCorrection = apCoeff * (-tensorAverage)

>>> densityEqn = TransientTerm(var=density) + VanLeerConvectionTerm(coeff=velocity.getFaceValue() - velocityCorrection, var=density) \
             + DiffusionTerm(coeff=-densityCoeff * faceDensity, var=potential)

velocityEqn = TransientTerm(coeff=density, var=velocity) + ConvectionTerm(coeff=[[1]] * faceDensity * velocity.getFaceValue(), var=velocity) \
              - DiffusionTerm(coeff=2 * viscosity, var=velocity) \
              + ConvectionTerm(coeff=density * [[1]], var=potential) - ImplicitSourceTerm(coeff=density.getGrad()[0], var=potential)

>>> potentialEqn = - ImplicitSourceTerm(coeff=potentialDerivative, var=density) + DiffusionTerm(coeff=epsilon, var=density) \
               + ImplicitSourceTerm(1., var=potential) - potentialConst

>>> coupledEqn = densityEqn & velocityEqn & potentialEqn

>>> viewers = Viewer(density), Viewer(velocity), Viewer(potential)
>>> for viewer in viewers:
...     viewer.plot()

>>> for timestep in xrange(1000):
... 
...     residual = 1.
...     dt *= 1.1
...     
...     density.updateOld()
...     velocity.updateOld()
... 
...     sweep = 0
... 
...     while residual > tolerance and sweep < 100:
...
... 
... 
...         dt = min(dt * 1.1, dx / max(abs(velocity)) * cfl)
... 
...         print 'sweep',sweep,'residual',residual, 'dt',dt
... 
...         coupledEqn.cacheMatrix()
...         residual = coupledEqn.sweep(dt=dt)
...         ap[:] = coupledEqn.getMatrix().takeDiagonal()[mesh.getNumberOfCells(): 2 * mesh.getNumberOfCells()]
... 
... ##        print 'residual',residual
... ##        raw_input('stopped')
...         if sweep == 0:
...             initialResidual = residual
... 
...         residual = residual / initialResidual
... 
...         sweep += 1
... 
...     for viewer in viewers:
...         viewer.plot()


"""

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
