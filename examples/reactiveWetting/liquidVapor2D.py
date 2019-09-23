r"""

A 2D version of the 1D example.

>>> molarWeight = 0.118
>>> ee = -0.455971
>>> gasConstant = 8.314
>>> temperature = 650.
>>> vbar = 1.3e-05

>>> liquidDensity = 7354.3402662299995
>>> vaporDensity = 82.855803327810008

>>> from fipy import CellVariable, Grid2D, TransientTerm, VanLeerConvectionTerm, DiffusionTerm, ImplicitSourceTerm, ConvectionTerm, CentralDifferenceConvectionTerm, Viewer
>>> from fipy.tools import numerix

>>> def f(rho):
...     return ee * rho**2 / molarWeight**2 + gasConstant * temperature * rho / molarWeight * \
...            numerix.log(rho / (molarWeight - vbar * rho))

>>> def mu(rho):
...     return 2 * ee * rho / molarWeight**2 + gasConstant * temperature / molarWeight * \
...            (numerix.log(rho / (molarWeight - vbar * rho)) + molarWeight / (molarWeight - vbar * rho))

>>> Lx = 1e-6
>>> nx = 100
>>> dx = Lx / nx
>>> mesh = Grid2D(nx=nx, ny=nx, dx=dx, dy=dx)

>>> density = CellVariable(mesh=mesh, hasOld=True, name=r'$\rho$')
>>> velocityX = CellVariable(mesh=mesh, hasOld=True, name=r'$u_x$')
>>> velocityY = CellVariable(mesh=mesh, hasOld=True, name=r'$u_y$')
>>> velocityVector = CellVariable(mesh=mesh, name=r'$\vec{u}$', rank=1)
>>> densityPrevious = density.copy()
>>> velocityXPrevious = velocityX.copy()
>>> velocityYPrevious = velocityY.copy()

>>> potentialNC = CellVariable(mesh=mesh, name=r'$\mu^{NC}$')

>>> epsilon = 1e-16
>>> freeEnergy = (f(density) + epsilon * temperature / 2 * density.grad.mag**2).cellVolumeAverage

>>> matrixDiagonal = CellVariable(mesh=mesh, name=r'$a_f$', value=1e+20, hasOld=True)
>>> correctionCoeff = mesh._faceAreas * mesh._cellDistances / matrixDiagonal.faceValue
>>> massEqn = TransientTerm(var=density) \
...           + VanLeerConvectionTerm(coeff=velocityVector.faceValue + correctionCoeff \
...                                         * (density * potentialNC.grad).faceValue, \
...                                   var=density) \
...           - DiffusionTerm(coeff=correctionCoeff * density.faceValue**2, var=potentialNC)

>>> viscosity = 1e-3
>>> ConvectionTerm = CentralDifferenceConvectionTerm
>>> ##matXX = numerix.array([[[2], [0]], [[0], [1]]])
>>> ##matYY = numerix.array([[[1], [0]], [[0], [2]]])
>>> ##matXY = numerix.array([[[0], [0.5]], [[0.5], [0]]])
>>> ##matYX = matXY
>>> matXX = 1
>>> matYY = 1
>>> matXY = 0
>>> matYX = 0
>>> momentumXEqn = TransientTerm(coeff=density, var=velocityX) \
...                + ConvectionTerm(coeff=density.faceValue * velocityVector.faceValue, var=velocityX) \
...                == DiffusionTerm(coeff=(viscosity * matXX,), var=velocityX) \
...                + DiffusionTerm(coeff=(viscosity * matXY,), var=velocityY) \
...                - ConvectionTerm(coeff=density.faceValue * [[1], [0]], var=potentialNC) \
...                + ImplicitSourceTerm(coeff=density.grad[0], var=potentialNC)

>>> momentumYEqn = TransientTerm(coeff=density, var=velocityY) \
...                + ConvectionTerm(coeff=density.faceValue * velocityVector.faceValue, var=velocityY) \
...                == DiffusionTerm(coeff=(viscosity * matYY,), var=velocityY) \
...                + DiffusionTerm(coeff=(viscosity * matYX,), var=velocityX) \
...                - ConvectionTerm(coeff=density.faceValue * [[0], [1]], var=potentialNC) \
...                + ImplicitSourceTerm(coeff=density.grad[1], var=potentialNC)

>>> velocityX.constrain(0, mesh.facesLeft & mesh.facesRight)
>>> velocityY.constrain(0, mesh.facesTop & mesh.facesBottom)

>>> potentialDerivative = 2 * ee / molarWeight**2 + \
...                       gasConstant * temperature * molarWeight / density / (molarWeight - vbar * density)**2

>>> potential = mu(density)

>>> potentialNCEqn = ImplicitSourceTerm(coeff=1, var=potentialNC) \
...                  == potential \
...                  + ImplicitSourceTerm(coeff=potentialDerivative, var=density) \
...                  - potentialDerivative * density \
...                  - DiffusionTerm(coeff=epsilon * temperature, var=density)

>>> potentialNC.faceGrad.constrain(value=[[0], [0]], where=mesh.exteriorFaces)

>>> coupledEqn = massEqn & momentumXEqn & momentumYEqn & potentialNCEqn

>>> numerix.random.seed(2012)
>>> density[:] = (liquidDensity + vaporDensity) / 2 * \
...    (1  + 0.01 * (2 * numerix.random.random(mesh.numberOfCells) - 1))

>>> from fipy import input
>>> if __name__ == '__main__':
...     viewers = Viewer(density), Viewer(velocityVector), Viewer(potentialNC)
...     for viewer in viewers:
...         viewer.plot()
...     input('arrange viewers')
...     for viewer in viewers:
...         viewer.plot()

>>> cfl = 0.1
>>> tolerance = 1e-1
>>> dt = 1e-14
>>> timestep = 0
>>> relaxation = 0.5
>>> sweeps = 0
>>> if __name__ == '__main__':
...     totalSteps = 1e+10
...     totalSweeps = 1e+10
... else:
...     totalSteps = 1
...     totalSweeps = 1

>>> while timestep < totalSteps:
... 
...     sweep = 0
...     dt *= 1.1
...     residual = 1.
...     initialResidual = None
... 
...     density.updateOld()
...     velocityX.updateOld()
...     velocityY.updateOld()
...     matrixDiagonal.updateOld()
... 
...     while residual > tolerance  and sweeps < totalSweeps:
...         sweeps += 1
...         densityPrevious[:] = density
...         velocityXPrevious[:] = velocityX
...         velocityYPrevious[:] = velocityY
...         previousResidual = residual
...         velocityVector[0] = velocityX
...         velocityVector[1] = velocityY
... 
...         dt = min(dt, dx / max(abs(velocityVector.mag)) * cfl)
... 
...         coupledEqn.cacheMatrix()
...         residual = coupledEqn.sweep(dt=dt)
... 
...         if initialResidual is None:
...             initialResidual = residual
... 
...         residual = residual / initialResidual
... 
...         if residual > previousResidual * 1.1 or sweep > 20:
...             density[:] = density.old
...             velocityX[:] = velocityX.old
...             velocityY[:] = velocityY.old
...             matrixDiagonal[:] = matrixDiagonal.old
...             dt = dt / 10.
...             if __name__ == '__main__':
...                 print('Recalculate the time step')
...             timestep -= 1
...             break
...         else:
...             matrixDiagonal[:] = coupledEqn.matrix.takeDiagonal()[mesh.numberOfCells:2 * mesh.numberOfCells]
...             density[:] = relaxation * density + (1 - relaxation) * densityPrevious
...             velocityX[:] = relaxation * velocityX + (1 - relaxation) * velocityXPrevious
...             velocityY[:] = relaxation * velocityY + (1 - relaxation) * velocityYPrevious
... 
...         sweep += 1
... 
...     if __name__ == '__main__' and timestep % 1 == 0:
...         print('timestep: %e / %e, dt: %1.5e, free energy: %1.5e' % (timestep, totalSteps, dt, freeEnergy))
...         for viewer in viewers:
...             viewer.plot()
... 
...     timestep += 1

>>> from fipy import input
>>> if __name__ == '__main__':
...     input('finished')

"""
from __future__ import unicode_literals

__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
