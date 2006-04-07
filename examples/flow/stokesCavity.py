#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "stokesCavity.py"
 #                                    created: 12/29/03 {3:23:47 PM}
 #                                last update: 4/7/06 {4:54:27 PM}
 # Stolen from:
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-10 JEG 1.0 original
 # ###################################################################
 ##

r"""

This example is an implementation of a rudimentary Stokes solver. It
solves the Navier-Stokes equation in the viscous limit,

.. raw:: latex

    $$ \nabla \mu \cdot \nabla \vec{u} = \nabla p $$ and the
    continuity equation, $$ \nabla \cdot \vec{u} = 0 $$ where $
    \vec{u} $ is the fluid velocity, $ p $ is the pressure and $\mu$
    is the viscosity.  The domain in this example is a square cavity
    of unit dimensions with a moving lid of unit speed.  This example
    uses the SIMPLE algorithm with Rhie-Chow interpolation to solve
    the pressure-momentum coupling. Some of the details of the
    algorithm will be highlighted below but a good reference for this
    material is Ferziger and Peri\'{c}~\cite{ferziger}. The
    solution has a high degree of error close to the corners of the
    domain for the pressure but does a reasonable job of predicting
    the velocities away from the boundaries. A number of aspects of
    \FiPy{} need to be improved to have a first class flow
    solver. These include, higher order spatial diffusion terms,
    proper wall boundary conditions, improved mass flux evaluation and
    extrapolation of cell values to the boundaries using gradients.

In the table below a comparison is made with the Dolfyn_ open source
code on a 100 by 100 grid. The table shows the frequency of values
that fall within the given error confidence bands. Dolfyn_ has the
added features described above. When these features are switched off
the results of Dolfyn_
    
.. _Dolfyn:    http://www.dolfyn.net/

.. raw:: latex

    and \FiPy{} are identical.

    \begin{tabular}{r|rrr}
    \hline
    \% frequency of cells & x-velocity error (\%) & y-velocity error (\%) & pressure error (\%) \\
    \hline
    90                    & $<$ 0.1               & $<$ 0.1               & $<$ 5               \\
    5                     & 0.1 - 0.6             & 0.1 - 0.3             & 5 - 11              \\
    4                     & 0.6 - 7               & 0.3 - 4               & 11 - 35             \\
    1                     & 7 - 96                & 4 - 80                & 35 - 179            \\
    0                     & $>$ 96                & $>$ 80                & $>$ 179             \\
    \hline
    %\caption{The frequency of cell values in \FiPy{} that are within the given error confidence.}
    \end{tabular}

    To start,

some parameters are declared. To run this example to steady
state increase `sweeps` to 200.

    >>> L = 1.0
    >>> N = 50
    >>> dL = L / N
    >>> viscosity = 1.
    >>> pressureRelaxation = 0.2
    >>> velocityRelaxation = 0.5
    >>> sweeps = 5

Build the mesh.

    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(nx = N, ny = N, dx = dL, dy = dL)

Declare the variables.

    >>> from fipy.variables.cellVariable import CellVariable
    >>> pressure = CellVariable(mesh = mesh, name = 'pressure')
    >>> pressureCorrection = CellVariable(mesh = mesh)
    >>> xVelocity = CellVariable(mesh = mesh, name = 'X velocity')
    >>> yVelocity = CellVariable(mesh = mesh, name = 'Y velocity')

The velocity is required as a `VectorFaceVariable` for calculating the
mass flux. This is a somewhat clumsy aspect of the

.. raw:: latex

    \FiPy{}

interface that needs improvement.

    >>> from fipy.variables.vectorFaceVariable import VectorFaceVariable
    >>> velocity = VectorFaceVariable(mesh = mesh)

Build the Stokes equations.

    >>> from fipy.terms.diffusionTerm import DiffusionTerm
    >>> xVelocityEq = DiffusionTerm(viscosity) - pressure.getGrad().dot([1,0])
    >>> yVelocityEq = DiffusionTerm(viscosity) - pressure.getGrad().dot([0,1])
    
.. raw:: latex

    In this example the SIMPLE algorithm is used to couple the
    pressure and momentum equations. Let us assume we have solved the
    discretized momentum equations using a guessed pressure field
    $p^{\ast}$ to obtain a velocity field $\vec{u}^{\ast}$. We would
    like to somehow correct these initial fields to satisfy both the
    discretized momentum and continuity equations. We now try to
    correct these initial fields with a correction such that $\vec{u}
    = \vec{u}^{\ast} + \vec{u}'$ and $p = p^{\ast} + p'$, where
    $\vec{u}$ and $p$ now satisfy the momentum and continuity
    equations. Substituting the exact solution into the equations we
    obtain, $$ \nabla \mu \cdot \nabla \vec{u}' = \vec{p}' $$ and $$
    \nabla \cdot \vec{u}^{\ast} + \nabla \cdot \vec{u}' = 0 $$ We now
    use the discretized form of the equations to write the velocity
    correction in terms of the pressure correction. The discretized
    form of the above equation is, $$ a_P \vec{u}'_P = \sum_f a_A
    \vec{u}'_A - V_P (\nabla p')_P $$ where notation from
    section~\ref{section:linear-equations} is used. The SIMPLE
    algorithm drops the second term in the above equation to leave, $$
    \vec{u}'_{P} = - \frac{ V_P (\nabla p')_P }{ a_P } $$ By
    substituting the above expression into the continuity equations we
    obtain the pressure correction equation, $$ \nabla \frac{V_P}{a_P}
    \cdot \nabla p' = \nabla \cdot \vec{u}^{\ast} $$ In the
    discretized version of the above equation $V_P / a_P$ is
    approximated at the face by $A_f d_{AP} / (a_P)_f.$ In \FiPy{} the
    pressure correction equation can be written as,

..

    >>> ap = CellVariable(mesh = mesh)
    >>> coeff = mesh._getFaceAreas() * mesh._getCellDistances() / ap.getArithmeticFaceValue()
    >>> pressureCorrectionEq = DiffusionTerm(coeff) - velocity.getDivergence()

Set up the no-slip boundary conditions

    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> bcs = (FixedValue(faces = mesh.getFacesLeft(), value = 0),
    ...        FixedValue(faces = mesh.getFacesRight(), value = 0),
    ...        FixedValue(faces = mesh.getFacesBottom(), value = 0),)
    >>> bcsX = bcs + (FixedValue(faces = mesh.getFacesTop(), value = 1),)
    >>> bcsY = bcs + (FixedValue(faces = mesh.getFacesTop(), value = 0),)

Set up the viewers,

    >>> if __name__ == '__main__':
    ...     from fipy import viewers
    ...     viewer = viewers.make(vars = (pressure, xVelocity, yVelocity, velocity))
    
Iterate for a set number of sweeps.

    >>> for sweep in range(sweeps):
    ...
    ...     ## solve the Stokes equations to get starred values
    ...     (xmat, xres) = xVelocityEq.solve(var = xVelocity,
    ...                                      boundaryConditions = bcsX,
    ...                                      returnItems = ['matrix', 'residual'],
    ...                                      underRelaxation = velocityRelaxation)
    ...     (yres,) = yVelocityEq.solve(var = yVelocity,
    ...                                 boundaryConditions = bcsY,
    ...                                 returnItems = ['residual'],
    ...                                 underRelaxation = velocityRelaxation)\
    ...
    ...     ## update the ap coefficient from the matrix diagonal
    ...     ap[:] = -xmat.takeDiagonal()
    ...
    ...     ## update the face velocities based on starred values
    ...     velocity[:,0] = xVelocity.getArithmeticFaceValue()
    ...     velocity[:,1] = yVelocity.getArithmeticFaceValue()
    ...     for id in mesh.getExteriorFaces():
    ...         velocity[id,:] = 0.
    ...
    ...     ## solve the pressure correction equation
    ...     pres, rhs = pressureCorrectionEq.solve(var = pressureCorrection,
    ...                                            returnItems = ['residual', 'RHSvector'])
    ...
    ...     ## update the pressure using the corrected value but hold one cell fixed
    ...     pressure.setValue(pressure + pressureRelaxation * \
    ...                                            (pressureCorrection - pressureCorrection[0]))
    ...     ## update the velocity using the corrected pressure
    ...     xVelocity.setValue(xVelocity - pressureCorrection.getGrad()[:,0] / \
    ...                                                ap * mesh.getCellVolumes())
    ...     yVelocity.setValue(yVelocity - pressureCorrection.getGrad()[:,1] / \
    ...                                                ap * mesh.getCellVolumes())
    ...
    ...     if __name__ == '__main__':
    ...         if sweep%1 == 0:
    ...             print 'sweep:',sweep,', x residual:',max(abs(xres)),', y residual',max(abs(yres)), \
    ...                                  ', p residual:',max(abs(pres)),', continuity:',max(abs(rhs))
    ...
    ... 	viewer.plot()

Test values in the last cell.

    >>> from fipy.tools import numerix
    >>> numerix.allclose(pressure[-1], 145.233883763)
    1
    >>> numerix.allclose(xVelocity[-1], 0.24964673696)
    1
    >>> numerix.allclose(yVelocity[-1], -0.164498041783)
    1

"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript(__name__))
    raw_input('finished')
