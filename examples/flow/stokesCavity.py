#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "stokesCavity.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
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
 # ###################################################################
 ##

r"""

This example is an implementation of a rudimentary Stokes solver. It
solves the Navier-Stokes equation in the viscous limit,

.. math::

   \nabla \mu \cdot \nabla \vec{u} = \nabla p
   
and the
continuity equation,

.. math::
    
   \nabla \cdot \vec{u} = 0 
   
where :math:`\vec{u}` is the fluid velocity, :math:`p` is the pressure and :math:`\mu`
is the viscosity.  The domain in this example is a square cavity
of unit dimensions with a moving lid of unit speed.  This example
uses the SIMPLE algorithm with Rhie-Chow interpolation to solve
the pressure-momentum coupling. Some of the details of the
algorithm will be highlighted below but a good reference for this
material is Ferziger and Peric [ferziger]_. The 
solution has a high degree of error close to the corners of the
domain for the pressure but does a reasonable job of predicting
the velocities away from the boundaries. A number of aspects of
:term:`FiPy` need to be improved to have a first class flow
solver. These include, higher order spatial diffusion terms,
proper wall boundary conditions, improved mass flux evaluation and
extrapolation of cell values to the boundaries using gradients.

In the table below a comparison is made with the Dolfyn_ open source
code on a 100 by 100 grid. The table shows the frequency of values
that fall within the given error confidence bands. Dolfyn_ has the
added features described above. When these features are switched off
the results of Dolfyn_
and :term:`FiPy` are identical.

.. _Dolfyn:    http://www.dolfyn.net/

.. tabularcolumns:: r|rrr
    
=====================  =====================  =====================  ===================
\% frequency of cells  x-velocity error (\%)  y-velocity error (\%)  pressure error (\%)
=====================  =====================  =====================  ===================
90                     :math:`<0.1`           :math:`< 0.1`          :math:`< 5`        
5                      0.1 to 0.6             0.1 to 0.3             5 to 11            
4                      0.6 to 7               0.3 to 4               11 to 35           
1                      7 to 96                4 to 80                35 to 179          
0                      :math:`> 96`           :math:`> 80`           :math:`> 179`      
=====================  =====================  =====================  ===================



To start, some parameters are declared.

>>> from fipy import *

>>> L = 1.0
>>> N = 50
>>> dL = L / N
>>> viscosity = 1.
>>> pressureRelaxation = 0.2
>>> velocityRelaxation = 0.5
>>> if __name__ == '__main__':
...     sweeps = 300
... else:
...     sweeps = 5

Build the mesh.

.. index:: Grid2D
   
>>> mesh = Grid2D(nx=N, ny=N, dx=dL, dy=dL)

Declare the variables.

.. index:: CellVariable
   
>>> pressure = CellVariable(mesh=mesh, name='pressure')
>>> pressureCorrection = CellVariable(mesh=mesh)
>>> xVelocity = CellVariable(mesh=mesh, name='X velocity')
>>> yVelocity = CellVariable(mesh=mesh, name='Y velocity')

The velocity is required as a rank-1
:class:`~fipy.variables.faceVariable.FaceVariable` for calculating the mass
flux. This is a somewhat clumsy aspect of the :term:`FiPy`
interface that needs improvement.

>>> velocity = FaceVariable(mesh=mesh, rank=1)

Build the Stokes equations.

>>> xVelocityEq = DiffusionTerm(coeff=viscosity) - pressure.getGrad().dot([1,0])
>>> yVelocityEq = DiffusionTerm(coeff=viscosity) - pressure.getGrad().dot([0,1])
    
In this example the SIMPLE algorithm is used to couple the
pressure and momentum equations. Let us assume we have solved the
discretized momentum equations using a guessed pressure field
:math:`p^{\ast}` to obtain a velocity field :math:`\vec{u}^{\ast}`. We would
like to somehow correct these initial fields to satisfy both the
discretized momentum and continuity equations. We now try to
correct these initial fields with a correction such that 
:math:`\vec{u} = \vec{u}^{\ast} + \vec{u}'` and :math:`p = p^{\ast} + p'`, where
:math:`\vec{u}` and :math:`p` now satisfy the momentum and continuity
equations. Substituting the exact solution into the equations we
obtain, 

.. math:: 
    
   \nabla \mu \cdot \nabla \vec{u}' = \vec{p}' 
   
and 

.. math::
   
   \nabla \cdot \vec{u}^{\ast} + \nabla \cdot \vec{u}' = 0 
   
We now
use the discretized form of the equations to write the velocity
correction in terms of the pressure correction. The discretized
form of the above equation is, 

.. math:: 
    
   a_P \vec{u}'_P = \sum_f a_A \vec{u}'_A - V_P (\nabla p')_P
   
where notation from :ref:`section:linear-equations` is used. The SIMPLE
algorithm drops the second term in the above equation to leave, 

.. math::
    
   \vec{u}'_{P} = - \frac{ V_P (\nabla p')_P }{ a_P }
   
By
substituting the above expression into the continuity equations we
obtain the pressure correction equation, 

.. math::
    
   \nabla \frac{V_P}{a_P} \cdot \nabla p' = \nabla \cdot \vec{u}^{\ast}
   
In the
discretized version of the above equation :math:`V_P / a_P` is
approximated at the face by :math:`A_f d_{AP} / (a_P)_f`. In :term:`FiPy` the
pressure correction equation can be written as, 

>>> ap = CellVariable(mesh=mesh)
>>> coeff = mesh._getFaceAreas() * mesh._getCellDistances() / ap.getArithmeticFaceValue()
>>> pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.getDivergence()

Set up the no-slip boundary conditions

.. index:: FixedValue
   
>>> bcs = (FixedValue(faces=mesh.getFacesLeft(), value=0),
...        FixedValue(faces=mesh.getFacesRight(), value=0),
...        FixedValue(faces=mesh.getFacesBottom(), value=0),)
>>> bcsX = bcs + (FixedValue(faces=mesh.getFacesTop(), value=1),)
>>> bcsY = bcs + (FixedValue(faces=mesh.getFacesTop(), value=0),)

Set up the viewers,

.. index::
   module: fipy.viewers
   
>>> if __name__ == '__main__':
...     viewer = Viewer(vars=(pressure, xVelocity, yVelocity, velocity))

Below, we iterate for a set number of sweeps. We use the :meth:`sweep`
method instead of :meth:`solve` because we require the residual for
output.  We also use the :meth:`cacheMatrix`, :meth:`getMatrix`,
:meth:`cacheRHSvector` and :meth:`getRHSvector` because both the matrix and
RHS vector are required by the SIMPLE algorithm. Additionally, the
:meth:`sweep` method is passed an ``underRelaxation`` factor to relax the
solution. This argument cannot be passed to :meth:`solve`.

.. index:: sweep, cacheMatrix, getMatrix, cacheRHSvector, getRHSvector
   
>>> for sweep in range(sweeps):
...
...     ## solve the Stokes equations to get starred values
...     xVelocityEq.cacheMatrix()
...     xres = xVelocityEq.sweep(var=xVelocity,
...                              boundaryConditions=bcsX,
...                              underRelaxation=velocityRelaxation)
...     xmat = xVelocityEq.getMatrix()
...
...     yres = yVelocityEq.sweep(var=yVelocity,
...                              boundaryConditions=bcsY,
...                              underRelaxation=velocityRelaxation)
...
...     ## update the ap coefficient from the matrix diagonal
...     ap[:] = -xmat.takeDiagonal()
...
...     ## update the face velocities based on starred values
...     velocity[0] = xVelocity.getArithmeticFaceValue()
...     velocity[1] = yVelocity.getArithmeticFaceValue()
...     velocity[..., mesh.getExteriorFaces().getValue()] = 0.
...
...     ## solve the pressure correction equation
...     pressureCorrectionEq.cacheRHSvector()
...     pres = pressureCorrectionEq.sweep(var=pressureCorrection)
...     rhs = pressureCorrectionEq.getRHSvector()
...
...     ## update the pressure using the corrected value but hold one cell fixed
...     pressure.setValue(pressure + pressureRelaxation * \
...                                            (pressureCorrection - pressureCorrection.getGlobalValue()[0]))
...     ## update the velocity using the corrected pressure
...     xVelocity.setValue(xVelocity - pressureCorrection.getGrad()[0] / \
...                                                ap * mesh.getCellVolumes())
...     yVelocity.setValue(yVelocity - pressureCorrection.getGrad()[1] / \
...                                                ap * mesh.getCellVolumes())
...
...     if __name__ == '__main__':
...         if sweep%1 == 0:
...             print 'sweep:',sweep,', x residual:',xres, \
...                                  ', y residual',yres, \
...                                  ', p residual:',pres, \
...                                  ', continuity:',max(abs(rhs))
...
...             viewer.plot()

.. image:: cavity.*
   :width: 90%
   :align: center

Test values in the last cell.

>>> print numerix.allclose(pressure.getGlobalValue()[...,-1], 145.233883763)
1
>>> print numerix.allclose(xVelocity.getGlobalValue()[...,-1], 0.24964673696)
1
>>> print numerix.allclose(yVelocity.getGlobalValue()[...,-1], -0.164498041783)
1
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript(__name__))
    raw_input('finished')
