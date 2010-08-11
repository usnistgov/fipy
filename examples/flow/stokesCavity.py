#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "stokesCavity.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: Benny Malengier <bm@cage.ugent.be>
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

This example is an implementation of a rudimentary Stokes solver on a 
colocated grid. 
It solves the Navier-Stokes equation in the viscous limit,

.. math::

   \nabla \mu \cdot \nabla \vec{u} = \nabla p
   
and the
continuity equation,

.. math::
    
   \nabla \cdot \vec{u} = 0 
   
where :math:`\vec{u}` is the fluid velocity, :math:`p` is the pressure and :math:`\mu`
is the viscosity.  The domain in this example is a square cavity
of unit dimensions with a moving lid of unit speed.  This example
uses the SIMPLE algorithm with Rhie-Chow interpolation for colocated grids to solve
the pressure-momentum coupling. Some of the details of the
algorithm will be highlighted below but a good reference for this
material is Ferziger and Peric [ferziger]_ and Rossow [rossow:2003]_. The 
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
>>> #from fipy.meshes.numMesh.grid2D import Grid2D

>>> L = 1.0
>>> N = 50
>>> dL = L / N
>>> viscosity = 1
>>> U = 1.
>>> #0.8 for pressure and 0.5 for velocity are typical relaxation values for SIMPLE
>>> pressureRelaxation = 0.8
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
flux. This is required by the Rhie-Chow correction to avoid pressure/velocity
decoupling.

>>> velocity = FaceVariable(mesh=mesh, rank=1)

Build the Stokes equations in the cell centers.

>>> xVelocityEq = DiffusionTerm(coeff=viscosity) - pressure.getGrad().dot([1.,0.])
>>> yVelocityEq = DiffusionTerm(coeff=viscosity) - pressure.getGrad().dot([0.,1.])
    
In this example the SIMPLE algorithm is used to couple the
pressure and momentum equations. Let us assume we have solved the
discretized momentum equations using a guessed pressure field
:math:`p^{\ast}` to obtain a velocity field :math:`\vec{u}^{\ast}`. That is 
:math:`\vec{u}^{\ast}` is found from

.. math::
   
   a_P \vec{u}^{\ast}_P = \sum_f a_A \vec{u}^{\ast}_A - V_P (\nabla p^{\ast})_P

We would like to somehow correct these initial fields to satisfy both the
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
   
We now use the discretized form of the equations to write the velocity
correction in terms of the pressure correction. The discretized
form of the above equation results in an equation for :math:`p = p'`, 

.. math:: 
    
   a_P \vec{u}'_P = \sum_f a_A \vec{u}'_A - V_P (\nabla p')_P
   
where notation from :ref:`section:linear-equations` is used. The SIMPLE
algorithm drops the second term in the above equation to leave, 

.. math::
    
   \vec{u}'_{P} = - \frac{ V_P (\nabla p')_P }{ a_P }
   
By substituting the above expression into the continuity equations we
obtain the pressure correction equation, 

.. math::
    
   \nabla \frac{V_P}{a_P} \cdot \nabla p' = \nabla \cdot \vec{u}^{\ast}
   
In the discretized version of the above equation :math:`V_P / a_P` is
approximated at the face by :math:`A_f d_{AP} / (a_P)_f`. In :term:`FiPy` the
pressure correction equation can be written as, 

>>> ap = CellVariable(mesh=mesh, value=1.)
>>> coeff = 1./ ap.getArithmeticFaceValue()*mesh._getFaceAreas() * mesh._getCellDistances() 
>>> pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.getDivergence()

Above would work good on a staggered grid, however, on a colocated grid as :term:`FiPy`
uses, the term :term:`velocity.getDivergence()` will cause oscillations in the
pressure solution as velocity is a face variable. 
We can apply the Rhie-Chow correction terms for this. In this an intermediate 
velocity term :math:`u^\Diamond` is considered which does not contain the pressure corrections:

.. math::
    
   \vec{u}^{\Diamond}_P = \vec{u}^{\ast}_P + \frac{V_P}{a_P} (\nabla p^{\ast})_P 
   = \sum_f \frac{a_A}{a_P} \vec{u}^{\ast}_A 

This velocity is interpolated at the edges, after which the pressure correction
term is added again, but now considered at the edge: 

.. math::
    
  \vec{u}_f = \frac{1}{2}(\vec{u}^{\Diamond}_L + \vec{u}^{\Diamond}_R)) 
  - \left(\frac{V}{a_P}\right)_{\mathrm{avg\ L,R}} (\nabla p^{\ast}_f) 

where :math:`\left(\frac{V}{a_P}\right)_{\mathrm{avg\ L,R}}`  is assumed a good approximation at the edge. Here L
and R denote the two cells adjacent to the face. Expanding the 
not calculated terms we arrive at

.. math::
    
  \vec{u}_f = \frac{1}{2}(\vec{u}^{\ast}_L + \vec{u}^{\ast}_R)) 
  + \frac{1}{2}\left(\frac{V}{a_P}\right)_{\mathrm{avg\ L,R}} (\nabla p^{\ast}_L+ \nabla p^{\ast}_R)
  - \left(\frac{V}{a_P}\right)_{\mathrm{avg\ L,R}} (\nabla p^{\ast}_f) 

where we have replaced the coefficients of the cell pressure gradients by 
an averaged value over the edge. 
This formula has the consequence that the velocity on a face depends not only
on the pressure of the adjacent cells, but also on the cells further away, which 
removes the unphysical pressure oscillations. We start by introducing needed
terms 

>>> from fipy.variables.faceGradVariable import _FaceGradVariable
>>> volume = CellVariable(mesh=mesh, value=mesh.getCellVolumes(), name='Volume')
>>> contrvolume=volume.getArithmeticFaceValue()

And set up the velocity with this formula in the SIMPLE loop. 
Now, set up the no-slip boundary conditions

.. index:: FixedValue
   
>>> bcs = (FixedValue(faces=mesh.getFacesLeft(), value=0.),
...        FixedValue(faces=mesh.getFacesRight(), value=0.),
...        FixedValue(faces=mesh.getFacesBottom(), value=0.),)
>>> bcsX = bcs + (FixedValue(faces=mesh.getFacesTop(), value=U ),)
>>> bcsY = bcs + (FixedValue(faces=mesh.getFacesTop(), value=0.),)
>>> bcsPC = (FixedValue(faces=mesh.getFacesLeft() & (mesh.getFaceCenters()[1]<0.9*dL), value=0.),)

Set up the viewers,

.. index::
   :module: fipy.viewers
   
>>> if __name__ == '__main__':
...     viewer = Viewer(vars=(pressure, xVelocity, yVelocity, velocity),
...                xmin=0., xmax=1., ymin=0., ymax=1., colorbar=True)

Below, we iterate for a set number of sweeps. We use the :meth:`sweep`
method instead of :meth:`solve` because we require the residual for
output.  We also use the :meth:`cacheMatrix`, :meth:`getMatrix`,
:meth:`cacheRHSvector` and :meth:`getRHSvector` because both the matrix and
RHS vector are required by the SIMPLE algorithm. Additionally, the
:meth:`sweep` method is passed an ``underRelaxation`` factor to relax the
solution. This argument cannot be passed to :meth:`solve`.


Should one want to use a different solver for the pressure correction, then 
that can be easily done by activatin solverpc.

>>> solverpc = DefaultAsymmetricSolver(tolerance=1e-7)

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
...     ## update the face velocities based on starred values with the 
...     ## Rhie-Chow correction. 
...     xvface = xVelocity.getArithmeticFaceValue()
...     yvface = yVelocity.getArithmeticFaceValue()
...     ## cell pressure gradient
...     presgrad = pressure.getGrad()
...     ## face pressure gradient
...     facepresgrad = _FaceGradVariable(pressure)
...
...     velocity[0] = xVelocity.getArithmeticFaceValue() \
...          + contrvolume / ap.getArithmeticFaceValue() * \
...            (presgrad[0].getArithmeticFaceValue()-facepresgrad[0])
...     velocity[1] = yVelocity.getArithmeticFaceValue() \
...          + contrvolume / ap.getArithmeticFaceValue() * \
...            (presgrad[1].getArithmeticFaceValue()-facepresgrad[1])
...     velocity[..., mesh.getExteriorFaces().getValue()] = 0.
...     velocity[0, mesh.getFacesTop().getValue()] = U
...
...     ## solve the pressure correction equation
...     pressureCorrectionEq.cacheRHSvector()
...     ## left bottom point must remain at pressure 0, so no correction
...     pres = pressureCorrectionEq.sweep(var=pressureCorrection, 
...                                       boundaryConditions=bcsPC,
...                                       #solver=solverpc
...                                       )
...     rhs = pressureCorrectionEq.getRHSvector()
...
...     ## update the pressure using the corrected value
...     pressure.setValue(pressure + pressureRelaxation * pressureCorrection )
...     ## update the velocity using the corrected pressure
...     xVelocity.setValue(xVelocity - pressureCorrection.getGrad()[0] / \
...                                                ap * mesh.getCellVolumes())
...     yVelocity.setValue(yVelocity - pressureCorrection.getGrad()[1] / \
...                                                ap * mesh.getCellVolumes())
...
...     if __name__ == '__main__':
...         if sweep%10 == 0:
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

>>> print numerix.allclose(pressure.getGlobalValue()[...,-1], 162.790867927)
1
>>> print numerix.allclose(xVelocity.getGlobalValue()[...,-1], 0.265072740929)
1
>>> print numerix.allclose(yVelocity.getGlobalValue()[...,-1], -0.150290488304)
1
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript(__name__))
    raw_input('finished')
