#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 8/8/05 {10:00:07 AM} { 5:14:21 PM}
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
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-17 JEG 1.0 original
 # ###################################################################
 ##

r"""

In this example we solve a coupled phase and temperature equation to model 
solidification, and eventually dendritic growth, from a circular seed in a 2D mesh:
    
    >>> numberOfCells = 40
    >>> Length = numberOfCells * 2.5 / 100.
    >>> nx = numberOfCells
    >>> ny = numberOfCells
    >>> dx = Length / nx
    >>> dy = Length / ny
    >>> radius = Length / 4.
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx,dy,nx,ny)
    
Dendritic growth will not be observed with this small test system. If
you wish to see dendritic growth reset the following parameters:
``numberOfCells = 200``, ``steps = 10000``, ``radius = Length / 80``.

The governing equation for the phase field is given by:

.. raw:: latex

    $$ \tau_{\phi} \frac{\partial \phi}{\partial t} 
    = \nabla \cdot \left[ D \nabla \phi + A \nabla \xi \right] +
      \phi ( 1 - \phi ) m ( \phi , T) $$

where

.. raw:: latex

    $$ m(\phi, T) 
    = \phi - \frac{1}{2} - \frac{ \kappa_1 }{ \pi } \arctan \left( \kappa_2 T \right). $$

    The coefficients $D$ and $A$ are given by,

    $$ D = \alpha^2 \left[ 1 + c \beta \right]^2 $$

    and

    $$ A = \alpha^2 c \left[ 1 + c \beta \right] \Phi_\psi $$

    where $ \beta = \frac{ 1 - \Phi^2 } { 1 + \Phi^2} $,
    $ \Phi = \tan \left( \frac{ \theta } { 2 } + \frac{ N } { 2 } \arctan \psi \right) $,
    $ \psi = \frac{ \phi_y } { \phi_x } $ and
    $ \xi_x = -\phi_y $ and $ \xi_y = \phi_x $.

    The governing equation for temperature is given by:

    $$ \frac{\partial T}{\partial t} = D_T \nabla^2 T + \frac{\partial \phi}{\partial t} $$

..  Further details of the numerical method for this problem can be found in
    "Extending Phase Field Models of Solidification to Polycrystalline
    Materials", J.A. Warren *et al.*, *Acta Materialia*, **51** (2003) 6035-6058.

Here the phase and temperature equations are solved with an explicit
and implicit technique, respectively.

The parameters for these equations are 

    >>> timeStepDuration = 5e-5
    >>> tau = 3e-4
    >>> alpha = 0.015
    >>> c = 0.02
    >>> N = 4.
    >>> kappa1 = 0.9
    >>> kappa2 = 20.    
    >>> tempDiffusionCoeff = 2.25
    >>> theta = 0

The `phase` variable is `0` for a liquid and `1` for a solid.  Here we
build an example `phase` variable, initialized as a liquid,

    >>> from fipy.variables.cellVariable import CellVariable
    >>> phase = CellVariable(
    ...     name = 'phaseField',
    ...     mesh = mesh,
    ...     value = 0.,
    ...     hasOld = 1)

The `hasOld` flag keeps the old value of the variable. This is
necessary for a transient solution. In this example we wish to set up
an interior region that is solid. A value of `1` is assigned to the
`phase` variable on a patch defined by the method:

    >>> def circleCells(cell,L = Length):
    ...     x = cell.getCenter()
    ...     r = radius
    ...     c = (L / 2., L / 2.)
    ...     if (x[0] - c[0])**2 + (x[1] - c[1])**2 < r**2:
    ...         return 1
    ...     else:
    ...         return 0
   
This method is passed to `mesh.getCells(filter = circleCells)` which
filters out the required cells.
   
    >>> interiorCells = mesh.getCells(filter = circleCells)           
    >>> phase.setValue(1.,interiorCells)

The temperature field is initialized to a value of `-0.4` throughout:

    >>> temperature = CellVariable(
    ...     name = 'temperature',
    ...     mesh = mesh,
    ...     value = -0.4,
    ...     hasOld = 1
    ...     )

.. raw:: latex
  
   The $m(\phi, T)$ is created from the `phase` and `temperature` variables.

..

    >>> from fipy.tools import numerix
    >>> mVar = phase - 0.5 - kappa1 / numerix.pi * \
    ...     numerix.arctan(kappa2 * temperature)

.. raw:: latex

    The following section of code builds up the $A$ and $D$ coefficients.

..

    >>> dPhiy = phase.getFaceGrad().dot((0, 1))
    >>> dPhix = phase.getFaceGrad().dot((1, 0))
    >>> arc = N * numerix.arctan2(dPhiy, dPhix) + theta
    >>> Phi = numerix.tan(arc / 2)
    >>> PhiSq = Phi**2
    >>> beta = (1. - PhiSq) / (1. + PhiSq)
    >>> dbdpsi = -N * 2 * Phi / (1 + PhiSq)
    >>> A = alpha**2 * c * (1.+ c * beta) * dbdpsi
    >>> D = alpha**2 * (1.+ c * beta)**2

    >>> dxi = phase.getFaceGrad()._take((1, 0), axis = 1) * (-1, 1)
    >>> anisotropySource = (A * dxi).getDivergence()
        
    >>> from fipy.terms.transientTerm import TransientTerm
    >>> from fipy.terms.explicitDiffusionTerm import ExplicitDiffusionTerm
    >>> from fipy.terms.implicitSourceTerm import ImplicitSourceTerm
    >>> phaseEq = TransientTerm(tau) == ExplicitDiffusionTerm(D) + \
    ...     ImplicitSourceTerm(mVar * ((mVar < 0) - phase)) + \
    ...     ((mVar > 0.) * mVar * phase + anisotropySource)

The temperature equation is built in the following way,

    >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
    >>> temperatureEq = TransientTerm() == \
    ...                 ImplicitDiffusionTerm(tempDiffusionCoeff) + \
    ...                 (phase - phase.getOld()) / timeStepDuration

If we are running this example interactively, we create viewers for
the phase and temperature fields

    >>> if __name__ == '__main__':
    ...     import fipy.viewers
    ...     phaseViewer = fipy.viewers.make(vars = phase)
    ...     temperatureViewer = fipy.viewers.make(vars = temperature,
    ...         limits = {'datamin': -0.5, 'datamax': 0.5})
    ...     phaseViewer.plot()
    ...     temperatureViewer.plot()

we iterate the solution in time, plotting as we go if running interactively,

    >>> steps = 10
    >>> for i in range(steps):
    ...     phase.updateOld()
    ...     temperature.updateOld()
    ...     phaseEq.solve(phase, dt = timeStepDuration)
    ...     temperatureEq.solve(temperature, dt = timeStepDuration)
    ...     if i%10 == 0 and __name__ == '__main__':
    ...         phaseViewer.plot()
    ...         temperatureViewer.plot()

The solution is compared with test data. The test data was created for
``steps = 10`` with a FORTRAN code written by Ryo Kobayashi for phase
field modeling. The following code opens the file ``test.gz`` extracts
the data and compares it with the `phase` variable.

   >>> import os
   >>> testFile = 'test.gz'
   >>> import examples.phase.anisotropy
   >>> import gzip
   >>> filepath = os.path.join(examples.phase.anisotropy.__path__[0], testFile)
   >>> filestream = gzip.open(filepath,'r')
   >>> import cPickle
   >>> testData = cPickle.load(filestream)
   >>> filestream.close()
   >>> phase =  numerix.array(phase)
   >>> testData = numerix.reshape(testData, phase.shape)
   >>> print testData.allclose(phase, rtol = 1e-10, atol = 1e-10)
   1
   
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input('finished')

