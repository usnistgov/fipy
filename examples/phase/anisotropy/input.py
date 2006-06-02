#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "input.py"
 #                                    created: 11/17/03 {10:29:10 AM} 
 #                                last update: 6/2/06 {1:53:15 PM} { 5:14:21 PM}
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
solidification, and eventually dendritic growth, based on the work of
Warren, Kobayashi, Lobkovsky and Carter
    
.. raw:: latex

   \cite{WarrenPolycrystal}
   
   We start from a circular seed in a 2D mesh:
   \IndexClass{Grid2D}

..

    >>> numberOfCells = 40
    >>> Length = numberOfCells * 2.5 / 100.
    >>> nx = numberOfCells
    >>> ny = numberOfCells
    >>> dx = Length / nx
    >>> dy = Length / ny
    >>> radius = Length / 4.
    >>> seedCenter = (Length / 2., Length / 2.)
    >>> initialTemperature = -0.4
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)
    
Dendritic growth will not be observed with this small test system. If
you wish to see dendritic growth reset the following parameters such
that ``numberOfCells = 500``, ``steps = 10000``, ``radius = dx * 5.``
``seedCenter = (0. , 0.)`` and ``initialTemperature = -0.5``.

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

    $$ A = \alpha^2 c \left[ 1 + c \beta \right] \beta_\psi $$

    where $ \beta = \frac{ 1 - \Phi^2 } { 1 + \Phi^2} $,
    $ \Phi = \tan \left( \frac{ N } { 2 } \psi \right) $,
    $ \psi = \theta + \arctan \frac{ \phi_y } { \phi_x } $ and
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
    >>> theta = 0.

The `phase` variable is `0` for a liquid and `1` for a solid.  Here,
the `phase` variable is initialized as a liquid,

.. raw:: latex

   \IndexClass{CellVariable}

..

    >>> from fipy.variables.cellVariable import CellVariable
    >>> phase = CellVariable(name='phase field', mesh=mesh, hasOld=1)

The `hasOld` flag keeps the old value of the variable. This is
necessary for a transient solution. In this example we wish to set up
an interior region that is solid. A value of `1` is assigned to the
`phase` variable on a patch defined by the method:

The domain is seeded with a circular solidified region with parameters
`seedCenter` and `radius` representing the center and radius of the
seed.
   
    >>> x, y = mesh.getCellCenters()[...,0], mesh.getCellCenters()[...,1]
    >>> phase.setValue(1., where=((x - seedCenter[0])**2 
    ...                           + (y - seedCenter[1])**2) < radius**2)

The temperature field is initialized to a value of `-0.4` throughout:

    >>> temperature = CellVariable(
    ...     name='temperature',
    ...     mesh=mesh,
    ...     value=initialTemperature,
    ...     hasOld=1
    ...     )

.. raw:: latex
  
   The $m(\phi, T)$ variable 
   
is created from the `phase` and `temperature` variables.

.. raw:: latex

   \IndexModule{numerix}
   \IndexConstant{\pi}{pi}
   \IndexFunction{arctan}
   \IndexFunction{arctan2}
   \IndexFunction{tan}

..

    >>> from fipy.tools.numerix import pi, arctan, arctan2, tan
    >>> mVar = phase - 0.5 - kappa1 / pi * arctan(kappa2 * temperature)

.. raw:: latex

    The following section of code builds up the $A$ and $D$ coefficients.

..

    >>> phaseY = phase.getFaceGrad().dot((0, 1))
    >>> phaseX = phase.getFaceGrad().dot((1, 0))
    >>> psi = theta + arctan2(phaseY, phaseX)
    >>> Phi = tan(N * psi / 2)
    >>> PhiSq = Phi**2
    >>> beta = (1. - PhiSq) / (1. + PhiSq)
    >>> betaPsi = -N * 2 * Phi / (1 + PhiSq)
    >>> A = alpha**2 * c * (1.+ c * beta) * betaPsi
    >>> D = alpha**2 * (1.+ c * beta)**2

.. raw:: latex

    The $\nabla \xi$ variable

(`dxi`),

.. raw:: latex

    given by $(\xi_x, \xi_y) = (-\phi_y, \phi_x)$,
    is constructed by first obtaining $\nabla \phi$
    
using `getFaceGrad()`. The axes are then swapped using
`_take((1,0))` and finally the vector is multiplied by `(-1, 1)` to
negate the x component.

    >>> dxi = phase.getFaceGrad()._take((1, 0), axis = 1) * (-1, 1)
    >>> anisotropySource = (A * dxi).getDivergence()

The phase equation can now be constructed.
    
.. raw:: latex

   \IndexClass{TransientTerm}
   \IndexClass{ExplicitDiffusionTerm}
   \IndexClass{ImplicitSourceTerm}

..

    >>> from fipy.terms.transientTerm import TransientTerm
    >>> from fipy.terms.explicitDiffusionTerm import ExplicitDiffusionTerm
    >>> from fipy.terms.implicitSourceTerm import ImplicitSourceTerm
    >>> phaseEq = TransientTerm(tau) == ExplicitDiffusionTerm(D) + \
    ...     ImplicitSourceTerm(mVar * ((mVar < 0) - phase)) + \
    ...     ((mVar > 0.) * mVar * phase + anisotropySource)

The temperature equation is built in the following way,

.. raw:: latex

   \IndexClass{ImplicitDiffusionTerm}

..

    >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
    >>> temperatureEq = TransientTerm() == \
    ...                 ImplicitDiffusionTerm(tempDiffusionCoeff) + \
    ...                 (phase - phase.getOld()) / timeStepDuration

If we are running this example interactively, we create viewers for
the phase and temperature fields

.. raw:: latex

   \IndexModule{viewers}

..

    >>> if __name__ == '__main__':
    ...     from fipy.viewers import make
    ...     phaseViewer = make(vars=phase)
    ...     temperatureViewer = make(vars=temperature,
    ...                              limits={'datamin': -0.5, 'datamax': 0.5})
    ...     phaseViewer.plot()
    ...     temperatureViewer.plot()

we iterate the solution in time, plotting as we go if running interactively,

    >>> steps = 10
    >>> for i in range(steps):
    ...     phase.updateOld()
    ...     temperature.updateOld()
    ...     phaseEq.solve(phase, dt=timeStepDuration)
    ...     temperatureEq.solve(temperature, dt=timeStepDuration)
    ...     if i%10 == 0 and __name__ == '__main__':
    ...         phaseViewer.plot()
    ...         temperatureViewer.plot()

The solution is compared with test data. The test data was created for
``steps = 10`` with a FORTRAN code written by Ryo Kobayashi for phase
field modeling. The following code opens the file ``test.gz`` extracts
the data and compares it with the `phase` variable.

.. raw:: latex

   \IndexModule{dump}
   \IndexModule{numerix}
   \IndexFunction{allclose}

..

   >>> import examples.phase.anisotropy
   >>> import os
   >>> filepath = os.path.join(examples.phase.anisotropy.__path__[0], 'test.gz')
   >>> from fipy.tools import dump
   >>> testData = dump.read(filepath)
   >>> from fipy.tools.numerix import allclose
   >>> print allclose(phase, testData)
   1
   
"""
__docformat__ = 'restructuredtext'

if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())

    raw_input('finished')

