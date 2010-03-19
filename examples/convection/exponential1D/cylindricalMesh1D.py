#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "cylindricalMesh1D.py"
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

This example solves the steady-state cylindrical convection-diffusion equation
given by

.. math::

   \nabla \cdot \left(D \nabla \phi + \vec{u} \phi \right) = 0

with coefficients :math:`D = 1` and :math:`\vec{u} = (10,)`, or

>>> diffCoeff = 1.
>>> convCoeff = (10.,)
    
We define a 1D cylindrical mesh representing an anulus

.. index:: Grid1D

>>> from fipy import *

>>> r0 = 1.
>>> r1 = 2.
>>> nr = 100
>>> mesh = CylindricalGrid1D(dr=(r1 - r0) / nr, nr=nr) + ((r0,),)

and impose the boundary conditions

.. math::

   \phi = \begin{cases}
   0& \text{at $r = r_0$,} \\
   1& \text{at $r = r_1$,}
   \end{cases}
   
or

>>> valueLeft = 0.
>>> valueRight = 1.
>>> boundaryConditions = (
...     FixedValue(faces=mesh.getFacesLeft(), value=valueLeft),
...     FixedValue(faces=mesh.getFacesRight(), value=valueRight),
...     )

The solution variable is initialized to ``valueLeft``:
    
>>> var = CellVariable(mesh=mesh, name = "variable")

The equation is created with the :class:`~fipy.terms.diffusionTerm.DiffusionTerm` and
:class:`~fipy.terms.exponentialConvectionTerm.ExponentialConvectionTerm`.

>>> eq = (DiffusionTerm(coeff=diffCoeff)
...       + ExponentialConvectionTerm(coeff=convCoeff))
   
More details of the benefits and drawbacks of each type of convection
term can be found in :ref:`sec:NumericalSchemes`.
Essentially, the :class:`~fipy.terms.exponentialConvectionTerm.ExponentialConvectionTerm` and :class:`~fipy.terms.powerLawConvectionTerm.PowerLawConvectionTerm` will
both handle most types of convection-diffusion cases, with the
:class:`~fipy.terms.powerLawConvectionTerm.PowerLawConvectionTerm` being more efficient.

We solve the equation

>>> eq.solve(var=var, boundaryConditions=boundaryConditions)
   
and test the solution against the analytical result

.. math::

   \phi = \exp{\frac{u}{D} \left(r_1 - r\right)} \left( \frac{ \ei{\frac{u r_0}{D}} - \ei{\frac{u r}{D}} }{ \ei{\frac{u r_0}{D}} - \ei{\frac{u r_1}{D}} } \right)

or

.. index:: exp

>>> axis = 0
>>> try:
...     from scipy.special import expi
...     r = mesh.getCellCenters()[axis]
...     AA = exp(convCoeff[axis] / diffCoeff * (r1 - r))
...     BB = expi(convCoeff[axis] * r0 / diffCoeff) - expi(convCoeff[axis] * r / diffCoeff)
...     CC = expi(convCoeff[axis] * r0 / diffCoeff) - expi(convCoeff[axis] * r1 / diffCoeff)
...     analyticalArray = AA * BB / CC
... except ImportError:
...     print "The SciPy library is unavailable. It is required for testing purposes."

>>> print var.allclose(analyticalArray, atol=1e-3)
1
   
If the problem is run interactively, we can view the result:

.. index::
   module: viewers

>>> if __name__ == '__main__':
...     viewer = viewers.make(vars=var)
...     viewer.plot()
"""
__docformat__ = 'restructuredtext'
     
if __name__ == '__main__':
    import fipy.tests.doctestPlus
    exec(fipy.tests.doctestPlus._getScript())
    
    raw_input('finished')
