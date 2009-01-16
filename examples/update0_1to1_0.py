#!/usr/bin/env python

## 
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
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

The syntax of |FiPy| scripts has changed considerably between version 0.1
and version 1.0.  We have incremented the full version-number to stress
that previous scripts are incompatible.  We strongly believe that these
changes are for the better, resulting in easier code to write and read as
well as slightly improved efficiency, but we realize that this represents
an inconvenience to our users that have already written scripts of their
own.  We will strive to avoid any such incompatible changes in the future.

All examples included with version 1.0 have been updated to use the new
syntax, but any scripts you have written for |FiPy| 0.1 will need to be
updated to work with |FiPy| 1.0.  As a tutorial for updating your scripts,
we will walk through updating the file
`examples/convection/exponential1D/input.py` from |FiPy| 0.1.  If you
attempt to run that script with |FiPy| 1.0, the script will fail and you
will see the errors shown below:

----

This example solves the steady-state convection-diffusion equation
given by:

.. raw:: latex

     $$ \nabla \cdot \left(D \nabla \phi + \vec{u} \phi \right) = 0 $$

with coefficients

.. raw:: latex

   $D = 1$ and $\vec{u} = (10, 0)$,
   
or

    >>> diffCoeff = 1.
    >>> convCoeff = (10.,0.)
    
We define a 1D mesh

.. raw:: latex

   \IndexClass{Grid2D}
   
..

    >>> L = 10.
    >>> nx = 1000
    >>> ny = 1
    >>> from fipy.meshes.grid2D import Grid2D
    >>> mesh = Grid2D(L / nx, L / ny, nx, ny)

and impose the boundary conditions

.. raw:: latex

   $$ \phi = \begin{cases}
   0& \text{at $x = 0$,} \\
   1& \text{at $x = L$,}
   \end{cases} $$ 
   or
   \IndexClass{FixedValue}
   \IndexClass{FixedFlux}
      
..


    >>> valueLeft = 0.
    >>> valueRight = 1.
    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> from fipy.boundaryConditions.fixedFlux import FixedFlux
    >>> boundaryConditions = (
    ...     FixedValue(mesh.getFacesLeft(), valueLeft),
    ...     FixedValue(mesh.getFacesRight(), valueRight),
    ...     FixedFlux(mesh.getFacesTop(), 0.),
    ...     FixedFlux(mesh.getFacesBottom(), 0.)
    ...     )

The solution variable is initialized to `valueLeft`:
    
.. raw:: latex

   \IndexClass{CellVariable}

..

    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(
    ...     name = "concentration",
    ...     mesh = mesh,
    ...     value = valueLeft)

The `SteadyConvectionDiffusionScEquation` object is
used to create the equation.  It needs to be passed a convection term
instantiator as follows:

.. raw:: latex

   \IndexClass{ExponentialConvectionTerm}
   \IndexClass{LinearLUSolver}
   \IndexClass{SteadyConvectionDiffusionScEquation}

..

   >>> from fipy.terms.exponentialConvectionTerm import ExponentialConvectionTerm
   >>> from fipy.solvers import *
   >>> from fipy.equations.stdyConvDiffScEquation import SteadyConvectionDiffusionScEquation
   Traceback (most recent call last):
   ...
   ImportError: No module named equations.stdyConvDiffScEquation
   >>> eq = SteadyConvectionDiffusionScEquation(
   ...      var = var,
   ...      diffusionCoeff = diffCoeff,
   ...      convectionCoeff = convCoeff,
   ...      solver = LinearLUSolver(tolerance = 1.e-15, steps = 2000),
   ...      convectionScheme = ExponentialConvectionTerm,
   ...      boundaryConditions = boundaryConditions
   ...      )
   Traceback (most recent call last):
   ...
   NameError: name 'SteadyConvectionDiffusionScEquation' is not defined

More details of the benefits and drawbacks of each type of convection
term can be found in the numerical section of the manual. Essentially
the `ExponentialConvectionTerm` and `PowerLawConvectionTerm` will both
handle most types of convection diffusion cases with the
`PowerLawConvectionTerm` being more efficient.

We iterate to equilibrium

.. raw:: latex

   \IndexClass{Iterator}

..

    >>> from fipy.iterators.iterator import Iterator
    >>> it = Iterator((eq,))
    Traceback (most recent call last):
    ...
    NameError: name 'eq' is not defined
    >>> it.timestep()
    Traceback (most recent call last):
    ...
    NameError: name 'it' is not defined

and test the solution against the analytical result

.. raw:: latex

   $$ \phi = \frac{1 - \exp(-u_x x / D)}{1 - \exp(-u_x L / D)} $$
   or 
   \IndexSoftware{numerix}

..

    >>> axis = 0
    >>> x = mesh.getCellCenters()[:,axis]
    >>> from fipy.tools import numerix
    >>> CC = 1. - numerix.exp(-convCoeff[axis] * x / diffCoeff)
    >>> DD = 1. - numerix.exp(-convCoeff[axis] * L / diffCoeff)
    >>> analyticalArray = CC / DD
    >>> numerix.allclose(analyticalArray, var, rtol = 1e-10, atol = 1e-10)
    0
   
If the problem is run interactively, we can view the result:

.. raw:: latex

   \IndexClass{Grid2DGistViewer}

..

    >>> if __name__ == '__main__':
    ...     from fipy.viewers.grid2DGistViewer import Grid2DGistViewer
    Traceback (most recent call last):
    ...
    ImportError: No module named grid2DGistViewer
    
::
    
    ...     viewer = Grid2DGistViewer(var)
    ...     viewer.plot()
    
----

We see that a number of errors are thrown:
    
    - `ImportError: No module named equations.stdyConvDiffScEquation`
    - `NameError: name 'SteadyConvectionDiffusionScEquation' is not defined`
    - `NameError: name 'eq' is not defined`
    - `NameError: name 'it' is not defined`
    - `ImportError: No module named grid2DGistViewer`
    
As is usually the case with computer programming, many of these errors are
caused by earlier errors.  Let us update the script, section by
section:

Although no error was generated by the use of `Grid2D`, |FiPy| 1.0 supports
a true 1D mesh class, so we instantiate the mesh as

.. raw:: latex

   \IndexClass{Grid1D}

..

    >>> L = 10.
    >>> nx = 1000
    >>> from fipy.meshes.grid1D import Grid1D
    >>> mesh = Grid1D(dx = L / nx, nx = nx)
    
The `Grid2D` class with `ny = 1` still works perfectly well for 1D
problems, but the `Grid1D` class is slightly more efficient, and it makes
the code clearer when a 1D geometry is actually desired.

Because the mesh is now 1D, we must update the convection coefficient
vector to be 1D as well

    >>> diffCoeff = 1.
    >>> convCoeff = (10.,)
    
The `FixedValue` boundary conditions at the left and right are unchanged,
but a `Grid1D` mesh does not even have top and bottom faces:

.. raw:: latex

   \IndexClass{FixedValue}

..

    >>> valueLeft = 0.
    >>> valueRight = 1.
    >>> from fipy.boundaryConditions.fixedValue import FixedValue
    >>> boundaryConditions = (
    ...     FixedValue(mesh.getFacesLeft(), valueLeft),
    ...     FixedValue(mesh.getFacesRight(), valueRight))

The creation of the solution variable is unchanged:
    
.. raw:: latex

   \IndexClass{CellVariable}

..

    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(name = "concentration",
    ...                    mesh = mesh,
    ...                    value = valueLeft)

The biggest change between |FiPy| 0.1 and |FiPy| 1.0 is that `Equation`
objects no longer exist at all.  Instead, `Term` objects can be simply
added, subtracted, and equated to assemble an equation.  Where before the
assembly of the equation occurred in the black-box of
`SteadyConvectionDiffusionScEquation`, we now assemble it directly:

.. raw:: latex

   \IndexClass{ImplicitDiffusionTerm}
   \IndexClass{ExponentialConvectionTerm}

..

    >>> from fipy.terms.implicitDiffusionTerm import ImplicitDiffusionTerm
    >>> diffTerm = ImplicitDiffusionTerm(coeff = diffCoeff)
    
    >>> from fipy.terms.exponentialConvectionTerm import ExponentialConvectionTerm
    >>> eq = diffTerm + ExponentialConvectionTerm(coeff = convCoeff, 
    ...                                           diffusionTerm = diffTerm)
    
One thing that `SteadyConvectionDiffusionScEquation` took care of
automatically was that a `ConvectionTerm` must know about any
`DiffusionTerm` in the equation in order to calculate a Peclet number.
Now, the `DiffusionTerm` must be explicitly passed to the `ConvectionTerm`
in the `diffusionTerm` parameter.

The `Iterator` class still exists, but it is no longer necessary.  Instead,
the solution to an implicit steady-state problem like this can simply be
obtained by telling the equation to solve itself (with an appropriate
`solver` if desired, although the default `LinearPCGSolver` is usually
suitable):
    
.. raw:: latex

   \IndexClass{LinearLUSolver}

..

    >>> from fipy.solvers import *
    >>> eq.solve(var = var, 
    ...          solver = LinearLUSolver(tolerance = 1.e-15, steps = 2000), 
    ...          boundaryConditions = boundaryConditions)

.. note:: In version 0.1, the `Equation` object had to be told about the
   `Variable`, `Solver`, and `BoundaryCondition` objects when it was created
   (and it, in turn, passed much of this information to the `Term` objects in
   order to create them).  In version 1.0, the `Term` objects (and the
   equation assembled from them) are abstract.  The `Variable`, `Solver`, and
   `BoundaryCondition` objects are only needed by the `solve()` method (and,
   in fact, the same equation could be used to solve different variables, with
   different solvers, subject to different boundary conditions, if desired).

The analytical solution is unchanged, and we can test as before

    >>> numerix.allclose(analyticalArray, var, rtol = 1e-10, atol = 1e-10)
    1
    
or we can use the slightly simpler syntax

    >>> print var.allclose(analyticalArray, rtol = 1e-10, atol = 1e-10)
    1

The `ImportError: No module named grid2DGistViewer` results because the
`Viewer` classes have been moved and renamed.  This error could be resolved
by changing the `import` statement appropriately:

.. raw:: latex

   \IndexClass{Gist1DViewer}

..

    >>> if __name__ == '__main__':
    ...     from fipy.viewers.gistViewer.gist1DViewer import Gist1DViewer
    ...     viewer = Gist1DViewer(vars = var)
    ...     viewer.plot()

Instead, rather than instantiating a particular `Viewer` (which you can
still do, if you desire), a generic "factory" method will return a `Viewer`
appropriate for the supplied `Variable` object(s):
    
.. raw:: latex

   \IndexModule{viewers}

..

    >>> if __name__ == '__main__':
    ...     import fipy.viewers
    ...     viewer = fipy.viewers.make(vars = var)
    ...     viewer.plot()

Please do not hesitate to contact us if this example does not help you
convert your existing scripts to |FiPy| 1.0.

.. |FiPy| raw:: latex

   \FiPy{}
"""
__docformat__ = 'restructuredtext'

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
    raw_input('finished')
