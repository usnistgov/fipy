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

|FiPy| 2.0 introduces several syntax changes from |FiPy| 1.0. We appreciate that
this is very inconvenient for our users, but we hope you'll agree that the new
sytax is easier to read and easier to use. We assure you that this is not
something we do casually; it has been over three years since our last
incompatible change (when |FiPy| 1.0 superceded |FiPy| 0.1).

All examples included with version 2.0 have been updated to use the new syntax,
but any scripts you have written for |FiPy| 1.0 will need to be updated. A
complete listing of the changes needed to take the |FiPy| examples scripts from
version 1.0 to version 2.0 can be found at
    
    http://www.matforge.org/fipy/wiki/upgrade1_0examplesTo2_0
    
but we summarize the necessary changes here. If these tips are not sufficient to
make your scripts compatible with |FiPy| 2.0, please don't hesitate to ask for
help on the `mailing list`_.


The following items **must** be changed in your scripts

 * The dimension axis of a `Variable` is now first, not last
   
      >>> x = mesh.getCellCenters()[0]

   instead of
   
      >>> x = mesh.getCellCenters()[...,0]
       
   This seemingly arbitrary change simplifies a great many things in |FiPy|, but
   the one most noticeable to user is that you can now write
   
       >>> x, y = mesh.getCellCenters()
  
   instead of
   
       >>> x = mesh.getCellCenters()[...,0]
       >>> y = mesh.getCellCenters()[...,1]

 * `VectorCellVariable` and `VectorFaceVariable` no longer exist. `CellVariable`
   and `FaceVariable` now both inherit from `MeshVariable`, and a `MeshVariable`
   can have arbitrary rank. A field of scalars (default) will have `rank=0`, a
   field of vectors will have `rank=1`, etc. You should write
   
       >>> vectorField = CellVariable(mesh=mesh, rank=1)

   instead of

       >>> vectorField = VectorCellVariable(mesh=mesh)

 * For internal reasons, |FiPy| now supports `CellVariable` and `FaceVariable`
   objects that contain integers, but it is not meaningful to solve a PDE for an
   integer field (|FiPy| should issue a warning if you try). As a result, when given,
   initial values must be specified as floating-point values:
 
       >>> var = CellVariable(mesh=mesh, value=1.)
 
   where they used to be quietly accepted as integers
   
       >>> var = CellVariable(mesh=mesh, value=1)
   
   If the `value` argument is not supplied, the `CellVariable` will contain
   floats, as before.
       
 * The `faces` argument to `BoundaryCondition` now takes a mask, instead of a
   list of `Face` IDs. Now you write
   
       >>> FixedValue(faces=mesh.getExteriorFaces() & (X**2 < 1e-6), value=...)
       
   instead of
   
       >>> FixedValue(faces=mesh.getExteriorFaces().where(X**2 < 1e-6), value=...)
       
   This simple case may not seem like much of an improvement, but with the old
   syntax, it was extremely difficult to achieve the effect of
   
       >>> FixedValue(faces=(mesh.getExteriorFaces() 
       ...                   & (((X**2 < 1e-6) 
       ...                       & (Y > 3.)) 
       ...                      | (phi.getArithmeticFaceValue() 
       ...                         < sin(gamma.getArithmeticFaceValue())))), value=...)

   although it probably could have been done with a rather convoluted (and
   slow!) `filter` function passed to `where`.

The remaining changes are not **required**, but they make scripts easier to read
and we recommend them. |FiPy| may issue a `DeprecationWarning` for some cases,
to indicate that we may not maintain the old syntax indefinitely.

 * All of the most commonly used classes and functions in |FiPy| are directly 
   accessible in the `fipy` namespace. For brevity, our examples now start with

       >>> from fipy import *

   instead of the explicit

       >>> from fipy.meshes.grid1D import Grid1D 
       >>> from fipy.terms.powerLawConvectionTerm import PowerLawConvectionTerm 
       >>> from fipy.variables.cellVariable import CellVariable 

   imports that we used to use. The explicit imports should continue to work, so
   you do not need to change them if you don't wish to, but we find our own scripts
   much easier to read without them.
   
   All of the `numerix` module is now imported into the `fipy` namespace, so you
   can call `numerix` functions a number of different ways, including:
       
       >>> from fipy import *
       >>> y = exp(x)

   or
   
       >>> from fipy import numerix
       >>> y = numerix.exp(x)
   
   or
   
       >>> from fipy.tools.numerix import exp
       >>> y = exp(x)
       
   We generally use the first, but you may see the others, as well, and should
   feel free to use whichever form you find most comfortable.

   .. note::
       
      Internally, |FiPy| uses explicit imports, as is considered 
      `best Python practice`_, but we feel that clarity trumps orthodoxy when it
      comes to the examples.

   .. _best Python practice: http://docs.python.org/howto/doanddont.html#from-module-import

 * The function `fipy.viewers.make()` has been renamed to
   `fipy.viewers.Viewer()`. All of the `limits` can now be supplied as direct
   arguments, as well (although this is not required). The result is a more
   natural syntax:
        
       >>> from fipy import Viewer
       >>> viewer = Viewer(vars=(alpha, beta, gamma), datamin=0, datamax=1)

   instead of

       >>> from fipy import viewers
       >>> viewer = viewers.make(vars=(alpha, beta, gamma), 
       ...                       limits={'datamin': 0, 'datamax': 1})

   With the old syntax, there was also a temptation to write

       >>> from fipy.viewers import make
       >>> viewer = make(vars=(alpha, beta, gamma))

   which can be very hard to understand after the fact (`make`? `make` what?).
   
 * A `ConvectionTerm` can now calculate its Peclet number automatically, so the
   `diffusionTerm` argument is no longer required

       >>> eq = (TransientTerm()  
       ...       == ImplicitDiffusionTerm(coeff=diffCoeff) 
       ...       + PowerLawConvectionTerm(coeff=convCoeff)) 

   instead of 
   
       >>> diffTerm = ImplicitDiffusionTerm(coeff=diffCoeff) 
       >>> eq = (TransientTerm() 
       ...       == diffTerm 
       ...       + PowerLawConvectionTerm(coeff=convCoeff, diffusionTerm=diffTerm)) 

 * An `ImplicitSourceTerm` now "knows" how to partion itself onto the solution
   matrix, so you can write

       >>> S0 = mXi * phase * (1 - phase) - phase * S1 
       >>> source = S0 + ImplicitSourceTerm(coeff=S1) 
   
   instead of
   
       >>> S0 = mXi * phase * (1 - phase) - phase * S1 * (S1 < 0) 
       >>> source = S0 + ImplicitSourceTerm(coeff=S1 * (S1 < 0)) 
   
   It is definitely still advantageous to hand-linearize your source terms, but
   it is no longer necessary to worry about putting the "wrong" sign on the
   diagonal of the matrix.


.. |FiPy| raw:: latex

   \FiPy{}
   
.. _mailing list:         http://www.ctcms.nist.gov/fipy/mail.html
"""
__docformat__ = 'restructuredtext'

