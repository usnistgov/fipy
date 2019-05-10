# -*- coding: utf-8 -*-


r"""
How to update scripts from version 1.0 to 2.0.

:term:`FiPy` 2.0 introduces several syntax changes from :term:`FiPy` 1.0. We appreciate that
this is very inconvenient for our users, but we hope you'll agree that the new
syntax is easier to read and easier to use. We assure you that this is not
something we do casually; it has been over three years since our last
incompatible change (when :term:`FiPy` 1.0 superceded :term:`FiPy` 0.1).

All examples included with version 2.0 have been updated to use the new syntax,
but any scripts you have written for :term:`FiPy` 1.0 will need to be updated. A
complete listing of the changes needed to take the :term:`FiPy` examples scripts from
version 1.0 to version 2.0 can be found with:

    $ git diff version-1_2 version-2_0 examples/

but we summarize the necessary changes here. If these tips are not sufficient to
make your scripts compatible with :term:`FiPy` 2.0, please don't hesitate to ask for
help on the `mailing list`_.


The following items **must** be changed in your scripts

 * The dimension axis of a :class:`~fipy.variables.variable.Variable` is now first, not last

   >>> x = mesh.getCellCenters()[0]

   instead of

   >>> x = mesh.getCellCenters()[..., 0]

   This seemingly arbitrary change simplifies a great many things in :term:`FiPy`, but
   the one most noticeable to the user is that you can now write

   >>> x, y = mesh.getCellCenters()

   instead of

   >>> x = mesh.getCellCenters()[..., 0]
   >>> y = mesh.getCellCenters()[..., 1]

   Unfortunately, we cannot reliably automate this conversion, but we find that
   searching for "``...,``" and "``:,``" finds almost everything. Please don't
   blindly "search & replace all" as that is almost bound to create more
   problems than it's worth.

   .. note::

      Any vector constants must be reoriented. For instance, in order to offset
      a :class:`~fipy.meshes.mesh.Mesh`, you must write

      >>> mesh = Grid2D(...) + ((deltax,), (deltay,))

      or

      >>> mesh = Grid2D(...) + [[deltax], [deltay]]

      instead of

      >>> mesh = Grid2D(...) + (deltax, deltay)

 * :class:`~fipy.variables.vectorCellVariable.VectorCellVariable` and
   :class:`~fipy.variables.vectorFaceVariable.VectorFaceVariable` no longer exist.
   :class:`~fipy.variables.cellVariable.CellVariable` and
   and :class:`~fipy.variables.faceVariable.FaceVariable` now both inherit from
   :class:`~fipy.varibles.meshVariable.MeshVariable`, which
   can have arbitrary rank. A field of scalars (default) will have ``rank=0``, a
   field of vectors will have ``rank=1``, etc. You should write

   >>> vectorField = CellVariable(mesh=mesh, rank=1)

   instead of

   >>> vectorField = VectorCellVariable(mesh=mesh)

   .. note::

      Because vector fields are properly supported, use vector operations to
      manipulate them, such as

      >>> phase.getFaceGrad().dot((( 0, 1),
      ...                          (-1, 0)))

      instead of the hackish

      >>> phase.getFaceGrad()._take((1, 0), axis=1) * (-1, 1)

 * For internal reasons, :term:`FiPy` now supports
   :class:`~fipy.variables.cellVariable.CellVariable` and
   :class:`~fipy.variables.faceVariable.FaceVariable`
   objects that contain integers, but it is not meaningful to solve a PDE for an
   integer field (:term:`FiPy` should issue a warning if you try). As a result, when given,
   initial values must be specified as floating-point values:

   >>> var = CellVariable(mesh=mesh, value=1.)

   where they used to be quietly accepted as integers

   >>> var = CellVariable(mesh=mesh, value=1)

   If the ``value`` argument is not supplied, the
   :class:`~fipy.variables.cellVariable.CellVariable` will contain floats, as before.

 * The ``faces`` argument to
   :class:`~fipy.boundaryConditions.boundaryCondition.BoundaryCondition` now takes a mask,
   instead of a list of :class:`~fipy.meshes.face.Face` IDs. Now you write

   >>> X, Y = mesh.getFaceCenters()
   >>> FixedValue(faces=mesh.getExteriorFaces() & (X**2 < 1e-6), value=...)

   instead of

   >>> exteriorFaces = mesh.getExteriorFaces()
   >>> X = exteriorFaces.getCenters()[..., 0]
   >>> FixedValue(faces=exteriorFaces.where(X**2 < 1e-6), value=...)

   With the old syntax, a different call to
   :meth:`~fipy.meshes.face.Face.getCenters` had to be made for each set
   of :class:`~fipy.meshes.face.Face` objects. It was also extremely
   difficult to specify boundary conditions that depended both on position in
   space and on the current values of any other :class:`~fipy.variables.variable.Variable`.

   >>> FixedValue(faces=(mesh.getExteriorFaces()
   ...                   & (((X**2 < 1e-6)
   ...                       & (Y > 3.))
   ...                      | (phi.getArithmeticFaceValue()
   ...                         < sin(gamma.getArithmeticFaceValue())))), value=...)

   although it probably could have been done with a rather convoluted (and
   slow!) ``filter`` function passed to ``where``. There no longer are any ``filter``
   methods used in :term:`FiPy`. You now would write

   >>> x, y = mesh.cellCenters
   >>> initialArray[(x < dx) | (x > (Lx - dx)) | (y < dy) | (y > (Ly - dy))] = 1.

   instead of the *much* slower

   >>> def cellFilter(cell):
   ...     return ((cell.center[0] < dx)
   ...             or (cell.center[0] > (Lx - dx))
   ...             or (cell.center[1] < dy)
   ...             or (cell.center[1] > (Ly - dy)))

   >>> positiveCells = mesh.getCells(filter=cellFilter)
   >>> for cell in positiveCells:
   ...     initialArray[cell.ID] = 1.

   Although they still exist, we find very little cause to ever call
   :meth:`~fipy.meshes.mesh.Mesh.getCells`
   or :meth:`fipy.meshes.mesh.Mesh.getFaces`.

 * Some modules, such as :mod:`fipy.solvers`, have been significantly rearranged.
   For example, you need to change

   >>> from fipy.solvers.linearPCGSolver import LinearPCGSolver

   to either

   >>> from fipy import LinearPCGSolver

   or

   >>> from fipy.solvers.pysparse.linearPCGSolver import LinearPCGSolver


 * The :func:`numerix.max` and :func:`numerix.min` functions no longer exist. Either
   call :func:`max` and :func:`min` or the :meth:`~fipy.variables.variable.Variable.max`
   and :meth:`~fipy.variables.variable.Variable.min` methods of a
   :class:`~fipy.variables.variable.Variable`.

 * The :mod:`Numeric` module has not been supported for a long time. Be sure to use

   >>> from fipy import numerix

   instead of

   >>> import Numeric

The remaining changes are not *required*, but they make scripts easier to read
and we recommend them. :term:`FiPy` may issue a :exc:`DeprecationWarning` for some cases,
to indicate that we may not maintain the old syntax indefinitely.

 * All of the most commonly used classes and functions in :term:`FiPy` are directly
   accessible in the :mod:`fipy` namespace. For brevity, our examples now start with

   >>> from fipy import *

   instead of the explicit

   >>> from fipy.meshes.grid1D import Grid1D
   >>> from fipy.terms.powerLawConvectionTerm import PowerLawConvectionTerm
   >>> from fipy.variables.cellVariable import CellVariable

   imports that we used to use. Most of the explicit imports should continue to
   work, so you do not need to change them if you don't wish to, but we find our
   own scripts much easier to read without them.

   All of the :mod:`~fipy.tools.numerix` module is now imported into the :mod:`fipy` namespace, so you
   can call :mod:`~fipy.tools.numerix` functions a number of different ways, including:

   >>> from fipy import *
   >>> y = exp(x)

   or

   >>> from fipy import numerix
   >>> y = numerix.exp(x)

   or

   >>> from fipy.tools.numerix import exp
   >>> y = exp(x)

   We generally use the first, but you may see us use the others, and should
   feel free to use whichever form you find most comfortable.

   .. note::

      Internally, :term:`FiPy` uses explicit imports, as is considered
      `best Python practice`_, but we feel that clarity trumps orthodoxy when it
      comes to the examples.

   .. _best Python practice: http://docs.python.org/howto/doanddont.html#from-module-import

 * The function :func:`fipy.viewers.make` has been renamed to
   :func:`fipy.viewers.Viewer`. All of the ``limits`` can now be supplied as direct
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

   which can be very hard to understand after the fact (``make``? ``make`` what?).

 * A :class:`~fipy.terms.convectionTerm.ConvectionTerm` can now calculate its
   Péclet number automatically, so the ``diffusionTerm`` argument is no longer required

   >>> eq = (TransientTerm()
   ...       == DiffusionTerm(coeff=diffCoeff)
   ...       + PowerLawConvectionTerm(coeff=convCoeff))

   instead of

   >>> diffTerm = DiffusionTerm(coeff=diffCoeff)
   >>> eq = (TransientTerm()
   ...       == diffTerm
   ...       + PowerLawConvectionTerm(coeff=convCoeff, diffusionTerm=diffTerm))

 * An :class:`~fipy.terms.term.implicitSourceTerm.ImplicitSourceTerm` now "knows"
   how to partition itself onto the solution matrix, so you can write

   >>> S0 = mXi * phase * (1 - phase) - phase * S1
   >>> source = S0 + ImplicitSourceTerm(coeff=S1)

   instead of

   >>> S0 = mXi * phase * (1 - phase) - phase * S1 * (S1 < 0)
   >>> source = S0 + ImplicitSourceTerm(coeff=S1 * (S1 < 0))

   It is definitely still advantageous to hand-linearize your source terms, but
   it is no longer necessary to worry about putting the "wrong" sign on the
   diagonal of the matrix.

 * To make clearer the distinction between iterations, timesteps, and sweeps
   (see FAQ :ref:`FAQ-IterationsTimestepsSweeps`)
   the ``steps`` argument to a :class:`~fipy.solvers.solver.Solver` object has
   been renamed ``iterations``.

 * :class:`~fipy.terms.implicitDiffusionTerm.ImplicitDiffusionTerm` has been
   renamed to :class:`~fipy.terms.diffusionTerm.DiffusionTerm`.

.. _mailing list:         http://www.ctcms.nist.gov/fipy/mail.html
"""
from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

