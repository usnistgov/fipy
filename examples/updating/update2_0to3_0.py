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
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

r"""
How to update scripts from version 2.0 to 3.0.

:term:`FiPy` 3.0 introduces several syntax changes from :term:`FiPy` 2.0. We appreciate that
this is very inconvenient for our users, but we hope you'll agree that the new
syntax is easier to read and easier to use. We assure you that this is not
something we do casually; it has been over two and a half years since our last
incompatible change (when :term:`FiPy` 2.0 superceded :term:`FiPy` 1.0).

All examples included with version 3.0 have been updated to use the new syntax,
but any scripts you have written for :term:`FiPy` 2.0 will need to be updated. A
complete listing of the changes needed to take the :term:`FiPy` examples scripts from
version 2.0 to version 3.0 can be found at

    http://www.matforge.org/fipy/wiki/upgrade2_0examplesTo3_0

but we summarize the necessary changes here. If these tips are not sufficient to
make your scripts compatible with :term:`FiPy` 3.0, please don't hesitate to ask for
help on the `mailing list`_.


The following items **must** be changed in your scripts

 * We have reconsidered the change in FiPy 2.0 that included all of the
   functions of the :mod:`~fipy.tools.numerix` module in the :mod:`fipy` namespace.
   You now must be more explicit when referring to any of these functions:

   >>> from fipy import *
   >>> y = numerix.exp(x)

   >>> from fipy.tools.numerix import exp
   >>> y = exp(x)

   We generally use the first, but you may see us import specific functions if
   we feel it improves readability. You should feel free to use whichever form
   you find most comfortable.

   .. note:: the old behavior can be obtained, at least for now, by setting
             the :envvar:`FIPY_INCLUDE_NUMERIX_ALL` environment variable.

 * If your equation contains a :class:`~fipy.terms.transientTerm.TransientTerm`,
   then you must specify the timestep by passing a ``dt=`` argument when calling
   :meth:`~fipy.terms.term.Term.solve` or :meth:`~fipy.terms.term.Term.sweep`.

The remaining changes are not *required*, but they make scripts easier to read
and we recommend them. :term:`FiPy` may issue a :exc:`DeprecationWarning` for some cases,
to indicate that we may not maintain the old syntax indefinitely.

 * "getter" and "setter" methods have been replaced with properties, e.g., use

   >>> x, y = mesh.cellCenters

   instead of

   >>> x, y = mesh.getCellCenters()

 * Boundary conditions are better applied with the
   :meth:`~fipy.variables.cellVariable.CellVariable.constrain` method than with
   the old :class:`~fipy.boundaryConditions.fixedValue.FixedValue` and
   :class:`~fipy.boundaryConditions.fixedFlux.FixedFlux` classes. See
   :ref:`BoundaryConditions`.

 * Individual :class:`~fipy.meshes.mesh.Mesh` classes should be imported
   directly from :mod:`fipy.meshes` and not :mod:`fipy.meshes.numMesh`.

 * The :term:`Gmsh` meshes now have simplified names:
   :class:`~fipy.meshes.gmshMesh.Gmsh2D` instead of
   :class:`~fipy.meshes.gmshMesh.GmshImporter2D`,
   :class:`~fipy.meshes.gmshMesh.Gmsh3D` instead of
   :class:`~fipy.meshes.gmshMesh.GmshImporter3D`, and
   :class:`~fipy.meshes.gmshMesh.Gmsh2DIn3DSpace` instead of
   :class:`~fipy.meshes.gmshMesh.GmshImporter2DIn3DSpace`.

.. _mailing list:         http://www.ctcms.nist.gov/fipy/mail.html
"""
__docformat__ = 'restructuredtext'
