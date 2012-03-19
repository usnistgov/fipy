#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "explicitDiffusionTerm.py"
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

__docformat__ = 'restructuredtext'

from fipy.terms.abstractDiffusionTerm import _AbstractDiffusionTerm

__all__ = ["ExplicitDiffusionTerm"]

class ExplicitDiffusionTerm(_AbstractDiffusionTerm):
    r"""
    The discretization for the `ExplicitDiffusionTerm` is given by

    .. math::

       \int_V \nabla \cdot (\Gamma\nabla\phi) dV \simeq \sum_f \Gamma_f
       \frac{\phi_A^\text{old}-\phi_P^\text{old}}{d_{AP}} A_f
       
    where :math:`\phi_A^\text{old}` and :math:`\phi_P^\text{old}` are the old values of the
    variable. The term is added to the RHS vector and makes no contribution to
    the solution matrix.

    """
    
    def _buildMatrix(self, var, SparseMatrix, boundaryConditions = (), dt=None, transientGeomCoeff=None, diffusionGeomCoeff=None):
        if hasattr(var, 'old'):
            varOld = var.old
        else:
            varOld = var
            
        varOld, L, b = _AbstractDiffusionTerm._buildMatrix(self, varOld, SparseMatrix, boundaryConditions = boundaryConditions, dt = dt,
                                                  transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)

        return (var, SparseMatrix(mesh=var.mesh), b - L * var.value)
        
    def _getNormals(self, mesh):
        return mesh._faceCellToCellNormals

    def _treatMeshAsOrthogonal(self, mesh):        
        return mesh._isOrthogonal

