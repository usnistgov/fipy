#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "upwindConvectionTerm.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.terms.abstractUpwindConvectionTerm import _AbstractUpwindConvectionTerm

__all__ = ["UpwindConvectionTerm"]

class UpwindConvectionTerm(_AbstractUpwindConvectionTerm):
    r"""
    The discretization for this :class:`~fipy.terms.term.Term` is given by

    .. math::
    
       \int_V \nabla \cdot (\vec{u} \phi)\,dV \simeq \sum_{f} (\vec{n}
       \cdot \vec{u})_f \phi_f A_f

    where :math:`\phi_f=\alpha_f \phi_P +(1-\alpha_f)\phi_A` and
    :math:`\alpha_f` is calculated using the upwind convection scheme.
    For further details see :ref:`sec:NumericalSchemes`.
    """

    def _getDefaultSolver(self, var, solver, *args, **kwargs):
        solver = solver or super(UpwindConvectionTerm, self)._getDefaultSolver(var, solver, *args, **kwargs)
        if solver and not solver._canSolveAsymmetric():
            import warnings
            warnings.warn("%s cannot solve assymetric matrices" % solver)
        from fipy.solvers import DefaultAsymmetricSolver
        return solver or DefaultAsymmetricSolver(*args, **kwargs)
    
    def _testPecletSign(self):
        r"""
            >>> from fipy import *
            >>> nx = 11
            >>> L = 1.
            >>> mesh = Grid1D(nx=nx, dx=L / nx)
            >>> var = CellVariable(mesh=mesh)

            >>> from fipy.variables.faceVariable import FaceVariable
            >>> convCoeff = FaceVariable(mesh=mesh, rank=1)
            >>> diffCoeff = FaceVariable(mesh=mesh, value=1e-20)
            >>> x = mesh.faceCenters[0]
            >>> D = 1000.
            >>> u = -1e+6
            >>> convCoeff.setValue((u,), where=x < L / 2)
            >>> diffCoeff.setValue(D, where=x > L / 2)

            >>> var.constrain(1., mesh.facesRight)
            >>> var.constrain(10, mesh.facesLeft)

            >>> dTerm = DiffusionTerm(diffCoeff)

            >>> x = mesh.cellCenters[0]
            >>> var0 = 2 * D / (-u * L + 2 * D)
            >>> analytical = 2 * (1 - var0) * x / L + 2 * var0 - 1
            >>> analytical = var0 * (x < L / 2) + analytical * (x >= L / 2)

            >>> var[:] = 0
            >>> eqn = UpwindConvectionTerm(coeff=convCoeff) == dTerm
            >>> eqn.solve(var)
            >>> print var.allclose(analytical)
            1

            >>> var[:] = 0
            >>> eqn = TransientTerm(1e-10) == UpwindConvectionTerm(coeff=-convCoeff) +  dTerm
            >>> eqn.solve(var, dt = 1e+10)
            >>> print var.allclose(analytical)
            1

            >>> var[:] = 0
            >>> eqn = 0 == UpwindConvectionTerm(coeff=-convCoeff) +  dTerm
            >>> eqn.solve(var)
            >>> print var.allclose(analytical)
            1

            >>> var[:] = 0
            >>> eqn = 0 == -UpwindConvectionTerm(coeff=convCoeff) +  dTerm
            >>> eqn.solve(var)
            >>> print var.allclose(analytical)
            1
            
        """
        pass

def _test(): 
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()
    
if __name__ == "__main__": 
    _test() 
