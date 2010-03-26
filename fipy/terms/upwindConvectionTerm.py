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


from fipy.terms.convectionTerm import ConvectionTerm
from fipy.variables.faceVariable import FaceVariable
from fipy.tools.dimensions.physicalField import PhysicalField
from fipy.tools import inline
from fipy.tools import numerix

class UpwindConvectionTerm(ConvectionTerm):
    r"""
    The discretization for this :class:`~fipy.terms.term.Term` is given by

    .. math::
    
       \int_V \nabla \cdot (\vec{u} \phi)\,dV \simeq \sum_{f} (\vec{n}
       \cdot \vec{u})_f \phi_f A_f

    where :math:`\phi_f=\alpha_f \phi_P +(1-\alpha_f)\phi_A` and
    :math:`\alpha_f` is calculated using the upwind convection scheme.
    For further details see :ref:`sec:NumericalSchemes`.
    """
    class _Alpha(FaceVariable):
        def __init__(self, P):
            FaceVariable.__init__(self, mesh = P.getMesh())
            self.P = self._requires(P)
            
        def _calcValuePy(self, P):
            alpha = numerix.where(P > 0., 1., 0.)
            return PhysicalField(value = alpha)

        def _calcValueIn(self, P):
            alpha = self._getArray().copy()
            inline._runInline("""
                alpha[i] = 0.5;
                
                if (P[i] > 0.) {
                    alpha[i] = 1.;
                } else {
                    alpha[i] = 0.;
                }
            """,
            alpha = alpha, P = P,
            ni = self.mesh._getNumberOfFaces()
            )

            return self._makeValue(value = alpha)

        def _calcValue(self):
            P  = self.P.getNumericValue()

            return inline._optionalInline(self._calcValueIn, self._calcValuePy, P)

    def _testPecletSign(self):
        r"""
            >>> from fipy import *
            >>> L = 10
            >>> mesh = Grid1D(nx=L, dx=1.)
            >>> var = CellVariable(mesh=mesh)

            >>> convCoeff = FaceVariable(mesh=mesh, rank=1)
            >>> diffCoeff = FaceVariable(mesh=mesh)
            >>> x = mesh.getFaceCenters()[0]
            >>> convCoeff.setValue((-1e+6,), where=x < 6)
            >>> diffCoeff.setValue(1e+3, where=x > 6)
            >>> print convCoeff
            [[-1000000. -1000000. -1000000. -1000000. -1000000. -1000000.        0.
                     0.        0.        0.        0.]]
            >>> print diffCoeff
            [    0.     0.     0.     0.     0.     0.     0.  1000.  1000.  1000.
              1000.]

            >>> BCs = (FixedValue(value=1., faces=mesh.getFacesRight()), 
            ...        FixedValue(value=0., faces=mesh.getFacesLeft()))
            >>> dTerm = DiffusionTerm(1000)

            >>> x = mesh.getCellCenters()[0]
            >>> boundary = 5.5
            >>> analytical = ((x - boundary) / (L - boundary)) * (x > boundary)

            >>> eqn = UpwindConvectionTerm(coeff=convCoeff) == dTerm
            >>> eqn.solve(var, boundaryConditions=BCs)
            >>> print var.allclose(analytical, atol=1e-3)
            1

            >>> eqn = TransientTerm(1e-10) == UpwindConvectionTerm(coeff=-convCoeff) +  dTerm
            >>> eqn.solve(var, dt = 1e+10, boundaryConditions=BCs)
            >>> print var.allclose(analytical, atol=1e-3)
            1

            >>> eqn = 0 == UpwindConvectionTerm(coeff=-convCoeff) +  dTerm
            >>> eqn.solve(var, boundaryConditions=BCs)
            >>> print var.allclose(analytical, atol=1e-3)
            1

            >>> eqn = 0 == -UpwindConvectionTerm(coeff=convCoeff) +  dTerm
            >>> eqn.solve(var, boundaryConditions=BCs)
            >>> print var.allclose(analytical, atol=1e-3)
            1
        """
        pass

def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
