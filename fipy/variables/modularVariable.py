#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "modularVariable.py"
 #                                    created: 12/8/03 {5:47:27 PM} 
 #                                last update: 12/28/05 {11:21:33 AM} 
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
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'
 
from fipy.variables.cellVariable import CellVariable

class ModularVariable(CellVariable):
    r"""
    The `ModularVariable` defines a variable that exisits on the circle between

    .. raw:: latex

        $-\pi$ and $\pi$

    The following examples show how `ModularVariable` works. When
    subtracting the answer wraps back around the circle.

        >>> from fipy.meshes.grid1D import Grid1D
        >>> mesh = Grid1D(nx = 2)
        >>> from fipy.tools import numerix
        >>> pi = numerix.pi
        >>> v1 = ModularVariable(mesh = mesh, value = (2*pi/3, -2*pi/3))
        >>> v2 = ModularVariable(mesh = mesh, value = -2*pi/3)
        >>> print numerix.allclose(v2 - v1, (2*pi/3, 0))
        1

    Obtaining the arithmetic face value.

        >>> print numerix.allclose(v1.getArithmeticFaceValue(), (2*pi/3, -pi, -2*pi/3))
        1

    Obtaining the gradient.

        >>> print numerix.allclose(v1.getGrad(), (pi/3, pi/3))
        1

    Obtaining the gradient at the faces.

        >>> print numerix.allclose(v1.getFaceGrad(), [[0], [2*pi/3], [0]])
        1
        
    Obtaining the gradient at the faces but without modular
    arithmetic.

        >>> print numerix.allclose(v1.getFaceGradNoMod(), [[0], [-4*pi/3], [0]])
        1
        
    """    
    def __init__(self, mesh, name = '', value=0., unit = None, hasOld = 0):
        """
        Creates a `ModularVariable` object.

        :Parameters:
          - `mesh`: The mesh that defines the geometry of this variable.
          - `name`: The user-readable name of the variable.
	  - `value`: The initial value.
	  - `unit`: The physical units of the variable.
	  - `hasOld`: Whether a variable keeps its old variable.
	
        """
	CellVariable.__init__(self, mesh = mesh, name = name, value = value, unit = unit, hasOld = hasOld)
        self.arithmeticFaceValue = None
        self.grad = None
        self.faceGradNoMod = None
        
    _modIn = """
    # define pi 3.141592653589793
    # define mod(x) (fmod(x + 3. * pi, 2. * pi) - pi)
    """

    def _setValue(self, value, unit = None, array = None):
        """
           >>> from fipy.meshes.grid1D import Grid1D
           >>> mesh = Grid1D(nx = 4)
           >>> from fipy.variables.modularVariable import ModularVariable
           >>> var = ModularVariable(mesh = mesh, value = 1, hasOld = 1)
           >>> print var
           [ 1., 1., 1., 1.,] 1
           >>> var.setValue(1)
           >>> print var
           [ 1., 1., 1., 1.,] 1

        """
        value = self._makeValue(value = value, unit = unit, array = array)
        from fipy.variables.modPhysicalField import _ModPhysicalField
	self.value = _ModPhysicalField(value = value, unit = unit, array = array)
	
    def updateOld(self):
        """
        Set the values of the previous solution sweep to the current values.
        Test case due to bug.

            >>> from fipy.meshes.grid1D import Grid1D
            >>> mesh = Grid1D(nx = 1)
            >>> var = ModularVariable(mesh = mesh, value = 1, hasOld = 1)
            >>> var.updateOld()
            >>> var[:] = 2
            >>> print var.getOld()
            [ 1.,] 1
            
        """
	self.setValue(self.getValue().mod(self()))
        if self.old is not None:
	    self.old.setValue(self.value.value.copy())

    def getGrad(self):
        r"""
        Return
        
        .. raw:: latex
        
           \( \nabla \phi \)
           
        as a `VectorCellVariable` (first-order gradient).
        Adjusted for a `ModularVariable`
        """
	if self.grad is None:

            from fipy.variables.modCellGradVariable import _ModCellGradVariable
            self.grad = _ModCellGradVariable(self, self._modIn, self.value.mod)

	return self.grad

    def getArithmeticFaceValue(self):
        r"""
        Returns a `FaceVariable` whose value corresponds to the arithmetic interpolation
        of the adjacent cells:
            
        .. raw:: latex
        
           \[ \phi_f = (\phi_1 - \phi_2) \frac{d_{f2}}{d_{12}} + \phi_2 \]

        Adjusted for a `ModularVariable`
           
        """
	if self.arithmeticFaceValue is None:
	    from modCellToFaceVariable import _ModCellToFaceVariable
	    self.arithmeticFaceValue = _ModCellToFaceVariable(self, self._modIn)

 	return self.arithmeticFaceValue

    def getFaceGrad(self):
        r"""
        Return
        
        .. raw:: latex
        
           \( \nabla \phi \)
           
        as a `VectorFaceVariable` (second-order gradient).
        Adjusted for a `ModularVariable`
        """
	if self.faceGrad is None:
	    from modFaceGradVariable import _ModFaceGradVariable
	    self.faceGrad = _ModFaceGradVariable(self, self._modIn)

	return self.faceGrad

    def getFaceGradNoMod(self):
        r"""
        
        .. raw:: latex
        
           \( \nabla \phi \)
           
        as a `VectorFaceVariable` (second-order gradient).
        Not adjusted for a `ModularVariable`
        """
        
        if self.faceGradNoMod is None:
            class NonModularTheta(CellVariable):
                def __init__(self, modVar):
                    CellVariable.__init__(self, mesh = modVar.getMesh())
                    self.modVar = self._requires(modVar)
                    
                def _calcValue(self):
                    return self.modVar[:]

            self.faceGradNoMod = NonModularTheta(self).getFaceGrad()

        return self.faceGradNoMod

    def __sub__(self, other):
        from fipy.terms.term import Term
        if isinstance(other, Term):
            return -other + self
        else:
            return self._getBinaryOperatorVariable(lambda a,b: a-b, other, canInline=False)
	
    def __rsub__(self, other):
	return self._getBinaryOperatorVariable(lambda a,b: b-a, other, canInline=False)


def _test(): 
    import doctest
    return doctest.testmod()
    
if __name__ == "__main__": 
    _test() 
