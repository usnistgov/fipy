#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "modPhysicalField.py"
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
 #  Description: 
 # 
 # Physical fields or quantities with units
 #
 # ###################################################################
 ##

from fipy.tools.dimensions.physicalField import PhysicalField

from fipy.tools.numerix import pi, fmod

class _ModPhysicalField(PhysicalField):

    def mod(self, argument):
        return fmod(argument + 3. * pi, 2. * pi) - pi

    def __sub__(self, other):
        if isinstance(other, _ModPhysicalField):
            return self.__class__(value = self.mod(self.value - other.value), unit = self.unit)
        else:
            return PhysicalField.__sub__(self, other)
    
    def __rsub__(self, other):
        if isinstance(other, _ModPhysicalField):
            return self.__class__(value = self.mod(argument = other.value - self.value), unit = self.unit)
        else:
            return PhysicalField.__rsub__(self, other)

