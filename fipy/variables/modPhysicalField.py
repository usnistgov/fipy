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

__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools.dimensions.physicalField import PhysicalField

from fipy.tools import numerix

class _ModPhysicalField(PhysicalField):

    def mod(argument):
        return numerix.fmod(argument + 3. * numerix.pi, 2. * numerix.pi) - numerix.pi
    mod = staticmethod(mod)

    def __sub__(self, other):
        if isinstance(other, _ModPhysicalField):
            return self.__class__(value= self.mod(self.inRadians() - other.inRadians()), unit="rad")
        else:
            return self._sum(other, sign2 = lambda b: -b)
    
    def __rsub__(self, other):
        if isinstance(other, _ModPhysicalField):
            return self.__class__(value = self.mod(argument=other.inRadians() - self.inRadians()), unit="rad")
        else:
            return self._sum(other, sign1 = lambda a: -a)

    def ravel(self):
        return self.value.ravel()
