#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "modPhysicalField.py"
 #                                    created: 12/28/03 {10:56:55 PM} 
 #                                last update: 2/2/04 {3:25:20 PM} 
 #  Author: Jonathan Guyer
 #  E-mail: guyer@nist.gov
 #  Author: Daniel Wheeler
 #  E-mail: daniel.wheeler@nist.gov
 #    mail: NIST
 #     www: http://ctcms.nist.gov
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  PFM is an experimental
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
 #  Copyright 1997-2004 by Konrad Hinsen, except as noted below.
 # 
 #  Permission to use, copy, modify, and distribute this software and its
 #  documentation for any purpose and without fee is hereby granted,
 #  provided that the above copyright notice appear in all copies and that
 #  both that copyright notice and this permission notice appear in
 #  supporting documentation.
 # 
 #  THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
 #  INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
 #  EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, INDIRECT OR
 #  CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
 #  USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 #  OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 #  PERFORMANCE OF THIS SOFTWARE.
 #  
 #  Description: 
 # 
 # Physical fields or quantities with units
 #
 # Based on PhysicalQuantities of the Scientific package, written by Konrad
 # Hinsen <hinsen@cnrs-orleans.fr>
 #
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-12-28 JEG 1.0 original
 #  1998/09/29 GPW     now supports conversions with offset 
 #                     (for temperature units)
 #  1998/09/28 GPW     now removes __args__ from local dict after eval
 # ###################################################################
 ##

from fipy.tools.dimensions.physicalField import PhysicalField

import Numeric

class ModPhysicalField(PhysicalField):

    def mod(self, argument):
        return Numeric.fmod(argument + 3. * Numeric.pi, 2. * Numeric.pi) - Numeric.pi

    def __sub__(self, other):
        if isinstance(other, ModPhysicalField):
            return self.__class__(value = self.mod(self.value - other.value), unit = self.unit)
        else:
            return self._sum(other, sign2 = lambda b: -b)
    
    def __rsub__(self, other):
        if isinstance(other, ModPhysicalField):
            return self.__class__(value = self.mod(argument = other.value - self.value), unit = self.unit)
        else:
            return self._sum(other, sign1 = lambda a: -a)
