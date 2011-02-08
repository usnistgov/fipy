#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
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
 
from abstractGridBuilder import AbstractGridBuilder

from fipy.meshes.builders.utilityClasses import (UniformNumPts,
                                                 NonuniformNumPts)

class Grid2DBuilder(AbstractGridBuilder):

    def _packOverlap(self, first, second):
        return {'left': 0, 'right': 0, 'bottom': first, 'top': second}  

    def _packOffset(self, arg):
        return (0, arg)

class NonuniformGrid2DBuilder(Grid2DBuilder):

    def __init__(self):
        self.NumPtsCalcClass = NonuniformNumPts

        super(NonuniformGrid2DBuilder, self).__init__()

class UniformGrid2DBuilder(Grid2DBuilder):

    def __init__(self):
        self.NumPtsCalcClass = UniformNumPts

        super(UniformGrid2DBuilder, self).__init__()
                                  
