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

__all__ = []

from fipy.meshes.builders.grid1DBuilder import _NonuniformGrid1DBuilder
 
class _PeriodicGrid1DBuilder(_NonuniformGrid1DBuilder):

    def buildGridData(self, *args, **kwargs):
        kwargs["cacheOccupiedNodes"] = True
        return super(_PeriodicGrid1DBuilder, self).buildGridData(*args, 
                                                                **kwargs)
          
    def _buildOverlap(self, overlap, procID, occupiedNodes):
        if occupiedNodes == 1:
            return super(_PeriodicGrid1DBuilder, self)._buildOverlap(overlap, 
                     procID, occupiedNodes)
        else:
            return (overlap, overlap, {'left': overlap, 'right': overlap})
            

