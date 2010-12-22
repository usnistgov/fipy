#!/usr/bin/env python

## 
 # -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "mesh.py"
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
 
def _notImplemented(self):
    raise NotImplementedError

class AbstractScaledMeshGeometry(object):

    scale                     = property(_notImplemented)
    scaledFaceAreas           = property(_notImplemented)
    scaledCellVolumes         = property(_notImplemented)
    scaledCellCenters         = property(_notImplemented)
    scaledFaceToCellDistances = property(_notImplemented)
    scaledCellDistances       = property(_notImplemented)
    scaledCellToCellDistances = property(_notImplemented)
    areaProjections           = property(_notImplemented)
    orientedAreaProjections   = property(_notImplemented)
    faceToCellDistanceRatio   = property(_notImplemented)
    faceAspectRatios          = property(_notImplemented)

class AbstractMeshGeometry(object):

    faceAreas                 = property(_notImplemented)           
    faceCenters               = property(_notImplemented)           
    cellCenters               = property(_notImplemented)         
    faceToCellDistances       = property(_notImplemented)
    cellToFaceDistanceVectors = property(_notImplemented)
    cellDistances             = property(_notImplemented)
    cellDistanceVectors       = property(_notImplemented)  
    faceNormals               = property(_notImplemented)          
    orientedFaceNormals       = property(_notImplemented)  
    cellVolumes               = property(_notImplemented)          
    cellCenters               = property(_notImplemented)           
    faceCellToCellNormals     = property(_notImplemented)
    faceTangents1             = property(_notImplemented) 
    faceTangents2             = property(_notImplemented)        
    cellToCellDistances       = property(_notImplemented)  
    cellAreas                 = property(_notImplemented)  
    cellNormals               = property(_notImplemented)  

    """Scaled business (mostly wrappers)"""
    scale                     = 1.
    scaledFaceAreas           = property(_notImplemented)
    scaledCellVolumes         = property(_notImplemented)
    scaledCellCenters         = property(_notImplemented)
    scaledFaceToCellDistances = property(_notImplemented)
    scaledCellDistances       = property(_notImplemented)
    scaledCellToCellDistances = property(_notImplemented)
    areaProjections           = property(_notImplemented)
    orientedAreaProjections   = property(_notImplemented)
    faceToCellDistanceRatio   = property(_notImplemented)
    faceAspectRatios          = property(_notImplemented)
     
