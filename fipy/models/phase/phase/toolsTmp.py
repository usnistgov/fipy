#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  PyFiVol - Python-based finite volume PDE solver
 # 
 #  FILE: "tools.py"
 #                                    created: 11/12/03 {10:39:23 AM} 
 #                                last update: 1/16/04 {11:57:58 AM}
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
 # 
 #  Description: 
 # 
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-11-12 JEG 1.0 original
 # ###################################################################
 ##

import Numeric

import fivol.tools.array as array
import fivol.inline.inline as inline

def _addOverFaces(faceGradient, faceVariable, mesh, NCells):
##     contributions = Numeric.sum(mesh.getAreaProjections() * faceGradient,1)
    contributions = array.sum(mesh.getAreaProjections() * faceGradient,1)
    
    contributions = contributions * faceVariable
    
##     NIntFac = len(mesh.getInteriorFaces())
##     NExtFac = len(mesh.getFaces()) - NIntFac
    
##     contributions = Numeric.concatenate((contributions[:NIntFac], Numeric.zeros(NExtFac,'d')))
    contributions[len(mesh.getInteriorFaces()):] = 0
    ids = mesh.getCellFaceIDs()
    
##     contributions = Numeric.take(contributions, ids)
    contributions = array.take(contributions, ids)
    
##     NMaxFac = mesh.getMaxFacesPerCell()
    
##     contributions = Numeric.reshape(contributions,(NCells,-1))

    contributions = array.reshape(contributions,(NCells,-1))
    
##     orientations = Numeric.reshape(mesh.getCellFaceOrientations(),(NCells,-1))
    orientations = array.reshape(mesh.getCellFaceOrientations(),(NCells,-1))
    
##     return Numeric.sum(orientations*contributions,1) / mesh.getCellVolumes()
    return array.sum(orientations*contributions,1) / mesh.getCellVolumes()

def _addOverFacesInline(faceGradient, faceVariable, mesh, NCells):

    returnValue = Numeric.zeros((NCells),'d')

    inline.runInline("""
    int i;

    for(i = 0; i < numberOfInteriorFaces; i++)
      {
        int j;
        for(j = 0; j < numberOfDimensions; j++)
          {
            contributions(i) += areaProjections(i,j) * faceGradient(i,j);
          }
        contributions(i) = contributions(i) * faceVariable(i);
      }

    for(i = 0; i < numberOfCells; i++)
      {
        int j;
        for(j = 0; j < numberOfCellFaces; j++)
          {
            returnValue(i) += orientations(i,j) * contributions(ids(i,j));
          }
        returnValue(i) = returnValue(i) / cellVolume(i);
      }
    """,numberOfInteriorFaces = len(mesh.getInteriorFaces()),
                     numberOfDimensions = mesh.getDim(),
                     numberOfCellFaces = mesh.getMaxFacesPerCell(),
                     numberOfCells = NCells,
                     contributions =  Numeric.zeros((len(mesh.getFaces())),'d'),
                     areaProjections = mesh.getAreaProjections().value[:],
                     faceGradient = faceGradient.getNumericValue()[:],
                     faceVariable = faceVariable[:],
                     ids = Numeric.array(mesh.getCellFaceIDs()[:]),
                     returnValue = returnValue,
                     orientations = mesh.getCellFaceOrientations()[:],
                     cellVolume = mesh.getCellVolumes()[:])

    return returnValue

def addOverFaces(faceGradient = None, faceVariable = None, mesh = None, NCells = None):

    return _addOverFaces(faceGradient, faceVariable, mesh, NCells)
##    return inline.optionalInline(_addOverFacesInline, _addOverFaces, faceGradient, faceVariable, mesh)



    


