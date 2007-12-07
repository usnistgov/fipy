#!/usr/bin/env python

## -*-Pyth-*-
 # #############################################################################
 # FiPy - a finite volume PDE solver in Python
 # 
 # FILE: "movingMesh.py"
 #                                     created: 11/25/07 {6:52:31 PM}
 #                                 last update: 12/6/07 {5:55:29 PM}
 # Author: Jonathan Guyer
 # E-mail: <jguyer@his.com>
 #   mail: Alpha Cabal
 #    www: <http://alphatcl.sourceforge.net>
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
 # History
 # 
 # modified   by  rev reason
 # ---------- --- --- -----------
 # 1941-11-23 JEG 1.0 original
 # 
 # #############################################################################
 ##

from fipy.tools import numerix, vector

from mesh import Mesh
from mesh1D import Mesh1D
from mesh2D import Mesh2D

from fipy.variables.cellVariable import CellVariable

def MovingMesh(mesh):
    if isinstance(mesh, Mesh1D):
        base = Mesh1D
    elif isinstance(mesh, Mesh2D):
        base = Mesh2D
    else:
        base = Mesh
        
    class _MovingMesh(base):
        def __init__(self, mesh):
            self.lagrangianMesh = mesh
            base.__init__(self, 
                          vertexCoords=mesh.getVertexCoords(), 
                          faceVertexIDs=mesh._getFaceVertexIDs(), 
                          cellFaceIDs=mesh._getCellFaceIDs())
                          
            self.xi = CellVariable(mesh=mesh, value=mesh.getCellCenters(), rank=1, hasOld=True)
            
            from fipy.boundaryConditions.fixedValue import FixedValue
            self.bcs = []
            for i in range(self.getDim()):
                if i is 0:
                    self.bcs.append((FixedValue(faces=mesh.getFacesLeft(),
                                                value=mesh.getFacesLeft().getCenters()[0]),
                                     FixedValue(faces=mesh.getFacesRight(),
                                                value=mesh.getFacesRight().getCenters()[0])))
                elif i is 1:
                    self.bcs.append((FixedValue(faces=mesh.getFacesTop(),
                                                value=mesh.getFacesTop().getCenters()[1]),
                                     FixedValue(faces=mesh.getFacesBottom(),
                                                value=mesh.getFacesBottom().getCenters()[1])))
                elif i is 2:
                    self.bcs.append((FixedValue(faces=mesh.getFacesFront(),
                                                value=mesh.getFacesFront().getCenters()[2]),
                                     FixedValue(faces=mesh.getFacesBack(),
                                                value=mesh.getFacesBack().getCenters()[2])))
            
        def move(self, monitorVariables=(), interpolants=None, beta=0.8, dt=None, sweeps=10):
            
            if interpolants is None:
                interpolants = monitorVariables
                
            lmesh = self.lagrangianMesh

            for sweep in range(sweeps):
                self.xi.updateOld()

                monitor = 0.
                lagrangianVar = CellVariable(mesh=lmesh, value=0.)
                magGradSqrt = lagrangianVar.getGrad().getMag().sqrt()
##                 magGrad = lagrangianVar.getGrad().getMag()
                for var in monitorVariables:
                    lagrangianVar.setValue(var.getValue())
##                     monitor += ((1.-beta) * (magGradSqrt.getCellVolumeAverage() 
##                                              * lmesh.getNumberOfCells())
##                                 + beta * magGradSqrt)
##                     monitor += ((1.-beta) * (magGrad.getCellVolumeAverage() 
##                                              * lmesh.getNumberOfCells())
##                                 + beta * magGrad)
##                     monitor += ((1.-beta) * magGrad.getCellVolumeAverage() 
##                                 + beta * magGrad)
                    monitor += ((1.-beta) * magGradSqrt.getCellVolumeAverage() 
                                + beta * magGradSqrt)

                
                from fipy.terms.transientTerm import TransientTerm
                from fipy.terms.diffusionTerm import DiffusionTerm
                if dt is None:
                    eq = DiffusionTerm(coeff=monitor) 
                else:
                    eq = TransientTerm(coeff=dt) == DiffusionTerm(coeff=monitor)
                
                for i in range(self.getDim()):
                    val = self.xi[i].copy()
                    eq.solve(var=val, boundaryConditions=self.bcs[i])
                    self.xi[i] = val

                oldVolumes = self.getCellVolumes().getValue().copy()
                oldVertexCoords = self.vertexCoords.copy()
                
                newVertexCoords = self.xi.getArithmeticFaceValue()._getArithmeticVertexValue().copy()

                M, F = self._getFaceVertexIDs().shape
                for i in range(self.getDim()):
                    for bc in self.bcs[i]:
                        for j in range(M):
                            ids = numerix.take(self._getFaceVertexIDs()[j], bc.faces)
                            newVertexCoords[i,ids] = bc.value

                displacement = newVertexCoords - oldVertexCoords # self.vertexCoords
                displacementMag = numerix.sqrtDot(displacement, displacement)
                
                maxVertexRadii = 0.1 * (mesh._getFaceToCellDistances().sum(0) / 2)._getArithmeticVertexValue()
                factor = displacementMag / maxVertexRadii

                factor = max(factor)

                if factor > 1.:
                    displacement /= factor
                    displacementMag = numerix.sqrtDot(displacement, displacement)
                    newVertexCoords = self.vertexCoords + displacement

                displacmentLInfNorm = max(displacementMag)

                print sweep, "max displacement:", displacmentLInfNorm
                
                if displacmentLInfNorm < 1e-8:
                    break
                    
                self.vertexCoords[:] = newVertexCoords
                
                from fipy.terms.implicitSourceTerm import ImplicitSourceTerm
                from fipy.terms.upwindConvectionTerm import UpwindConvectionTerm
                from fipy.solvers import LinearLUSolver
                
                for var in interpolants: #self.getSubscribedVariables():
##                     var = var()
                    if var._isSolvable():
                        oldVar = var.copy()
                        updateEq = (ImplicitSourceTerm(coeff=-1.) + (oldVolumes / self.getCellVolumes()) * oldVar
                                    + UpwindConvectionTerm(coeff=displacement._getArithmeticFaceValue()) == 0)

                        res = 1e100
                        while res > 1e-12:
                            res = updateEq.sweep(var=var, solver=LinearLUSolver(tolerance=1.e-15, iterations=2000))
    ##                         res = updateEq.sweep(var=var, solver=LinearCGSSolver())
                            print "%d: %s move residual: %g" % (sweep, var.name, res)

    return _MovingMesh(mesh)
    
