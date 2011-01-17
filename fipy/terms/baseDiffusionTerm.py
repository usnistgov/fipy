#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "baseDiffusionTerm.py"
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
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

import os

from fipy.terms.unaryTerm import _UnaryTerm
from fipy.tools import numerix

class _BaseDiffusionTerm(_UnaryTerm):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, coeff = (1.,), var=None):
        """
        Create a `DiffusionTerm`.

        :Parameters:
          - `coeff`: `Tuple` or `list` of `FaceVariables` or numbers.
          
        """
        if self.__class__ is _BaseDiffusionTerm:
            raise NotImplementedError, "can't instantiate abstract base class"
        
        if type(coeff) not in (type(()), type([])):
            coeff = [coeff]

        self.order = len(coeff) * 2


        if len(coeff) > 0:
            self.nthCoeff = coeff[0]

            from fipy.variables.variable import Variable
            if not isinstance(self.nthCoeff, Variable):
                self.nthCoeff = Variable(value = self.nthCoeff)

            from fipy.variables.cellVariable import CellVariable
            if isinstance(self.nthCoeff, CellVariable):
                self.nthCoeff = self.nthCoeff.getArithmeticFaceValue()

        else:
            self.nthCoeff = None

        _UnaryTerm.__init__(self, coeff=coeff, var=var)
        
        if self.order > 0:
            self.lowerOrderDiffusionTerm = self.__class__(coeff = coeff[1:])

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            self.coeff[0] = other * self.coeff[0] 
            return self.__class__(coeff=self.coeff, var=self.var)
        else:
            raise Exception, "Must multiply terms by int or float."

    __rmul__ = __mul__
    
    def __neg__(self):
        negatedCoeff = list(self.coeff)
        negatedCoeff[0] = -negatedCoeff[0]
        return self.__class__(coeff=negatedCoeff, var=self.var)
            
    def __getBoundaryConditions(self, boundaryConditions):
        higherOrderBCs = []
        lowerOrderBCs = []

        for bc in boundaryConditions:
            bcDeriv = bc._getDerivative(self.order - 2)
            if bcDeriv:
                higherOrderBCs.append(bcDeriv)
            else:
                lowerOrderBCs.append(bc)
                
        return higherOrderBCs, lowerOrderBCs

    def __getRotationTensor(self, mesh):
        if not hasattr(self, 'rotationTensor'):

            from fipy.variables.faceVariable import FaceVariable
            rotationTensor = FaceVariable(mesh=mesh, rank=2)
            
            rotationTensor[:, 0] = self._getNormals(mesh)

            if mesh.getDim() == 2:
                rotationTensor[:,1] = rotationTensor[:,0].dot((((0, 1), (-1, 0))))
            elif mesh.getDim() ==3:
                epsilon = 1e-20

                div = numerix.sqrt(1 - rotationTensor[2,0]**2)
                flag = numerix.resize(div > epsilon, (mesh.getDim(), mesh._getNumberOfFaces()))

                rotationTensor[0, 1] = 1
                rotationTensor[:, 1] = numerix.where(flag,
                                                     rotationTensor[:,0].dot((((0, 1, 0), (-1, 0, 0), (0, 0, 0)))) / div,
                                                     rotationTensor[:, 1])

                rotationTensor[1, 2] = 1
                rotationTensor[:, 2] = numerix.where(flag,
                                                     rotationTensor[:,0] * rotationTensor[2,0] / div,
                                                     rotationTensor[:, 2])
                rotationTensor[2, 2] = -div
                
            self.rotationTensor = rotationTensor

        return self.rotationTensor

    def __calcAnisotropySource(self, coeff, mesh, var):

        if not hasattr(self, 'anisotropySource'):
            if len(coeff) > 1:
                if hasattr(var.getArithmeticFaceValue(), 'constraints'):                
                    varNoConstraints = var.copy()
                else:
                    varNoConstraints = var
                gradients = varNoConstraints.getGrad().getHarmonicFaceValue().dot(self.__getRotationTensor(mesh))
                from fipy.variables.addOverFacesVariable import _AddOverFacesVariable
                self.anisotropySource = _AddOverFacesVariable(gradients[1:].dot(coeff[1:])) * mesh.getCellVolumes()

    def _calcGeomCoeff(self, mesh):
        if self.nthCoeff is not None:
          
            coeff = self.nthCoeff
            shape = numerix.getShape(coeff)

            from fipy.variables.faceVariable import FaceVariable
            if isinstance(coeff, FaceVariable):
                rank = coeff.getRank()
            else:
                rank = len(shape)

            if rank == 0 and self._treatMeshAsOrthogonal(mesh):
                tmpBop = (coeff * mesh._getFaceAreas() / mesh._getCellDistances())[numerix.newaxis, :]
            else:

                if rank == 1 or rank == 0:
                    coeff = coeff * numerix.identity(mesh.getDim())

                if rank > 0:
                    shape = numerix.getShape(coeff)
                    if mesh.getDim() != shape[0] or mesh.getDim() != shape[1]:
                        raise IndexError, 'diffusion coefficent tensor is not an appropriate shape for this mesh'          

                faceNormals = FaceVariable(mesh=mesh, rank=1, value=mesh._getFaceNormals())
                rotationTensor = self.__getRotationTensor(mesh)
                rotationTensor[:,0] = rotationTensor[:,0] / mesh._getCellDistances()
                
                tmpBop = faceNormals.dot(coeff).dot(rotationTensor) * mesh._getFaceAreas()

            return tmpBop

        else:

            return None

    def _getCoefficientMatrixForTests(self, SparseMatrix, mesh, coeff):
        """
        This method was introduced because __getCoefficientMatrix is private, but
        the tests in DiffusionTerm need to call it.
        """
        return self.__getCoefficientMatrix(SparseMatrix, mesh, coeff)

    def __getCoefficientMatrix(self, SparseMatrix, mesh, coeff):
        interiorCoeff = numerix.array(coeff)
        
        interiorCoeff[mesh.getExteriorFaces().getValue()] = 0
        
        interiorCoeff = numerix.take(interiorCoeff, mesh._getCellFaceIDs())

        coefficientMatrix = SparseMatrix(mesh=mesh, bandwidth = mesh._getMaxFacesPerCell() + 1)
        
        coefficientMatrix.addAtDiagonal(numerix.sum(interiorCoeff, 0))
        del interiorCoeff

        interiorFaces = mesh.getInteriorFaceIDs()
        interiorFaceCellIDs = mesh.getInteriorFaceCellIDs()

        interiorCoeff = -numerix.take(coeff, interiorFaces, axis=-1)
        coefficientMatrix.addAt(interiorCoeff, interiorFaceCellIDs[0], interiorFaceCellIDs[1])
        interiorCoeff = -numerix.take(coeff, interiorFaces, axis=-1)
        coefficientMatrix.addAt(interiorCoeff, interiorFaceCellIDs[1], interiorFaceCellIDs[0])
        
        return coefficientMatrix
        
    def __bcAdd(self, coefficientMatrix, boundaryB, LL, bb):
        coefficientMatrix += LL
        boundaryB += bb
        
    def __doBCs(self, SparseMatrix, higherOrderBCs, N, M, coeffs, coefficientMatrix, boundaryB):
        for boundaryCondition in higherOrderBCs:
            LL, bb = boundaryCondition._buildMatrix(SparseMatrix, N, M, coeffs)
            if os.environ.has_key('FIPY_DISPLAY_MATRIX'):
                self._viewer.title = r"%s %s" % (boundaryCondition.__class__.__name__, self.__class__.__name__)
                self._viewer.plot(matrix=LL, RHSvector=bb)
                from fipy import raw_input
                raw_input()
            self.__bcAdd(coefficientMatrix, boundaryB, LL, bb) 
            
        return coefficientMatrix, boundaryB

    def _buildMatrix(self, var, SparseMatrix, boundaryConditions=(), dt=1., transientGeomCoeff=None, diffusionGeomCoeff=None):

        if var is self.var or self.var is None:
        
            var, L, b = self.__higherOrderbuildMatrix(var, SparseMatrix, boundaryConditions=boundaryConditions, dt=dt, transientGeomCoeff=transientGeomCoeff, diffusionGeomCoeff=diffusionGeomCoeff)

            if self.order == 2:
                if not hasattr(self, 'constraintB'):

                    mesh = var.getMesh()
                    from fipy.variables.faceVariable import FaceVariable
    ##                normalsDotCoeff = FaceVariable(mesh=mesh, rank=1, value=mesh._getOrientedFaceNormals()) * self.nthCoeff
                    normalsDotCoeff = FaceVariable(mesh=mesh, rank=1, value=mesh._getOrientedFaceNormals()) * self.nthCoeff

                    self.constraintB = 0
                    self.constraintL = 0

                    if var.getFaceGrad().getConstraintMask() is not None:
                        self.constraintB -= (var.getFaceGrad().getConstraintMask() * self.nthCoeff * var.getFaceGrad()).getDivergence() * mesh.getCellVolumes()

                    if var.getArithmeticFaceValue().getConstraintMask() is not None:
                        constrainedNormalsDotCoeffOverdAP = var.getArithmeticFaceValue().getConstraintMask() * normalsDotCoeff / mesh._getCellDistances()
                        self.constraintB -= (constrainedNormalsDotCoeffOverdAP * var.getArithmeticFaceValue()).getDivergence() * mesh.getCellVolumes()
                        self.constraintL -= constrainedNormalsDotCoeffOverdAP.getDivergence() * mesh.getCellVolumes()

                L.addAtDiagonal(self.constraintL)
                b += self.constraintB

            return (var, L, b)
        else:
            return (var, SparseMatrix(mesh=var.getMesh()), 0)
        
    def __higherOrderbuildMatrix(self, var, SparseMatrix, boundaryConditions = (), dt = 1., transientGeomCoeff=None, diffusionGeomCoeff=None):
        mesh = var.getMesh()
        
        N = mesh.getNumberOfCells()
        M = mesh._getMaxFacesPerCell()

        if self.order > 2:

            higherOrderBCs, lowerOrderBCs = self.__getBoundaryConditions(boundaryConditions)
            
            var, lowerOrderL, lowerOrderb = self.lowerOrderDiffusionTerm._buildMatrix(var = var, SparseMatrix=SparseMatrix,
                                                                                      boundaryConditions = lowerOrderBCs, 
                                                                                      dt = dt, transientGeomCoeff=transientGeomCoeff,
                                                                                      diffusionGeomCoeff=diffusionGeomCoeff)
            del lowerOrderBCs
            
            lowerOrderb = lowerOrderb / mesh.getCellVolumes()
            volMatrix = SparseMatrix(mesh=var.getMesh(), bandwidth = 1)
            
            volMatrix.addAtDiagonal(1. / mesh.getCellVolumes() )
            lowerOrderL = volMatrix * lowerOrderL
            del volMatrix

            if not hasattr(self, 'coeffDict'):

                coeff = self._getGeomCoeff(mesh)[0]
                minusCoeff = -coeff
                
                coeff.dontCacheMe()
                minusCoeff.dontCacheMe()

                self.coeffDict = {
                    'cell 1 diag':     minusCoeff,
                    'cell 1 offdiag':  coeff
                    }
                del coeff
                del minusCoeff

                self.coeffDict['cell 2 offdiag'] = self.coeffDict['cell 1 offdiag']
                self.coeffDict['cell 2 diag'] = self.coeffDict['cell 1 diag']


            mm = self.__getCoefficientMatrix(SparseMatrix, mesh, self.coeffDict['cell 1 diag'])
            L, b = self.__doBCs(SparseMatrix, higherOrderBCs, N, M, self.coeffDict, 
                               mm, numerix.zeros(N,'d'))
                               
            del higherOrderBCs
            del mm

            b = L * lowerOrderb + b
            del lowerOrderb

            L = L * lowerOrderL
            del lowerOrderL

        elif self.order == 2:

            if not hasattr(self, 'coeffDict'):

                coeff = self._getGeomCoeff(mesh)
                minusCoeff = -coeff[0]

                coeff[0].dontCacheMe()
                minusCoeff.dontCacheMe()

                self.coeffDict = {
                    'cell 1 diag':    minusCoeff,
                    'cell 1 offdiag':  coeff[0]
                    }

                self.coeffDict['cell 2 offdiag'] = self.coeffDict['cell 1 offdiag']
                self.coeffDict['cell 2 diag'] = self.coeffDict['cell 1 diag']

                self.__calcAnisotropySource(coeff, mesh, var)

                del coeff
                del minusCoeff
                
            higherOrderBCs, lowerOrderBCs = self.__getBoundaryConditions(boundaryConditions)
            del lowerOrderBCs

            L, b = self.__doBCs(SparseMatrix, higherOrderBCs, N, M, self.coeffDict, 
                               self.__getCoefficientMatrix(SparseMatrix, mesh, self.coeffDict['cell 1 diag']), numerix.zeros(N,'d'))

            if hasattr(self, 'anisotropySource'):
                b -= self.anisotropySource
                               
            del higherOrderBCs


        else:
            
            L = SparseMatrix(mesh=mesh)
            L.addAtDiagonal(mesh.getCellVolumes())
            b = numerix.zeros((N),'d')
            
        return (var, L, b)

    def _getDiffusionGeomCoeff(self, mesh):
        return self._getGeomCoeff(mesh)

    def _getTransientGeomCoeff(self, mesh):
        return None  

    def _getDefaultSolver(self, solver, *args, **kwargs):
        return None

def _test(): 
    import doctest
    return doctest.testmod()

if __name__ == "__main__":
    _test()
