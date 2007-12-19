#!/usr/bin/env python


from fipy.tools import numerix, vector

from fipy.meshes.numMesh.mesh import Mesh
from fipy.meshes.numMesh.mesh1D import Mesh1D
from fipy.meshes.numMesh.mesh2D import Mesh2D

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
            
        def move(self, monitorVariables=(), interpolants=None, beta=0.8, dt=None, sweeps=100):
            
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
                
                maxVertexRadii = 1. * (mesh._getFaceToCellDistances().sum(0) / 2)._getArithmeticVertexValue()
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

                TSVViewer(vars=interpolants + (self.getCellVolumes(),)).plot("interim%d.txt" % sweep)

    return _MovingMesh(mesh)
    

from fipy import *

from fipy.meshes.numMesh.grid2D import Grid2D as NonUniformGrid2D
from fipy.meshes.numMesh.grid1D import Grid1D as NonUniformGrid1D

dx = 5e-5 # cm
## nx = 40
## dx = 2e-5 # cm
nx = 40
L = nx * dx
ny = 40
dy = dx

mesh1 = Grid2D(nx=nx, ny=ny, dx=dx, dy=dy)
mesh2 = NonUniformGrid2D(nx=nx, ny=ny, dx=dx, dy=dy)
mesh3 = MovingMesh(mesh1)
mesh4 = NonUniformGrid1D(nx=nx, dx=dx)
mesh = MovingMesh(mesh2)


x, y = mesh.getCellCenters()
## x,  = mesh.getCellCenters()

phase = CellVariable(name="phase", mesh=mesh, hasOld=1)
phase.setValue(1.)
phase.setValue(0., where=(x > L/2) & (y > L/2))

Lv = 2350 # J / cm**3
Tm = 1728. # K
T = Variable(value=Tm)
enthalpy = Lv * (T - Tm) / Tm # J / cm**3

delta = 1.5 * dx
sigma = 3.7e-5 # J / cm**2
beta = 0.33 # cm / (K s)
kappa = 6 * sigma * delta # J / cm
W = 6 * sigma / delta # J / cm**3
Mphi = Tm * beta / (6. * Lv * delta) # cm**3 / (J s)



mPhi = -((1 - 2 * phase) * W + 30 * phase * (1 - phase) * enthalpy)
dmPhidPhi = 2 * W - 30 * (1 - 2 * phase) * enthalpy
S1 = dmPhidPhi * phase * (1 - phase) + mPhi * (1 - 2 * phase)
S0 = mPhi * phase * (1 - phase) - S1 * phase
eq = TransientTerm(coeff=1/Mphi) == ImplicitDiffusionTerm(coeff=kappa) \
                        + S0 + ImplicitSourceTerm(coeff = S1)


## rho = CellVariable(mesh=mesh, 
##                    value=1 + 10. * exp(-50 * (y - 1./2 - 1./4 * sin(2*pi*x))**2))

## t = 0.75
## rho = CellVariable(mesh=mesh, 
##                    name="rho",
##                    value=1 + 10. * exp(-50 * abs((x - 1./2 - 1./4 * cos(2*pi*t))**2
##                                                  + (y - 1./2 - 1./4 * sin(2*pi*t))**2
##                                                  - (1./10)**2)))

## TSVViewer(vars=(rho, rho.getGrad().dot(rho.getGrad()), mesh.getCellVolumes())).plot("test0.txt")
## TSVViewer(vars=(phase, mesh.getCellVolumes())).plot("test%d-%d.txt" % (0, 0))

print "0-0", phase.getCellVolumeAverage()

TSVViewer(vars=(phase, mesh.getCellVolumes())).plot("test%d-before.txt" % 0)

mesh.move(monitorVariables=(phase,), interpolants=(phase,), beta=0.99, sweeps=1)

TSVViewer(vars=(phase, mesh.getCellVolumes())).plot("test%d-after.txt" % 0)

print 0, phase.getCellVolumeAverage()

timeStep = 1e-6
for step in range(1,10):
    
    TSVViewer(vars=(phase, mesh.getCellVolumes())).plot("test%d-before.txt" % step)

    mesh.move(monitorVariables=(phase,), interpolants=(phase,), beta=0.99, sweeps=1)

    TSVViewer(vars=(phase, mesh.getCellVolumes())).plot("test%d-after.txt" % step)

    print step, phase.getCellVolumeAverage()
        
    phase.updateOld()
    res = 1e+10
    while res > 1e-9:
        res = eq.sweep(var=phase, dt=timeStep)
        print "solve residual:", res


T.setValue(T() - 1)

velocity = beta * abs(Tm - T()) # cm / s
timeStep = .1 * dx / velocity # s
elapsed = 0
while elapsed < 0.5 * L / velocity:
    step += 1
    
    TSVViewer(vars=(phase, mesh.getCellVolumes())).plot("test%d-before.txt" % step)

    mesh.move(monitorVariables=(phase,), interpolants=(phase,), beta=0.99, sweeps=1)

    TSVViewer(vars=(phase, mesh.getCellVolumes())).plot("test%d-after.txt" % step)

    print step, phase.getCellVolumeAverage()

    
    phase.updateOld()
    res = 1e+10
    while res > 1e-9:
        res = eq.sweep(var=phase, dt=timeStep)
        print "solve residual:", res
    elapsed += timeStep

