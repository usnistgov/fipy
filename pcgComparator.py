import os
import time

import fipy as fp
from fipy.tools import numerix
from fipy.tools.parser import parse

numberOfElements = parse('--numberOfElements', action='store',
                         type='int', default=10000)
N = int(numerix.sqrt(numberOfElements))

tolerance =  parse('--tolerance', action='store',
                   type='float', default=1e-10)

iterations =  parse('--iterations', action='store',
                   type='int', default=1000)

mesh = fp.Grid2D(nx=N, Lx=1., ny=N, Ly=1.)

var = fp.CellVariable(mesh=mesh, value=1., hasOld=True)
var.constrain(1., where=mesh.facesLeft)
var.constrain(0., where=mesh.facesRight)

eq = fp.TransientTerm() == fp.DiffusionTerm(coeff=var)

solver = fp.LinearPCGSolver(tolerance=tolerance, iterations=iterations, precon=None)

try:
    os.mkdir(fp.solvers.solver)
except:
    pass

start = time.clock()

for sweep in range(10):
    eq.cacheMatrix()
    eq.cacheRHSvector()
    res = eq.sweep(var=var, dt=1., solver=solver)
    eq.matrix.exportMmf(os.path.join(fp.solvers.solver, "sweep{0}.mmf".format(sweep)))
    fp.tools.dump.write((var, eq.RHSvector), 
                        filename=os.path.join(fp.solvers.solver,
                                              "sweep{0}.tar.gz".format(sweep)))

with open(os.path.join(fp.solvers.solver, "elapsed.txt"), 'w') as f:
    f.write(str(time.clock() - start))