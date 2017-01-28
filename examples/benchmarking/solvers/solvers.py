import time
import uuid 

import datreant.core as dtr

import fipy as fp
from fipy.tools import numerix
from fipy.tools.parser import parse

output = parse('--output', action='store',
               type='string', default=str(uuid.uuid4()))

print "storing results in {0}".format(output)

data = dtr.Treant(output)

numberOfElements = parse('--numberOfElements', action='store',
                         type='int', default=10000)
data.categories['numberOfElements'] = numberOfElements

tolerance =  parse('--tolerance', action='store',
                   type='float', default=1e-10)
data.categories['tolerance'] = tolerance

iterations =  parse('--iterations', action='store',
                   type='int', default=1000)
data.categories['iterations'] = iterations

writeFiles = parse('--writeFiles', action='store_true', default=False)

data.categories['solver'] = fp.solvers.solver


N = int(numerix.sqrt(numberOfElements))
mesh = fp.Grid2D(nx=N, Lx=1., ny=N, Ly=1.)

var = fp.CellVariable(mesh=mesh, value=1., hasOld=True)
var.constrain(1., where=mesh.facesLeft)
var.constrain(0., where=mesh.facesRight)

eq = fp.TransientTerm() == fp.DiffusionTerm(coeff=var)

solver = fp.LinearPCGSolver(tolerance=tolerance, iterations=iterations, precon=None)
if fp.solvers.solver == "trilinos":
    # PySparse does b-normalization for (P)CG
    solver.convergenceCheck = solver.AZ_rhs

start = time.clock()

for sweep in range(10):
    eq.cacheMatrix()
    eq.cacheRHSvector()
    res = eq.sweep(var=var, dt=1., solver=solver)
    
    if writeFiles:
        eq.matrix.exportMmf(data["sweep{0}.mtx".format(sweep)].make().abspath)
        fp.tools.dump.write((var, eq.RHSvector), 
                            filename=data["sweep{0}.tar.gz".format(sweep)].make().abspath)

data.categories['elapsed'] = time.clock() - start
