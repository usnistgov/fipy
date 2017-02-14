import argparse
import time
import uuid 

import datreant.core as dtr

import fipy as fp
from fipy.tools import numerix
from fipy.tools import parallelComm

parser = argparse.ArgumentParser()
parser.add_argument("--output", help="directory to store results in",
                    default=str(uuid.uuid4()))
parser.add_argument("--numberOfElements", help="number of total cells in a Grid2D",
                    type=int, default=10000)
parser.add_argument("--solver", help="solver class to use",
                    choices=("cg", "pcg", "cgs", "gmres", "lu"), default="cg")
parser.add_argument("--sweeps", help="number of nonlinear sweeps to take",
                    type=int, default=10)
parser.add_argument("--iterations", help="maximum number of linear iterations to take for each sweep",
                    type=int, default=1000)
parser.add_argument("--tolerance", help="linear solver tolerance",
                    type=float, default=1e-10)
parser.add_argument("--writeFiles", help="whether to write solution values and matrix to OUTPUT",
                    action='store_true')

args, unknowns = parser.parse_known_args()

if parallelComm.procID == 0:
    print "storing results in {0}".format(args.output)
    data = dtr.Treant(args.output)
else:
    class dummyTreant(object):
        categories = dict()
        
    data = dummyTreant()

data.categories['processes'] = parallelComm.Nproc
data.categories['numberOfElements'] = args.numberOfElements
data.categories['sweeps'] = args.sweeps
data.categories['iterations'] = args.iterations
data.categories['tolerance'] = args.tolerance
data.categories['solver'] = args.solver
data.categories['library'] = fp.solvers.solver
data.categories['script'] = __file__

N = int(numerix.sqrt(args.numberOfElements))
mesh = fp.Grid2D(nx=N, Lx=1., ny=N, Ly=1.)

var = fp.CellVariable(mesh=mesh, value=1., hasOld=True)
var.constrain(1., where=mesh.facesLeft)
var.constrain(0., where=mesh.facesRight)

eq = fp.TransientTerm() == fp.DiffusionTerm(coeff=var)

if fp.solvers.solver != "scipy" and args.solver in ("pcg", "cgs", "gmres"):
    precon = fp.JacobiPreconditioner()
else:
    precon = None
    
if args.solver in ("cg", "pcg"):
    solver = fp.LinearPCGSolver(tolerance=args.tolerance, iterations=args.iterations, precon=precon)
    if fp.solvers.solver == "trilinos":
        # PySparse does b-normalization for (P)CG
        solver.convergenceCheck = solver.AZ_rhs
elif args.solver == "cgs":
    solver = fp.LinearCGSSolver(tolerance=args.tolerance, iterations=args.iterations, precon=precon)
elif args.solver == "gmres":
    solver = fp.LinearGMRESSolver(tolerance=args.tolerance, iterations=args.iterations, precon=precon)
elif args.solver == "lu":
    solver = fp.LinearLUSolver(tolerance=args.tolerance, iterations=args.iterations)
else:
    raise Exception("Unknown solver: {0}".format(args.solver))
        
start = time.clock()

for sweep in range(args.sweeps):
    eq.cacheMatrix()
    eq.cacheRHSvector()
    res = eq.sweep(var=var, dt=1., solver=solver)
    
    data.categories['sweep {0} - iterations'.format(sweep)] = solver.status['iterations']
    
    if args.writeFiles and parallelComm.procID == 0:
        eq.matrix.exportMmf(data["sweep{0}.mtx".format(sweep)].make().abspath)
        fp.tools.dump.write((var, eq.RHSvector), 
                            filename=data["sweep{0}.tar.gz".format(sweep)].make().abspath)

data.categories['elapsed'] = time.clock() - start
