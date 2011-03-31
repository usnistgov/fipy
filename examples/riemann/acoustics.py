from fipy import *

L = 2.
X0 = -1.
nx = 800
cfl = 10.0
K = 4.
rho = 1.

dx = L / nx
m = Grid1D(nx=nx, dx=dx) + X0
x, = m.getCellCenters()

q = CellVariable(mesh=m, rank=1)
q[0,:] = numerix.exp(-50 * (x - 0.3)**2) * numerix.cos(20 * (x - 0.3))
q[0, x > 0.3] = 0.

##p = CellVariable(mesh=m, value=numerix.exp(-50 * (x - 0.3)**2) * numerix.cos(20 * (x - 0.3)))
##p[x > 0.3] = 0
##u = CellVariable(mesh=m)

##ConvectionTerm = CentralDifferenceConvectionTerm


Ax = ((0, K), (1 / rho, 0))

eqn = TransientTerm() + ConvectionTerm((Ax,)) == 0

##eqn = (TransientTerm(var=p) + CentralDifferenceConvectionTerm(coeff=K, var=u)) & \
##      (TransientTerm(var=u) + CentralDifferenceConvectionTerm(coeff=1. / rho, var=p))

vi = Viewer(q)
vi.plot()
raw_input('stopped')
for step in range(100000):
    eqn.solve(q, dt=cfl * dx / 2.)
    print 'step',step
    if step % 10 ==  0:
        vi.plot()
##    raw_input('stopped')
