from fipy import *

L = 8.
X0 = -4.
nx = 3200
cfl = 10.0
K = 4.
rho = 1.

dx = L / nx
m = Grid1D(nx=nx, dx=dx) + X0
x, = m.getCellCenters()

p = CellVariable(mesh=m, value=numerix.exp(-50 * (x - 0.3)**2) * numerix.cos(20 * (x - 0.3)))
p[x > 0.3] = 0
u = CellVariable(mesh=m)

ConvectionTerm = CentralDifferenceConvectionTerm

TransientTerm(var=p)

eqn = (TransientTerm(var=p) + CentralDifferenceConvectionTerm(coeff=K, var=u)) & \
      (TransientTerm(var=u) + CentralDifferenceConvectionTerm(coeff=1. / rho, var=p))

vi = Viewer((p, u))
vi.plot()
raw_input('stopped')
for step in range(100000):
    eqn.solve(dt=cfl * dx / 2.)
    print 'step',step
    if step % 10 ==  0:
        vi.plot()
##    raw_input('stopped')
