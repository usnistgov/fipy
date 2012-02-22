from fipy import CellVariable, Viewer, Grid2D, numerix
from pylsmlib import solveEikonalEquation2d, computeDistanceFunction2d

CFL = 10.0
N = 1000
L = 1.
m = Grid2D(nx=N, ny=N, Lx=L, Ly=L)
phi = CellVariable(mesh=m, name='distance function')
slow = .5
vel = CellVariable(mesh=m, name='velocity', value=slow)
phi[:] = numerix.sqrt((m.x - L / 2)**2 + (m.y - L / 2)**2) - 0.25
vel.setValue(1., where=m.x > (L / 8))
vel.setValue(slow, where=m.x > (2 * L / 8))
vel.setValue(1., where=m.x > (3 * L / 8))
vel.setValue(slow, where=m.x > (4 * L / 8))
vel.setValue(1., where=m.x > (5 * L / 8))
vel.setValue(slow, where=m.x > (6 * L / 8))
vel.setValue(1., where=m.x > (7 * L / 8))



viVel = Viewer(vel)
viPhi = Viewer(phi, datamax=1e-10, datamin=-1e-10)

viPhi.plot()
viVel.plot()


viPhi.plot()
viVel.plot()

raw_input('stopped')
phi[:] = solveEikonalEquation2d(phi.value, vel.value, nx=m.nx, ny=m.ny, dx=m.dx, dy=m.dy)

viPhi.plot()
raw_input('stopped')
dt = CFL * m.dx / max(vel)

for step in range(100000):

    phi[:] = phi[:] - dt

    print 'step',step
    viPhi.plot()
    raw_input('stopped')

