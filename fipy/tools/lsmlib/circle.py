from fipy import CellVariable, Viewer, Grid2D
from pylsmlib import computeExtensionFields2d

CFL = 0.1
N = 100
L = 1.
m = Grid2D(nx=N, ny=N, Lx=L, Ly=L)
phi = CellVariable(mesh=m, name='distance function')
vel = CellVariable(mesh=m, rank=1, elementshape=(1,), name='velocity')
extVel = CellVariable(mesh=m, name=r'extension velocity')
phi[:] = (m.x - L / 2)**2 + (m.y - L / 2)**2 - 0.25**2
vel[:] = m.x**2 + m.y**2
extVel[:] = vel[:]
viVel = Viewer(vel[0])
viPhi = Viewer(phi, datamax=0.1, datamin=-0.1)
viExt = Viewer(extVel)
viPhi.plot()
viVel.plot()
viExt.plot()
raw_input('stopped')

phi[:], extVel[:] = computeExtensionFields2d(phi.value, vel.value, nx=m.nx, ny=m.ny, dx=m.dx, dy=m.dy)
viPhi.plot()
viVel.plot()
viExt.plot()
raw_input('stop')

for step in range(100000):
##    phi[:], vel[:] = computeDistanceFunction2d(phi.value, nx=m.nx, ny=m.ny, dx=m.dx, dy=m.dy)
##    print vel
##    print vel.value.shape
    tmpphi, extVel[:] = computeExtensionFields2d(phi.value, vel.value, nx=m.nx, ny=m.ny, dx=m.dx, dy=m.dy)
    dt = CFL * L / N / max(extVel[:])
    phi[:] = phi[:] - dt * extVel[:]

    if step % 10 == 0:
        print 'step',step
        viPhi.plot()
        viVel.plot()
        viExt.plot()
        raw_input('stopped')

