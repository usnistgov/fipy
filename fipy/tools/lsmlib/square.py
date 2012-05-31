from fipy import CellVariable, Viewer, Grid2D
from pylsmlib import computeDistanceFunction2d 

Nlsm = 999
Llsm = 2.

N = Nlsm + 1
L = Llsm + Llsm / Nlsm
dx = L / N
m = Grid2D(nx=N, ny=N, Lx=L, Ly=L) - ((dx / 2,), (dx / 2,))
print m.x
print dx
phi = CellVariable(mesh=m, value=1.e-2)
mask = (m.x < (3 * Llsm / 4)) & (m.x > (Llsm / 4)) & (m.y < (3 * Llsm / 4)) & (m.y > (Llsm / 4))
print mask
##phi[mask] = -1.e-2
phi.setValue(-1.e-2, where=mask)
##vi = Viewer(phi)

##vi.plot()
##raw_input('stopped')
phi[:] = computeDistanceFunction2d(phi.value, nx=m.nx, ny=m.ny, dx=m.dx, dy=m.dy)

##vi.plot()
##raw_input('stopped')

import numpy as np
import pylab
r = np.reshape(phi.value, (N, N))
print phi.value
pylab.contourf(r)
pylab.show()
raw_input('stopped')
