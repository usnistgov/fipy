from fipy import *

N = 80
L = 2.
dx = L / N 
origin =[[-1], [-1]]

mesh = Grid2D(nx=N, ny=N, dx=dx, dy=dx) + origin
x, y = mesh.cellCenters

var = CellVariable(mesh=mesh, hasOld=True)

var[(0.1 < x) & (x < 0.6) & (-0.25 < y) & (y < 0.25)] = 1.
r = numerix.sqrt((x + 0.45)**2 + y**2)
var.setValue(1 - r / 0.35, where=r < 0.35)
##var[r < 0.35] = 1 - r / 0.35

vel = CellVariable(mesh=mesh, rank=1)

vel[0] = 2 * y
vel[1] = -2 * x


##print 'dt',dt
eqn = TransientTerm() + RoeConvectionTerm(vel)
##print  roeConvectionTerm.maxeigenvalue(var)

viewer = Viewer(var)

##def run():
for i in range(100000):
    print 'i',i
    var.updateOld()
    eqn.solve(var, dt=0.0009)
    if i % 10 == 0:
        viewer.plot()

##import cProfile
##cProfile.runctx("run()", globals(), locals(), "Profile.prof")



##raw_input('stopped')
