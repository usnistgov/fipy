from fipy import CellVariable, Gmsh2D, DiffusionTerm, Viewer
from fipy.tools import numerix
from matplotlib import cm

monitor = None

geo = """\
c1_x = 6;
c1_y = 6;
c2_x = 4;
c2_y = 4;
r = 0.2;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {10, 0, 0, 1.0};
Point(3) = {10, 10, 0, 1.0};
Point(4) = {1, 10, 0, 1.0};
Point(5) = {0, 10, 0, 1.0};
Point(6) = {0, 9, 0, 1.0};
Point(7) = {c1_x, c1_y, 0, 1.0};
Point(8) = {c1_x, c1_y + r, 0, 1.0};
Point(9) = {c1_x + r, c1_y, 0, 1.0};
Point(10) = {c1_x, c1_y - r, 0, 1.0};
Point(11) = {c1_x - r, c1_y, 0, 1.0};
Point(12) = {c2_x, c2_y, 0, 1.0};
Point(13) = {c2_x, c2_y + r, 0, 1.0};
Point(14) = {c2_x + r, c2_y, 0, 1.0};
Point(15) = {c2_x, c2_y - r, 0, 1.0};
Point(16) = {c2_x - r, c2_y, 0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Circle(7) = {8, 7, 11};
Circle(8) = {11, 7, 10};
Circle(9) = {10, 7, 9};
Circle(10) = {9, 7, 8};
Circle(11) = {13, 12, 16};
Circle(12) = {16, 12, 15};
Circle(13) = {15, 12, 14};
Circle(14) = {14, 12, 13};
Line Loop(1) = {1, 2, 3, 4, 5, 6};
Line Loop(2) = {7, 8, 9, 10};
Line Loop(3) = {11, 12, 13, 14};
Plane Surface(1) = {1, 2, 3};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

Physical Line("Ground") = {4, 5, 6};
Physical Surface("Field") = {1};
Physical Surface("Anode") = {2};
Physical Surface("Cathode") = {3};
"""
for refinement in range(10):
    mesh = Gmsh2D(geo, background=monitor)

    charge = CellVariable(mesh=mesh, name=r"$\rho$", value=0.)
    charge.setValue(+1, where=mesh.physicalCells["Anode"])
    charge.setValue(-1, where=mesh.physicalCells["Cathode"])

    potential = CellVariable(mesh=mesh, name=r"$\psi$")
    potential.constrain(0., where=mesh.physicalFaces["Ground"])

    eq = DiffusionTerm(coeff=1.) == -charge

    res0 = eq.sweep(var=potential)

    res = eq.justResidualVector(var=potential)

    res1 = numerix.L2norm(res)
    res1a = CellVariable(mesh=mesh, value=abs(res))

    res = CellVariable(mesh=mesh, name="residual", value=abs(res) / mesh.cellVolumes**(1./mesh.dim) / 1e-3)

    # want cells no bigger than 1 and no smaller than 0.001
    maxSize = 1.
    minSize = 0.001
    monitor = CellVariable(mesh=mesh, name="monitor", value= 1. / (res + maxSize) +  minSize)

    viewer = Viewer(vars=potential, xmin=3.5, xmax=4.5, ymin=3.5, ymax=4.5)
#     viewer = Viewer(vars=(potential, charge))
    viewer.plot()

#     resviewer = Viewer(vars=res1a, log=True, datamin=1e-6, datamax=1e-2, cmap=cm.gray)
#     monviewer = Viewer(vars=monitor, log=True, datamin=1e-3, datamax=1)

    raw_input("refinement %d, res0: %g, res: %g:%g, N: %d, min: %g, max: %g, avg: %g. Press <return> to proceed..." \
              % (refinement, res0, res1, res1a.cellVolumeAverage, mesh.numberOfCells, numerix.sqrt(min(mesh.cellVolumes)), numerix.sqrt(max(mesh.cellVolumes)), numerix.mean(numerix.sqrt(mesh.cellVolumes))))
