"""

The `gapFillMesh` function glues 3 meshes together to form a composite
mesh. The first mesh is a `Grid2D` object that is fine and deals with
the area around the trench or via. The second mesh is a `Gmsh2D`
object that forms a transition mesh from a fine to a course
region. The third mesh is another `Grid2D` object that forms the
boundary layer. This region consists of very large elements and is
only used for the diffusion in the boundary layer.

"""
from __future__ import division
from __future__ import unicode_literals

__docformat__ = 'restructuredtext'

from fipy.meshes import Gmsh2D
from fipy.meshes.gmshMesh import _gmshVersion
from fipy.meshes import Grid2D
from fipy.tools import serialComm
from fipy.tools import parallelComm
from fipy.tools.version import Version

class GapFillMesh(Gmsh2D):
    """
    The following test case tests for diffusion across the domain.
    >>> domainHeight = 5.
    >>> mesh = GapFillMesh(transitionRegionHeight = 2.,
    ...                    cellSize = 0.1,
    ...                    desiredFineRegionHeight = 1.,
    ...                    desiredDomainHeight = domainHeight,
    ...                    desiredDomainWidth = 1.) # doctest: +GMSH

    >>> import fipy.tools.dump as dump
    >>> (f, filename) = dump.write(mesh) # doctest: +GMSH
    >>> if parallelComm.Nproc == 1:
    ...     mesh = dump.read(filename, f) # doctest: +GMSH

    >>> print(136 < mesh.globalNumberOfCells < 300) # doctest: +GMSH
    True

    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(mesh = mesh) # doctest: +GMSH

    >>> from fipy.terms.diffusionTerm import DiffusionTerm
    >>> eq = DiffusionTerm()

    >>> var.constrain(0., mesh.facesBottom) # doctest: +GMSH
    >>> var.constrain(domainHeight, mesh.facesTop) # doctest: +GMSH

    >>> eq.solve(var) # doctest: +GMSH

    Evaluate the result:

    >>> centers = mesh.cellCenters[1].copy() # doctest: +GMSH

    .. note:: the copy makes the array contiguous for inlining

    >>> localErrors = (centers - var)**2 / centers**2 # doctest: +GMSH
    >>> from fipy.tools import numerix
    >>> globalError = numerix.sqrt(numerix.sum(localErrors) / mesh.numberOfCells) # doctest: +GMSH
    >>> argmax = numerix.argmax(localErrors) # doctest: +GMSH

    >>> print(numerix.sqrt(localErrors[argmax]) < 0.1) # doctest: +GMSH
    1
    >>> print(globalError < 0.05) # doctest: +GMSH
    1

    """
    def __init__(self,
                 cellSize=None,
                 desiredDomainWidth=None,
                 desiredDomainHeight=None,
                 desiredFineRegionHeight=None,
                 transitionRegionHeight=None,
                 communicator=parallelComm):

        """
        Arguments:

        `cellSize` - The cell size in the fine grid around the trench.

        `desiredDomainWidth` - The desired domain width.

        `desiredDomainHeight` - The total desired height of the
        domain.

        `desiredFineRegionHeight` - The desired height of the in the
        fine region around the trench.

        `transitionRegionHeight` - The height of the transition region.
        """

        # Calculate the fine region cell counts.
        nx = int(desiredDomainWidth / cellSize)
        ny = int(desiredFineRegionHeight / cellSize)

        # Calculate the actual mesh dimensions
        actualFineRegionHeight = ny * cellSize
        actualDomainWidth = nx * cellSize
        boundaryLayerHeight = desiredDomainHeight - actualFineRegionHeight - transitionRegionHeight
        numberOfBoundaryLayerCells = int(boundaryLayerHeight / actualDomainWidth)

        # Build the fine region mesh.
        self.fineMesh = Grid2D(nx=nx, ny=ny, dx=cellSize, dy=cellSize, communicator=serialComm)

        if _gmshVersion() < Version("2.7"):
            # kludge: must offset cellSize by `eps` to work properly
            eps = float(cellSize)/(nx * 10)
        else:
            eps = 0.

        super(GapFillMesh, self).__init__("""
        ny       = %(ny)g;
        cellSize = %(cellSize)g - %(eps)g;
        height   = %(actualFineRegionHeight)g;
        width    = %(actualDomainWidth)g;
        boundaryLayerHeight = %(boundaryLayerHeight)g;
        transitionRegionHeight = %(transitionRegionHeight)g;
        numberOfBoundaryLayerCells = %(numberOfBoundaryLayerCells)g;

        Point(1) = {0, 0, 0, cellSize};
        Point(2) = {width, 0, 0, cellSize};
        Line(3) = {1, 2};

        Point(10) = {0, height, 0, cellSize};
        Point(11) = {width, height, 0, cellSize};
        Point(12) = {0, height + transitionRegionHeight, 0, width};
        Point(13) = {width, height + transitionRegionHeight, 0, width};
        Line(14) = {10,11};
        Line(15) = {11,13};
        Line(16) = {13,12};
        Line(17) = {12,10};
        Line Loop(18) = {14, 15, 16, 17};
        Plane Surface(19) = {18};

        Extrude{0, height, 0} {
            Line{3}; Layers{ ny }; Recombine;}

        Line(100) = {12, 13};
        Extrude{0, boundaryLayerHeight, 0} {
            Line{100}; Layers{ numberOfBoundaryLayerCells }; Recombine;}
        """ % locals(), communicator=communicator)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
