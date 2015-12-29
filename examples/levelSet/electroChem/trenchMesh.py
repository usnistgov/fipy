#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 #
 #  FILE: "trenchMesh.py"
 #
 #  Author: Jonathan Guyer   <guyer@nist.gov>
 #  Author: Daniel Wheeler   <daniel.wheeler@nist.gov>
 #  Author: James Warren     <jwarren@nist.gov>
 #  Author: Andrew Acquaviva <andrewa@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  setup.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 #
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

__docformat__ = 'restructuredtext'

from fipy.tools import numerix
from fipy.variables.cellVariable import CellVariable
from fipy.tools import parallelComm
from gapFillMesh import GapFillMesh

class TrenchMesh(GapFillMesh):

    """
    The following test case tests for diffusion across the domain.

    >>> cellSize = 0.05e-6
    >>> trenchDepth = 0.5e-6
    >>> boundaryLayerDepth = 50e-6
    >>> domainHeight = 10 * cellSize + trenchDepth + boundaryLayerDepth

    >>> mesh = TrenchMesh(trenchSpacing = 1e-6,
    ...                   cellSize = cellSize,
    ...                   trenchDepth = trenchDepth,
    ...                   boundaryLayerDepth = boundaryLayerDepth,
    ...                   aspectRatio = 1.) # doctest: +GMSH

    >>> import fipy.tools.dump as dump
    >>> (f, filename) = dump.write(mesh) # doctest: +GMSH
    >>> if parallelComm.Nproc == 1:
    ...     mesh = dump.read(filename, f) # doctest: +GMSH
    >>> print mesh.globalNumberOfCells - len(numerix.nonzero(mesh.electrolyteMask)[0]) # doctest: +GMSH, +SERIAL
    150
    >>> print 400 < mesh.globalNumberOfCells < 800 # doctest: +GMSH
    True

    >>> from fipy.variables.cellVariable import CellVariable
    >>> var = CellVariable(mesh = mesh, value = 0.) # doctest: +GMSH

    >>> from fipy.terms.diffusionTerm import DiffusionTerm
    >>> eq = DiffusionTerm() # doctest: +GMSH

    >>> var.constrain(0., mesh.facesBottom) # doctest: +GMSH
    >>> var.constrain(domainHeight, mesh.facesTop) # doctest: +GMSH

    >>> eq.solve(var) # doctest: +GMSH

    Evaluate the result:

    >>> centers = mesh.cellCenters[1].copy() # doctest: +GMSH

    .. note:: the copy makes the array contiguous for inlining

    >>> localErrors = (centers - var)**2 / centers**2 # doctest: +GMSH
    >>> globalError = numerix.sqrt(numerix.sum(localErrors) / mesh.numberOfCells) # doctest: +GMSH
    >>> argmax = numerix.argmax(localErrors) # doctest: +GMSH
    >>> print numerix.sqrt(localErrors[argmax]) < 0.051 # doctest: +GMSH
    1
    >>> print globalError < 0.02 # doctest: +GMSH
    1

    """

    def __init__(self,
                 trenchDepth=None,
                 trenchSpacing=None,
                 boundaryLayerDepth=None,
                 cellSize=None,
                 aspectRatio=None,
                 angle=0.,
                 communicator=parallelComm):
        """

        `trenchDepth` - Depth of the trench.

        `trenchSpacing` - The distance between the trenches.

        `boundaryLayerDepth` - The depth of the hydrodynamic boundary
        layer.

        `cellSize` - The cell Size.

        `aspectRatio` - trenchDepth / trenchWidth

        `angle` - The angle for the taper of the trench.

        The trench mesh takes the parameters generally used to define
        a trench region and recasts then for the general
        `GapFillMesh`.

        """

        heightBelowTrench = cellSize * 10.

        heightAboveTrench = trenchDepth / 1.

        fineRegionHeight = heightBelowTrench + trenchDepth + heightAboveTrench
        transitionHeight = fineRegionHeight * 3.
        domainWidth = trenchSpacing / 2.
        domainHeight = heightBelowTrench + trenchDepth + boundaryLayerDepth

        super(TrenchMesh, self).__init__(cellSize=cellSize,
                                         desiredDomainWidth=domainWidth,
                                         desiredDomainHeight=domainHeight,
                                         desiredFineRegionHeight=fineRegionHeight,
                                         transitionRegionHeight=transitionHeight,
                                         communicator=parallelComm)

        trenchWidth = trenchDepth / aspectRatio

        x, y = self.cellCenters
        Y = (y - (heightBelowTrench + trenchDepth / 2))
        taper = numerix.tan(angle) * Y
        self.electrolyteMask = numerix.where(y > trenchDepth + heightBelowTrench,
                                             1,
                                             numerix.where(y < heightBelowTrench,
                                                           0,
                                                           numerix.where(x > trenchWidth / 2 + taper,
                                                                         0,
                                                                         1)))

    def __getstate__(self):
        dict = super(TrenchMesh, self).__getstate__()
        dict['electrolyteMask'] = self.electrolyteMask
        return dict

    def __setstate__(self, dict):
        self.electrolyteMask = dict['electrolyteMask']
        del dict['electrolyteMask']
        super(TrenchMesh, self).__setstate__(dict)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
