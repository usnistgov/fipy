#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 # 
 #  FILE: "vtkViewer.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Stiles  <daniel.stiles@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this software is not subject to copyright
 # protection and is in the public domain.  FiPy is an experimental
 # system.  NIST assumes no responsibility whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 # 
 # This software can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 # ###################################################################
 ##


__docformat__ = 'restructuredtext'

from fipy.variables.faceVariable import FaceVariable
from fipy.viewers.viewer import _Viewer

from vtkViewer import _VTKViewer

class VTKFaceViewer(_VTKViewer):
    """Renders `_MeshVariable` data in VTK format
    """
    def _makeDataSet(self, mesh):
        return mesh.getVTKFaceDataSet()
        
    def _getData(self):
        return self.dataset.point_data

    def _getVariableClass(self):
        return FaceVariable
        
    def _test(self):
        """
        >>> import os
        >>> from tempfile import mkstemp
        >>> f, fname = mkstemp(".vtk")
        >>> os.close(f)

        >>> from enthought.mayavi.sources.vtk_file_reader import VTKFileReader

        >>> from fipy import *
        >>> from fipy.viewers.vtkViewer import VTKFaceViewer

        >>> m = Grid1D(nx=10)
        >>> x, = m.getCellCenters()
        >>> v1 = CellVariable(mesh=m, value=x*x, name="x*x")
        >>> v2 = CellVariable(mesh=m, value=x)
        >>> v3 = v1.getFaceGrad()
        >>> v3.name = "v1.getFaceGrad()"
        >>> v4 = v1.getHarmonicFaceValue()
        >>> v4.name = "v1.getHarmonicFaceValue()"
        >>> v5 = v1.getArithmeticFaceValue()
        >>> v5.name = "v1.getArithmeticFaceValue()"
        >>> VTKFaceViewer(vars=(v3, v4, v5)).plot(fname)
        >>> VTKFileReader().initialize(fname)

        >>> m = Grid2D(nx=1, ny=2)
        >>> x, y = m.getCellCenters()
        >>> v1 = CellVariable(mesh=m, value=x*y, name="x*y")
        >>> v2 = CellVariable(mesh=m, value=x*x) #, name="v2")
        >>> v3 = v1.getFaceGrad()
        >>> v3.name = "v1.getFaceGrad()"
        >>> v4 = v1.getHarmonicFaceValue()
        >>> v4.name = "v1.getHarmonicFaceValue()"
        >>> v5 = v1.getArithmeticFaceValue()
        >>> v5.name = "v1.getArithmeticFaceValue()"
        >>> VTKFaceViewer(vars=(v3, v4, v5)).plot(fname)
        >>> VTKFileReader().initialize(fname)

        >>> m = (Grid2D(nx=5, ny=10, dx=0.1, dy=0.1)
        ...      + (Tri2D(nx=5, ny=5, dx=0.1, dy=0.1))
        ...      + ((0.5,), (0.2,)))
        >>> x, y = m.getCellCenters()
        >>> v1 = CellVariable(mesh=m, value=x*y, name="x*y")
        >>> v2 = CellVariable(mesh=m, value=x*x) #, name="v2")
        >>> v3 = v1.getFaceGrad()
        >>> v3.name = "v1.getFaceGrad()"
        >>> v4 = v1.getHarmonicFaceValue()
        >>> v4.name = "v1.getHarmonicFaceValue()"
        >>> v5 = v1.getArithmeticFaceValue()
        >>> v5.name = "v1.getArithmeticFaceValue()"
        >>> VTKFaceViewer(vars=(v3, v4, v5)).plot(fname)
        >>> VTKFileReader().initialize(fname)

        >>> m = Grid3D(nx=2, ny=1, nz=1)
        >>> x, y, z = m.getCellCenters()
        >>> v1 = CellVariable(mesh=m, value=x*y*z, name="x*y*z")
        >>> v2 = CellVariable(mesh=m, value=x*y*y, name="x*y*y")
        >>> v3 = v1.getFaceGrad()
        >>> v3.name = "v1.getFaceGrad()"
        >>> v4 = v1.getHarmonicFaceValue()
        >>> v4.name = "v1.getHarmonicFaceValue()"
        >>> v5 = v1.getArithmeticFaceValue()
        >>> v5.name = "v1.getArithmeticFaceValue()"
        >>> VTKFaceViewer(vars=(v3, v4, v5)).plot(fname)
        >>> VTKFileReader().initialize(fname)

        >>> os.remove(fname)
        """


if __name__ == "__main__": 
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
