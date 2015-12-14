#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "mayaviClient.py"
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

import os
import subprocess
import tempfile
import time

from fipy.viewers.viewer import AbstractViewer

__all__ = ["MayaviClient"]

class MayaviClient(AbstractViewer):
    """
    The `MayaviClient` uses the Mayavi_ python plotting package.

    .. _Mayavi: http://code.enthought.com/projects/mayavi

    """
    __doc__ += AbstractViewer._test1D(viewer="MayaviClient")
    __doc__ += AbstractViewer._test2D(viewer="MayaviClient")
    __doc__ += AbstractViewer._test2Dirregular(viewer="MayaviClient")
    __doc__ += AbstractViewer._test3D(viewer="MayaviClient")

    def __init__(self, vars, title=None, daemon_file=None, fps=1.0, **kwlimits):
        """
        Create a `MayaviClient`.

        :Parameters:
          vars
            a `CellVariable` or tuple of `CellVariable` objects to plot
          title
            displayed at the top of the `Viewer` window
          xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax
            displayed range of data. A 1D `Viewer` will only use `xmin` and
            `xmax`, a 2D viewer will also use `ymin` and `ymax`, and so on. All
            viewers will use `datamin` and `datamax`. Any limit set to a
            (default) value of `None` will autoscale.
          daemon_file
            the path to the script to run the separate MayaVi viewer process.
            Defaults to "fipy/viewers/mayaviViewer/mayaviDaemon.py"
          fps
            frames per second to attempt to display
        """
        self.fps = fps

        self.vtkdir = tempfile.mkdtemp()
        self.vtkcellfname = os.path.join(self.vtkdir, "cell.vtk")
        self.vtkfacefname = os.path.join(self.vtkdir, "face.vtk")
        self.vtklockfname = os.path.join(self.vtkdir, "lock")

        from fipy.viewers.vtkViewer import VTKCellViewer, VTKFaceViewer

        try:
            self.vtkCellViewer = VTKCellViewer(vars=vars)
            cell_vars = self.vtkCellViewer.vars
        except TypeError:
            self.vtkCellViewer = None
            cell_vars = []

        try:
            self.vtkFaceViewer = VTKFaceViewer(vars=vars)
            face_vars = self.vtkFaceViewer.vars
        except TypeError:
            self.vtkFaceViewer = None
            face_vars = []

        AbstractViewer.__init__(self, vars=cell_vars + face_vars, title=title, **kwlimits)

        self.plot()

        from pkg_resources import Requirement, resource_filename
        daemon_file = (daemon_file
                       or resource_filename(Requirement.parse("FiPy"),
                                            "fipy/viewers/mayaviViewer/mayaviDaemon.py"))

        cmd = ["python",
               daemon_file,
               "--lock",
               self.vtklockfname,
               "--fps",
               str(self.fps)]

        if self.vtkCellViewer is not None:
            cmd += ["--cell", self.vtkcellfname]

        if self.vtkFaceViewer is not None:
            cmd += ["--face", self.vtkfacefname]


        cmd += self._getLimit('xmin')
        cmd += self._getLimit('xmax')
        cmd += self._getLimit('ymin')
        cmd += self._getLimit('ymax')
        cmd += self._getLimit('zmin')
        cmd += self._getLimit('zmax')
        cmd += self._getLimit('datamin')
        cmd += self._getLimit('datamax')

        self.daemon = subprocess.Popen(cmd)

    def __del__(self):
        for fname in [self.vtkcellfname, self.vtkfacefname, self.vtklockfname]:
            if fname and os.path.isfile(fname):
                os.unlink(fname)
        os.rmdir(self.vtkdir)

    def _getLimit(self, key, default=None):
        """
        Return the limit associated with the key

        .. Note::

           `MayaviClient` does not need the generality of multiple keys
           because it is always 3D

        :Parameters:
          key
            a key string that identifies the limit of interest

        :Returns:
          the value of the limit or `None`
        """
        lim = AbstractViewer._getLimit(self, key, default=None)
        if lim is not None:
            return ["--%s" % key, str(lim)]
        else:
            return []

    def plot(self, filename=None):
        start = time.time()
        plotted = False
        while not plotted:
            if not os.path.isfile(self.vtklockfname):
                if self.vtkCellViewer is not None:
                    self.vtkCellViewer.plot(filename=self.vtkcellfname)
                if self.vtkFaceViewer is not None:
                    self.vtkFaceViewer.plot(filename=self.vtkfacefname)
                lock = file(self.vtklockfname, 'w')
                if filename is not None:
                    lock.write(filename)
                lock.close()
                plotted = True

            if (time.time() - start > 30. / self.fps) and not plotted:
                print "viewer: NOT READY"
                start = time.time()
        if not plotted:
            print "viewer: SKIPPED"

    def _validFileExtensions(self):
        return [".png",".jpg",".bmp",".tiff",".ps",".eps",".pdf",".rib",".oogl",".iv",".vrml",".obj"]

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
