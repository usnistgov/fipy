from __future__ import print_function
from __future__ import unicode_literals
from builtins import str
__docformat__ = 'restructuredtext'

import os
import subprocess
import sys
import tempfile
import time

from fipy.viewers.viewer import AbstractViewer

__all__ = ["MayaviClient"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

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
        """Create a `MayaviClient`.

        Parameters
        ----------
        vars : ~fipy.variables.cellVariable.CellVariable or list
            `CellVariable` objects to plot
        title : str, optional
            displayed at the top of the `Viewer` window
        float xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax : float, optional
            displayed range of data. A 1D `Viewer` will only use `xmin` and
            `xmax`, a 2D viewer will also use `ymin` and `ymax`, and so on. All
            viewers will use `datamin` and `datamax`. Any limit set to a
            (default) value of `None` will autoscale.
        daemon_file : str, optional
            the path to the script to run the separate Mayavi viewer process.
            Defaults to `fipy/viewers/mayaviViewer/mayaviDaemon.py`
        fps : float, optional
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

        pyth = sys.executable or "python"

        cmd = [pyth,
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

        Parameters
        ----------
        key : str
            dictionary key that identifies the limit of interest

        Returns
        -------
        float or None
            the value of the limit
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
                print("viewer: NOT READY")
                start = time.time()
        if not plotted:
            print("viewer: SKIPPED")

    def _validFileExtensions(self):
        return [".png", ".jpg", ".bmp", ".tiff", ".ps", ".eps", ".pdf", ".rib", ".oogl", ".iv", ".vrml", ".obj"]

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
