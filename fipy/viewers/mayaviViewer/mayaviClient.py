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

    def __init__(self, vars, title=None, daemon_file=None, fps=1.0, **kwlimits):
        """Create a `MayaviClient`.

        Parameters
        ----------
        vars : ~fipy.variables.cellVariable.CellVariable or list
            `CellVariable` objects to plot
        title : str, optional
            displayed at the top of the `Viewer` window
        xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax : float, optional
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
        self._fps = fps

        self._vtkdir = tempfile.mkdtemp()
        self._vtkcellfname = os.path.join(self._vtkdir, "cell.vtk")
        self._vtkfacefname = os.path.join(self._vtkdir, "face.vtk")
        self._vtklockfname = os.path.join(self._vtkdir, "lock")

        from fipy.viewers.vtkViewer import VTKCellViewer, VTKFaceViewer

        try:
            self._vtkCellViewer = VTKCellViewer(vars=vars)
            cell_vars = self._vtkCellViewer.vars
        except TypeError:
            self._vtkCellViewer = None
            cell_vars = []

        try:
            self._vtkFaceViewer = VTKFaceViewer(vars=vars)
            face_vars = self._vtkFaceViewer.vars
        except TypeError:
            self._vtkFaceViewer = None
            face_vars = []

        AbstractViewer.__init__(self, vars=cell_vars + face_vars, title=title, **kwlimits)

        self.plot()

        try:
            from importlib import resources
            
            ref = resources.files("fipy") / "viewers/mayaviViewer/mayaviDaemon.py"
            with resources.as_file(ref) as path:
                builtin_daemon = path.as_posix()
        except ImportError:
            from pkg_resources import Requirement, resource_filename
            builtin_daemon = resource_filename(Requirement.parse("FiPy"),
                                               "fipy/viewers/mayaviViewer/mayaviDaemon.py")

        daemon_file = (daemon_file or builtin_daemon)

        pyth = sys.executable or "python"

        cmd = [pyth,
               daemon_file,
               "--lock",
               self._vtklockfname,
               "--fps",
               str(self._fps)]

        if self._vtkCellViewer is not None:
            cmd += ["--cell", self._vtkcellfname]

        if self._vtkFaceViewer is not None:
            cmd += ["--face", self._vtkfacefname]


        cmd += self._getLimit('xmin')
        cmd += self._getLimit('xmax')
        cmd += self._getLimit('ymin')
        cmd += self._getLimit('ymax')
        cmd += self._getLimit('zmin')
        cmd += self._getLimit('zmax')
        cmd += self._getLimit('datamin')
        cmd += self._getLimit('datamax')

        self._daemon = subprocess.Popen(cmd)

    @property
    def fps(self):
        """The frames per second to attempt to display."""
        return self._fps

    def __del__(self):
        for fname in [self._vtkcellfname, self._vtkfacefname, self._vtklockfname]:
            if fname and os.path.isfile(fname):
                os.unlink(fname)
        os.rmdir(self._vtkdir)

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
            if not os.path.isfile(self._vtklockfname):
                if self._vtkCellViewer is not None:
                    self._vtkCellViewer.plot(filename=self._vtkcellfname)
                if self._vtkFaceViewer is not None:
                    self._vtkFaceViewer.plot(filename=self._vtkfacefname)
                with open(self._vtklockfname, 'w') as lock:
                    if filename is not None:
                        lock.write(filename)
                plotted = True

            if (time.time() - start > 30. / self._fps) and not plotted:
                print("viewer: NOT READY")
                start = time.time()
        if not plotted:
            print("viewer: SKIPPED")

    def _validFileExtensions(self):
        return [".png", ".jpg", ".bmp", ".tiff", ".ps", ".eps", ".pdf", ".rib", ".oogl", ".iv", ".vrml", ".obj"]

    @classmethod
    def _doctest_body(cls):
        return (cls._test1D()
                + cls._test2D()
                + cls._test2Dirregular()
                + cls._test3D())

    @classmethod
    def _doctest_extra(cls):
        return ""

if __name__ == "__main__":
    import fipy.tests.doctestPlus
    fipy.tests.doctestPlus.execButNoTest()
