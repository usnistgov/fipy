#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "viewer.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
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

__all__ = ["AbstractViewer"]

import sys

class AbstractViewer(object):
    """
    .. attention:: This class is abstract. Always create one of its subclasses.
    """
    def __init__(self, vars, title=None, **kwlimits):
        """Create a `AbstractViewer` object.

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
        """
        if self.__class__ is AbstractViewer:
            raise NotImplementedError, "can't instantiate abstract base class"

        self.vars = self._getSuitableVars(vars)

        self.limits = kwlimits

        if title is None:
            if len(self.vars) == 1:
                title = self.vars[0].name
            else:
                title = ''

        self.title = title

    def _getSuitableVars(self, vars):
        if type(vars) not in [type([]), type(())]:
            vars = [vars]

        return [var for var in vars]

    def setLimits(self, limits={}, **kwlimits):
        """
        Update the limits.

        :Parameters:
          limits : dict
            a (deprecated) alternative to limit keyword arguments
          xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax
            displayed range of data. A 1D `Viewer` will only use `xmin` and
            `xmax`, a 2D viewer will also use `ymin` and `ymax`, and so on. All
            viewers will use `datamin` and `datamax`. Any limit set to a
            (default) value of `None` will autoscale.

        """
        self.limits.update(limits)
        self.limits.update(kwlimits)

    def _getLimit(self, keys, default=None):
        """
        Return the limit associated with the first available key in `keys`

        :Parameters:
          keys
            a `tuple`, `list`, or single key string that identifies
            the limit of interest

        :Returns:
          the value of the limit or `None`
        """
        if not (isinstance(keys, list) or isinstance(keys, tuple)):
            keys = (keys,)
        limit = default
        for key in keys:
            limit = self.limits.get(key, limit)
            if limit is not None:
                break

        return limit

    def plot(self, filename=None):
        """
        Update the display of the viewed variables.

        :Parameters:
          filename
            If not `None`, the name of a file to save the image into.
        """

        raise NotImplementedError

    def plotMesh(self, filename=None):
        """
        Display a representation of the mesh

        :Parameters:
          filename
            If not `None`, the name of a file to save the image into.
        """
        pass

    def _autoscale(self, vars, datamin=None, datamax=None):
        from fipy.tools import numerix

        if datamin is None:
            datamin = 1e300
            for var in vars:
                datamin = min(datamin, numerix.nanmin(var))

        if datamax is None:
            from fipy.tools import numerix
            datamax = -1e300
            for var in vars:
                datamax = max(datamax, numerix.nanmax(var))

        return datamin, datamax

    def _validFileExtensions(self):
        return []

    _saved_stdout = sys.stdout

    @staticmethod
    def _serial_doctest_raw_input(prompt):
        """Replacement for `raw_input()` that works in doctests
        """
        AbstractViewer._saved_stdout.write("\n")
        AbstractViewer._saved_stdout.write(prompt)
        AbstractViewer._saved_stdout.flush()
        return sys.stdin.readline()

    @staticmethod
    def _doctest_raw_input(prompt):
        """Replacement for `raw_input()` that works in doctests

        This routine attempts to be savvy about running in parallel.
        """
        try:
            from fipy.tools import parallelComm
            parallelComm.Barrier()
            AbstractViewer._saved_stdout.flush()
            if parallelComm.procID == 0:
                txt = AbstractViewer._serial_doctest_raw_input(prompt)
            else:
                txt = ""
            parallelComm.Barrier()
        except ImportError:
            txt = AbstractViewer._serial_doctest_raw_input(prompt)
        return txt


    def _promptForOpinion(self, prompt="Describe any problems with this figure or hit Return: "):
        # This method is usually invoked from a test, which can have a weird
        # state; In particular, it may have a special `raw_input` to prevent user
        # interaction during the test.

        opinion = self._doctest_raw_input(self.__class__.__name__ + ": " + prompt).strip()
        if len(opinion) > 0:
            extensions = ", ".join(self._validFileExtensions())
            if len(extensions) > 0:
                extensions = " (%s)" % extensions
            snapshot = self._doctest_raw_input("Enter a filename%s to save a snapshot (leave blank to skip): " % extensions).strip()
            self.plot(snapshot)
            print opinion


    def _test1D(**kwargs):
        s = """
            >>> from fipy import *
            >>> mesh = Grid1D(nx=100)
            >>> x, = mesh.cellCenters
            >>> xVar = CellVariable(mesh=mesh, name="x", value=x)
            >>> k = Variable(name="k", value=0.)
            >>> viewer = VIEWERSUBSTITUTION(vars=(numerix.sin(k * xVar), numerix.cos(k * xVar / numerix.pi)),
            ...                 limits={'xmin': 10, 'xmax': 90},
            ...                 datamin=-0.9, datamax=2.0,
            ...                 title="VIEWERSUBSTITUTION test")
            >>> for kval in numerix.arange(0,0.3,0.03):
            ...     k.setValue(kval)
            ...     viewer.plot()
            >>> viewer._promptForOpinion()
        """
        s = s.replace("VIEWERSUBSTITUTION", kwargs["viewer"])
        return s
    _test1D = staticmethod(_test1D)

    def _test2Dbase(**kwargs):
        s = """
            >>> from fipy import *
            >>> mesh = MESHSUBSTITUTION
            >>> x, y = mesh.cellCenters
            >>> xyVar = CellVariable(mesh=mesh, name="x y", value=x * y)
            >>> k = Variable(name="k", value=0.)
            >>> viewer = VIEWERSUBSTITUTION(vars=numerix.sin(k * xyVar),
            ...                 limits={'ymin': 0.1, 'ymax': 0.9},
            ...                 datamin=-0.9, datamax=2.0,
            ...                 title="VIEWERSUBSTITUTION test")
            >>> for kval in range(10):
            ...     k.setValue(kval)
            ...     viewer.plot()
            >>> viewer._promptForOpinion()
        """
        s = s.replace("VIEWERSUBSTITUTION", kwargs["viewer"])
        s = s.replace("MESHSUBSTITUTION", kwargs["mesh"])
        return s
    _test2Dbase = staticmethod(_test2Dbase)

    def _test2D(**kwargs):
        return AbstractViewer._test2Dbase(mesh="Grid2D(nx=50, ny=100, dx=0.1, dy=0.01)",
                                      **kwargs)
    _test2D = staticmethod(_test2D)

    def _test2Dirregular(**kwargs):
        """"""
        return AbstractViewer._test2Dbase(mesh="""(Grid2D(nx=5, ny=10, dx=0.1, dy=0.1)
            ...         + (Tri2D(nx=5, ny=5, dx=0.1, dy=0.1)
            ...          + ((0.5,), (0.2,))))""", **kwargs)
    _test2Dirregular = staticmethod(_test2Dirregular)

    def _test2DvectorBase(**kwargs):
        s = """
            >>> from fipy import *
            >>> mesh = MESHSUBSTITUTION
            >>> x, y = mesh.cellCenters
            >>> xyVar = CellVariable(mesh=mesh, name="x y", value=x * y)
            >>> k = Variable(name="k", value=1.)
            >>> viewer = VIEWERSUBSTITUTION(vars=numerix.sin(k * xyVar).grad,
            ...                 title="VIEWERSUBSTITUTION test")
            >>> for kval in numerix.arange(1, 10):
            ...     k.setValue(kval)
            ...     viewer.plot()
            >>> viewer._promptForOpinion()

            >>> viewer = VIEWERSUBSTITUTION(vars=numerix.sin(k * xyVar).faceGrad,
            ...                 title="VIEWERSUBSTITUTION test")
            >>> for kval in numerix.arange(1, 10):
            ...     k.setValue(kval)
            ...     viewer.plot()
            >>> viewer._promptForOpinion()
        """
        s = s.replace("VIEWERSUBSTITUTION", kwargs["viewer"])
        s = s.replace("MESHSUBSTITUTION", kwargs["mesh"])
        return s

    _test2DvectorBase = staticmethod(_test2DvectorBase)

    def _test2Dvector(**kwargs):
        return AbstractViewer._test2DvectorBase(mesh="Grid2D(nx=50, ny=100, dx=0.1, dy=0.01)",
                                            **kwargs)
    _test2Dvector = staticmethod(_test2Dvector)

    def _test2DvectorIrregular(**kwargs):
        """"""
        return AbstractViewer._test2DvectorBase(mesh="""(Grid2D(nx=5, ny=10, dx=0.1, dy=0.1)
            ...         + (Tri2D(nx=5, ny=5, dx=0.1, dy=0.1)
            ...          + ((0.5,), (0.2,))))""", **kwargs)
    _test2DvectorIrregular = staticmethod(_test2DvectorIrregular)


    def _test3D(**kwargs):
        s = """
            >>> from fipy import *
            >>> mesh = Grid3D(nx=50, ny=100, nz=10, dx=0.1, dy=0.01, dz=0.1)
            >>> x, y, z = mesh.cellCenters
            >>> xyzVar = CellVariable(mesh=mesh, name=r"x y z", value=x * y * z)
            >>> k = Variable(name="k", value=0.)
            >>> viewer = VIEWERSUBSTITUTION(vars=numerix.sin(k * xyzVar),
            ...                     limits={'ymin': 0.1, 'ymax': 0.9},
            ...                     datamin=-0.9, datamax=2.0,
            ...                     title="VIEWERSUBSTITUTION test")
            >>> for kval in range(10):
            ...     k.setValue(kval)
            ...     viewer.plot()
            >>> viewer._promptForOpinion()
        """
        s = s.replace("VIEWERSUBSTITUTION", kwargs["viewer"])
        return s

    _test3D = staticmethod(_test3D)
