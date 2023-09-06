from __future__ import print_function
from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext'

__all__ = ["AbstractViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

import sys

from future.utils import with_metaclass

# Pattern for accessing classmethods on construction
# https://stackoverflow.com/a/13901161/2019542
# adapted for future
# https://python-future.org/compatible_idioms.html#metaclasses
class _MetaViewer(type):
    def __new__(mcl, classname, bases, classdict):
        """Create new viewer with class-appropriate doctests"""
        cls = type.__new__(mcl, classname, bases, classdict)

        doc = cls._doctest_body()
        if cls.__doc__ is not None:
            doc = cls.__doc__ + doc
        cls.__doc__ = doc

        return cls

class AbstractViewer(with_metaclass(_MetaViewer, object)):
    """Base class for FiPy Viewers

    .. attention:: This class is abstract. Always create one of its subclasses.
    """

    def __init__(self, vars, title=None, **kwlimits):
        """Create a `AbstractViewer` object.

        Parameters
        ----------
        vars : ~fipy.variables.cellVariable.CellVariable or list
            the `CellVariable` objects to display.
        title : str, optional
            displayed at the top of the `Viewer` window
        xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax : float, optional
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
        """
        if self.__class__ is AbstractViewer:
            raise NotImplementedError("can't instantiate abstract base class")

        self._vars = self._getSuitableVars(vars)

        self._limits = kwlimits

        if title is None:
            if len(self.vars) == 1:
                title = self.vars[0].name
            else:
                title = ''

        self.title = title

    @property
    def limits(self):
        return self._limits

    @property
    def title(self):
        """The text appearing at the top center.

        (default: if ``len(self.vars) == 1``, the name of
        the only :class:`~fipy.variables.variable.Variable`,
        otherwise ``""``.)
        """
        return self._title

    @title.setter
    def title(self, value):
        self._title = value

    @property
    def vars(self):
        """The :class:`~fipy.variables.variable.Variable` or list
        of :class:`~fipy.variables.variable.Variable` objects to
        display."""
        return self._vars

    def _getSuitableVars(self, vars):
        if type(vars) not in [type([]), type(())]:
            vars = [vars]

        return [var for var in vars]

    def setLimits(self, limits={}, **kwlimits):
        """
        Update the limits.

        Parameters
        ----------
        limits : dict, optional
            a (deprecated) alternative to limit keyword arguments
        xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax : float, optional
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

        Parameters
        ----------
        keys : str or :obj:`tuple` of :obj:`str`
            dictionary keys that identify the limits of interest

        Returns
        -------
        float or None
            the value of the limit
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

        Parameters
        ----------
        filename : str
            If not `None`, the name of a file to save the image into.
        """

        raise NotImplementedError

    def plotMesh(self, filename=None):
        """
        Display a representation of the mesh

        Parameters
        ----------
        filename : str
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

    @classmethod
    def _doctest_body(cls):
        return ""

    @classmethod
    def _doctest_extra(cls):
        return """
            >>> viewer.title = "{cls.__name__} changed"
            >>> viewer.plot()
            >>> viewer._promptForOpinion()
        """.format(**locals())

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
            print(opinion)

    @classmethod
    def _test1D(cls):
        return """
            >>> import fipy as fp
            >>> mesh = fp.Grid1D(nx=100)
            >>> x, = mesh.cellCenters
            >>> xVar = fp.CellVariable(mesh=mesh, name="x", value=x)
            >>> k = fp.Variable(name="k", value=0.)
            >>> viewer = {cls.__name__}(vars=(fp.numerix.sin(k * xVar) + 2,
            ...                               fp.numerix.cos(k * xVar / fp.numerix.pi) + 2),
            ...                 xmin=10, xmax=90,
            ...                 datamin=1.1, datamax=4.0,
            ...                 title="{cls.__name__} test")
            >>> for kval in fp.numerix.arange(0, 0.3, 0.03):
            ...     k.setValue(kval)
            ...     viewer.plot()
            >>> viewer._promptForOpinion()
        """.format(**locals()) + cls._doctest_extra()

    @classmethod
    def _test2Dbase(cls, mesh):
        return """
            >>> import fipy as fp
            >>> mesh = {mesh}
            >>> x, y = mesh.cellCenters
            >>> xyVar = fp.CellVariable(mesh=mesh, name="x y", value=x * y)
            >>> k = fp.Variable(name="k", value=0.)
            >>> viewer = {cls.__name__}(vars=fp.numerix.sin(k * xyVar) * 1000 + 1002,
            ...                 ymin=0.1, ymax=0.9,
            ...                 # datamin=1.1, datamax=4.0,
            ...                 title="{cls.__name__} test")
            >>> from builtins import range
            >>> for kval in range(10):
            ...     k.setValue(kval)
            ...     viewer.plot()
            >>> viewer._promptForOpinion()
        """.format(**locals()) + cls._doctest_extra()

    @classmethod
    def _test2D(cls):
        return cls._test2Dbase(mesh="fp.Grid2D(nx=50, ny=100, dx=0.1, dy=0.01)")

    @classmethod
    def _test2Dirregular(cls):
        """"""
        return cls._test2Dbase(mesh="""(fp.Grid2D(nx=5, ny=10, dx=0.1, dy=0.1)
            ...         + (fp.Tri2D(nx=5, ny=5, dx=0.1, dy=0.1)
            ...          + ((0.5,), (0.2,))))""")

    @classmethod
    def _test2DvectorBase(cls, mesh):
        return """
            >>> import fipy as fp
            >>> mesh = {mesh}
            >>> x, y = mesh.cellCenters
            >>> xyVar = fp.CellVariable(mesh=mesh, name="x y", value=x * y)
            >>> k = fp.Variable(name="k", value=1.)
            >>> viewer = {cls.__name__}(vars=fp.numerix.sin(k * xyVar).grad,
            ...                 title="{cls.__name__} test")
            >>> for kval in fp.numerix.arange(1, 10):
            ...     k.setValue(kval)
            ...     viewer.plot()
            >>> viewer._promptForOpinion()

            >>> viewer = {cls.__name__}(vars=fp.numerix.sin(k * xyVar).faceGrad,
            ...                 title="{cls.__name__} test")
            >>> from builtins import range
            >>> for kval in range(10):
            ...     k.setValue(kval)
            ...     viewer.plot()
            >>> viewer._promptForOpinion()
        """.format(**locals()) + cls._doctest_extra()

    @classmethod
    def _test2Dvector(cls):
        return cls._test2DvectorBase(mesh="fp.Grid2D(nx=50, ny=100, dx=0.1, dy=0.01)")

    @classmethod
    def _test2DvectorIrregular(cls):
        return cls._test2DvectorBase(mesh="""(fp.Grid2D(nx=5, ny=10, dx=0.1, dy=0.1)
            ...         + (fp.Tri2D(nx=5, ny=5, dx=0.1, dy=0.1)
            ...          + ((0.5,), (0.2,))))""")

    @classmethod
    def _test3D(cls):
        return """
            >>> import fipy as fp
            >>> mesh = fp.Grid3D(nx=50, ny=100, nz=10, dx=0.1, dy=0.01, dz=0.1)
            >>> x, y, z = mesh.cellCenters
            >>> xyzVar = fp.CellVariable(mesh=mesh, name=r"x y z", value=x * y * z)
            >>> k = fp.Variable(name="k", value=0.)
            >>> viewer = {cls.__name__}(vars=fp.numerix.sin(k * xyzVar) + 2,
            ...                     ymin=0.1, ymax=0.9,
            ...                     datamin=1.1, datamax=4.0,
            ...                     title="{cls.__name__} test")
            >>> from builtins import range
            >>> for kval in range(10):
            ...     k.setValue(kval)
            ...     viewer.plot()
            >>> viewer._promptForOpinion()
        """.format(**locals()) + cls._doctest_extra()
