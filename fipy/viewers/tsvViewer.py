from __future__ import unicode_literals
from builtins import range
__docformat__ = 'restructuredtext'

import sys

from fipy.tools import numerix
from fipy.viewers.viewer import AbstractViewer
from fipy.variables.cellVariable import CellVariable
from fipy.variables.faceVariable import FaceVariable

__all__ = ["TSVViewer"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class TSVViewer(AbstractViewer):
    """
    "Views" one or more variables in tab-separated-value format.

    Output is a list of coordinates and variable values at each cell center.

    File contents will be, e.g.::

        title
        x	y	...	var0	var2	...
        0.0	0.0	...	3.14	1.41	...
        1.0	0.0	...	2.72	0.866	...
        :
        :

    """
    _axis = ["x", "y", "z"]

    def __init__(self, vars, title=None, limits={}, **kwlimits):
        """
        Creates a `TSVViewer`.

        Any cell centers that lie outside the limits provided will not be included.
        Any values that lie outside the `datamin` or `datamax` will be
        replaced with `nan`.

        All variables must have the same mesh.

        It tries to do something reasonable with rank-1 `CellVariable` and `FaceVariable` objects.

        Parameters
        ----------
        vars : ~fipy.variables.cellVariable.CellVariable or ~fipy.variables.faceVariable.FaceVariable or list
            the :class:`~fipy.variables.meshVariable.MeshVariable` objects to display.
        title : str, optional
            displayed at the top of the `Viewer` window
        limits : dict, optional
            a (deprecated) alternative to limit keyword arguments
        xmin, xmax, ymin, ymax, zmin, zmax, datamin, datamax : float, optional
            displayed range of data. Any limit set to
            a (default) value of `None` will autoscale.
        """
        kwlimits.update(limits)
        AbstractViewer.__init__(self, vars=vars, title=title, **kwlimits)

        mesh = self.vars[0].mesh

        for var in self.vars:
            assert mesh is var.mesh


    def _plot(self, values, f, dim):
        for index in range(values.shape[-1]):
            lineValues = values[..., index]

            # omit any elements whose cell centers lie outside of the specified limits
            skip = False
            for axis in range(dim):
                mini = self._getLimit("%smin" % self._axis[axis])
                maxi = self._getLimit("%smax" % self._axis[axis])

                if (mini and lineValues[axis] < mini) or (maxi and lineValues[axis] > maxi):
                    skip = True
                    break

            if skip:
                continue

            # replace any values that lie outside of the specified datalimits with 'nan'
            for valIndex in range(dim, len(lineValues)):
                mini = self._getLimit("datamin")
                maxi = self._getLimit("datamax")

                if (mini and lineValues[valIndex] < mini) or (maxi and lineValues[valIndex] > maxi):
                    lineValues[valIndex] = float("NaN")
                    break

            line = ["%.15g" % value for value in lineValues]
            f.write("\t".join(line))
            f.write("\n")

    def plot(self, filename=None):
        """
        "plot" the coordinates and values of the variables to `filename`.
        If `filename` is not provided, "plots" to `stdout`.

        >>> from fipy.meshes import Grid1D
        >>> m = Grid1D(nx = 3, dx = 0.4)
        >>> from fipy.variables.cellVariable import CellVariable
        >>> v = CellVariable(mesh = m, name = "var", value = (0, 2, 5))
        >>> TSVViewer(vars = (v, v.grad)).plot() #doctest: +NORMALIZE_WHITESPACE
        x       var     var_gauss_grad_x
        0.2     0       2.5
        0.6     2       6.25
        1       5       3.75

        >>> from fipy.meshes import Grid2D
        >>> m = Grid2D(nx = 2, dx = .1, ny = 2, dy = 0.3)
        >>> v = CellVariable(mesh = m, name = "var", value = (0, 2, -2, 5))
        >>> TSVViewer(vars = (v, v.grad)).plot() #doctest: +NORMALIZE_WHITESPACE
        x       y       var     var_gauss_grad_x        var_gauss_grad_y
        0.05    0.15    0       10      -3.33333333333333
        0.15    0.15    2       10      5
        0.05    0.45    -2      35      -3.33333333333333
        0.15    0.45    5       35      5

        Parameters
        ----------
        filename : str
            If not `None`, the name of a file to save the image into.
        """

        mesh = self.vars[0].mesh
        dim = mesh.dim

        if filename is not None:
            import os
            if mesh.communicator.procID == 0:
                if os.path.splitext(filename)[1] == ".gz":
                    import gzip
                    f = gzip.GzipFile(filename = filename, mode = 'w', fileobj = None)
                else:
                    f = open(filename, "w")
            else:
                f = open(os.devnull, mode='w')
        else:
            f = sys.stdout

        if self.title and len(self.title) > 0:
            f.write(self.title)
            f.write("\n")

        headings = []
        for index in range(dim):
            headings.extend(self._axis[index])

        for var in self.vars:
            name = var.name
            if (isinstance(var, CellVariable) or isinstance(var, FaceVariable)) and var.rank == 1:
                for index in range(dim):
                    headings.extend(["%s_%s" % (name, self._axis[index])])
            else:
                headings.extend([name])

        f.write("\t".join(headings))
        f.write("\n")

        cellVars = [var for var in self.vars if isinstance(var, CellVariable)]
        faceVars = [var for var in self.vars if isinstance(var, FaceVariable)]

        if len(cellVars) > 0:
            values = mesh.cellCenters.globalValue
            for var in self.vars:
                if isinstance(var, CellVariable) and var.rank == 1:
                    values = numerix.concatenate((values, numerix.array(var.globalValue)))
                else:
                    values = numerix.concatenate((values, (numerix.array(var.globalValue),)))

            self._plot(values, f, dim)

        if len(faceVars) > 0:
            values = mesh.faceCenters.globalValue
            for var in self.vars:
                if isinstance(var, FaceVariable) and var.rank == 1:
                    values = numerix.concatenate((values, numerix.array(var.globalValue)))
                else:
                    values = numerix.concatenate((values, (numerix.array(var.globalValue),)))

            self._plot(values, f, dim)

        if f is not sys.stdout:
            f.close()

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
