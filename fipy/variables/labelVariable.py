from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable
from fipy.tools import numerix

__all__ = ["LabelVariable"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LabelVariable(CellVariable):
    """Label features in `var` using scipy.ndimage.label
    
    >>> import fipy as fp

    Create a 1D domain with two distinct non-zero regions.

    >>> mesh = fp.Grid1D(nx=5)

    >>> features = fp.CellVariable(mesh=mesh, name="features")
    >>> features.setValue(1., where=(mesh.x > 1) & (mesh.x < 3))
    >>> features.setValue(0.5, where=(mesh.x > 4) & (mesh.x < 5))
    >>> print(features.value)
    [ 0.   1.   1.   0.   0.5]

    Label the non-zero regions.

    >>> labels = fp.LabelVariable(features, name="labels")
    >>> print(labels.num_features)
    2
    >>> print(labels)
    [0 1 1 0 2]

    To have no connectivity between cells, we can assign a structure with
    no neighbors.

    >>> labels = fp.LabelVariable(features, structure=[0,1,0],
    ...                           name="labels")
    >>> print(labels.num_features)
    3
    >>> print(labels)
    [0 1 2 0 3]

    Similarly, create a 2D domain with two distinct non-zero regions.

    >>> mesh = fp.Grid2D(nx=5, ny=5)

    >>> features = fp.CellVariable(mesh=mesh, name="features")
    >>> features.setValue(1., where=((mesh.x > 1) & (mesh.x < 3)
    ...                              & (mesh.y > 1) & (mesh.y < 3)))
    >>> features.setValue(0.5, where=((mesh.x > 3) & (mesh.x < 4)
    ...                               & (mesh.y > 3) & (mesh.y < 4)))

    Note that FiPy arranges its gridded cells in a right-handed Cartesian
    fashion, with x increasing to the right, y increasing up, and z
    increasing toward you.  Conversely, NumPy arrays and ndimages are
    arranged in stacks of increasing z, each consisting of rows of
    increasing y, each element of which increases in x.  As a result, we
    reverse FiPy's (x,y) shape to NumPy's (y, x):

    >>> print(features.value.reshape(mesh.shape[::-1]))
    [[ 0.   0.   0.   0.   0. ]
     [ 0.   1.   1.   0.   0. ]
     [ 0.   1.   1.   0.   0. ]
     [ 0.   0.   0.   0.5  0. ]
     [ 0.   0.   0.   0.   0. ]]

    >>> labels = fp.LabelVariable(features, name="labels")
    >>> print(labels.num_features)
    2
    >>> print(labels.value.reshape(mesh.shape[::-1]))
    [[0 0 0 0 0]
     [0 1 1 0 0]
     [0 1 1 0 0]
     [0 0 0 2 0]
     [0 0 0 0 0]]

    By default, the two domains are seen as unconnected because there is no
    overlap of cells along either the x or y axis.  The following structure
    creates connectivity along one diagonal, but not the other.  Note that
    the structure array follows NumPy (y, x) ordering, rather than FiPy (x,
    y) ordering.

    >>> labels = fp.LabelVariable(features, structure=[[1,1,0],
    ...                                                [1,1,1],
    ...                                                [0,1,1]],
    ...                             name="labels")
    >>> print(labels.num_features)
    1
    >>> print(labels.value.reshape(mesh.shape[::-1]))
    [[0 0 0 0 0]
     [0 1 1 0 0]
     [0 1 1 0 0]
     [0 0 0 1 0]
     [0 0 0 0 0]]

    Similarly, create a 3D domain with three distinct non-zero regions.

    >>> mesh = fp.Grid3D(nx=3, ny=3, nz=3)

    >>> features = fp.CellVariable(mesh=mesh, name="features")
    >>> features.setValue(1., where=((mesh.x > 0) & (mesh.x < 2)
    ...                              & (mesh.y > 0) & (mesh.y < 2)
    ...                              & (mesh.z > 0) & (mesh.z < 2)))
    >>> features.setValue(0.7, where=((mesh.x > 2) & (mesh.x < 3)
    ...                               & (mesh.y > 2) & (mesh.y < 3)
    ...                               & (mesh.z > 0) & (mesh.z < 1)))
    >>> features.setValue(0.5, where=((mesh.x > 2) & (mesh.x < 3)
    ...                               & (mesh.y > 2) & (mesh.y < 3)
    ...                               & (mesh.z > 2) & (mesh.z < 3)))

    We reverse FiPy's (x,y, z) shape to NumPy's (z, y, x)

    >>> print(features.value.reshape(mesh.shape[::-1]))
    [[[ 1.   1.   0. ]
      [ 1.   1.   0. ]
      [ 0.   0.   0.7]]
    <BLANKLINE>
     [[ 1.   1.   0. ]
      [ 1.   1.   0. ]
      [ 0.   0.   0. ]]
    <BLANKLINE>
     [[ 0.   0.   0. ]
      [ 0.   0.   0. ]
      [ 0.   0.   0.5]]]

    >>> labels = fp.LabelVariable(features, name="labels")
    >>> print(labels.num_features)
    3
    >>> print(labels.value.reshape(mesh.shape[::-1]))
    [[[1 1 0]
      [1 1 0]
      [0 0 2]]
    <BLANKLINE>
     [[1 1 0]
      [1 1 0]
      [0 0 0]]
    <BLANKLINE>
     [[0 0 0]
      [0 0 0]
      [0 0 3]]]

    By default, the three domains are seen as unconnected because there is
    no overlap of cells along any of the x, y, or z axes.  The following
    structure creates connectivity along one major diagonal, but not the
    other two.  Note that the structure array follows NumPy (z, y, x)
    ordering, rather than FiPy (x, y, z) ordering.

    >>> labels = fp.LabelVariable(features, structure=[[[1,0,0],
    ...                                                 [0,1,0],
    ...                                                 [0,0,0]],
    ...                                                [[0,1,0],
    ...                                                 [1,1,1],
    ...                                                 [0,1,0]],
    ...                                                [[0,0,0],
    ...                                                 [0,1,0],
    ...                                                 [0,0,1]]],
    ...                         name="labels")
    >>> print(labels.num_features)
    2
    >>> print(labels.value.reshape(mesh.shape[::-1]))
    [[[1 1 0]
      [1 1 0]
      [0 0 2]]
    <BLANKLINE>
     [[1 1 0]
      [1 1 0]
      [0 0 0]]
    <BLANKLINE>
     [[0 0 0]
      [0 0 0]
      [0 0 1]]]

    Parameters
    ----------
    var : ~fipy.variables.cellVariable.CellVariable
        Field to be labeled.  Any non-zero values in input are counted as
        features and zero values are considered the background.
        
        .. important:
           Only sensible if `var` is defined on a `...Grid...` Mesh.
    structure : array_like, optional
        A structuring element that defines feature connections.
        `structure` must be centrosymmetric
        (see ```scipy.ndimage.label`` Notes
        <https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.label.html#scipy.ndimage.label>`_).
        If no structuring element is provided,
        one is automatically generated with a squared connectivity equal to
        one.  That is, for a 2-D `input` array, the default structuring element
        is::
            [[0,1,0],
             [1,1,1],
             [0,1,0]]
    dtype : date-type, optional
        The desired data-type for the labels. Note that the type must be able
        to store the largest label, or this Variable will raise an Exception.
        Default: int.
    """
    def __init__(self, var, name="", structure=None, dtype=int):
        # We want our value to hold dtype,
        # but if we pass an array, the CellVariable
        # will probably be wonky
        value = numerix.array(0.).astype(dtype).item()
        CellVariable.__init__(self,
                              mesh=var.mesh,
                              name=name,
                              value=value,
                              elementshape=var.shape[:-1])
        self.var = self._requires(var)
        self.structure = structure
        self.dtype = dtype
        self._num_features = None
    
    def _calcValue(self):
        """Label features of `var`
        
        Side-effect: sets self._num_features
        """
        from scipy import ndimage

        feat = self.var.globalValue
        feat = feat.reshape(self.var.mesh.shape[::-1])

        arr = numerix.empty(self.var.mesh.shape[::-1], dtype=self.dtype)
        self._num_features = ndimage.label(input=feat,
                                           structure=self.structure,
                                           output=arr)
        return arr.flatten()
        
    @property
    def num_features(self):
        """How many objects were found
        """
        if self.stale or not self._isCached() or self._num_features is None:
            self._getValue()

        return self._num_features

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
