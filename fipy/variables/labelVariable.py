from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

from fipy.variables.cellVariable import CellVariable
from fipy.tools import numerix

__all__ = ["LabelVariable"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

class LabelVariable(CellVariable):
    """Label features in `var` using scipy.ndimage.label
    
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
        arr = self.var.globalValue.astype(self.dtype)
        shape = (self.var.mesh.args['nx'], self.var.mesh.args['ny'])
        arr = arr.reshape(shape)
        self._num_features = ndimage.label(input=arr,
                                           structure=self.structure,
                                           output=arr)
        return arr.flat
        
    @property
    def num_features(self):
        """How many objects were found
        """
        if self.stale or not self._isCached() or self._num_features is None:
            self._getValue()

        return self._num_features
