"""Vector utility functions that are inexplicably absent from Numeric
"""
from __future__ import unicode_literals

from builtins import range
from builtins import zip
from fipy.tools import inline, numerix

__all__ = ["putAdd", "prune"]
from future.utils import text_to_native_str
__all__ = [text_to_native_str(n) for n in __all__]

# Factored out for fipy.variables.surfactantConvectionVariable._ConvectionCoeff
# for some reason
def _putAdd(vector, ids, additionVector, mask=False):
    """This is a temporary replacement for `Numeric.put` as it was not doing
    what we thought it was doing.
    """
    additionVector = numerix.array(additionVector)

    if numerix.any(mask):
        if len(vector.shape) < len(additionVector.shape):
            for j in range(vector.shape[0]):
                for id, value, masked in zip(ids.flat, additionVector[j].flat, mask.flat):
                    if not masked:
                        vector[j].flat[id] += value
        else:
            for id, value, masked in zip(ids.flat, additionVector.flat, mask.flat):
                if not masked:
                    vector.flat[id] += value

    else:
        if len(vector.shape) < len(additionVector.shape):
            for j in range(vector.shape[0]):
                for id, value in zip(ids.flat, additionVector[j].flat):
                    vector[j].flat[id] += value
        else:
            for id, value in zip(ids.flat, additionVector.flat):
                vector.flat[id] += value

if inline.doInline:
    ## FIXME: inline version doesn't account for all of the conditions that Python
    ## version does.
    def putAdd(vector, ids, additionVector):
        """ This is a temporary replacement for `Numeric.put` as it was not doing
        what we thought it was doing.
        """
        inline._runInline("""
                              int ID = ids[i];
                              vector[ID] += additionVector[i];
                          """,
                          vector=vector,
                          ids=ids,
                          additionVector=numerix.array(additionVector),
        ni = len(ids.flat))
else:
    def putAdd(vector, ids, additionVector):
        """ This is a temporary replacement for `Numeric.put` as it was not doing
        what we thought it was doing.
        """
        _putAdd(vector, ids, additionVector)

def prune(array, shift, start=0, axis=0):
    """
    removes elements with indices i = start + shift * n
    where n = 0, 1, 2, ...

        >>> prune(numerix.arange(10), 3, 5)
        array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        >>> prune(numerix.arange(10), 3, 2)
        array([0, 1, 3, 4, 6, 7, 9])
        >>> prune(numerix.arange(10), 3)
        array([1, 2, 4, 5, 7, 8])
        >>> prune(numerix.arange(4, 7), 3)
        array([5, 6])

    """

    takeArray = numerix.nonzero(numerix.arange(array.shape[-1]) % shift != start)[0]
    return numerix.take(array, takeArray, axis=axis)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
