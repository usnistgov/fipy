#!/usr/bin/env python

##
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "tools.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed at the National Institute of Standards
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # United States Code this software is not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.
 # NIST assumes no responsibility whatsoever for its use by whatsoever for its use by
 # other parties, and makes no guarantees, expressed or implied, about
 # its quality, reliability, or any other characteristic.  We would
 # appreciate acknowledgement if the software is used.
 #
 # To the extent that NIST may hold copyright in countries other than the
 # United States, you are hereby granted the non-exclusive irrevocable and
 # unconditional right to print, publish, prepare derivative works and
 # distribute this software, in any medium, or authorize others to do so on
 # your behalf, on a royalty-free basis throughout the world.
 #
 # You may improve, modify, and create derivative works of the software or
 # any portion of the software, and you may copy and distribute such
 # modifications or works.  Modified works should carry a notice stating
 # that you changed the software and should note the date and nature of any
 # such change.  Please explicitly acknowledge the National Institute of
 # Standards and Technology as the original source.
 #
 # This software can be redistributed and/or modified freely provided that
 # any derivative works bear some notice that they are derived from it, and
 # any modified versions bear some notice that they have been modified.
 # ========================================================================
 #
 # ###################################################################
 ##

"""Vector utility functions that are inexplicably absent from Numeric
"""

from fipy.tools import inline, numerix

__all__ = ["putAdd", "prune"]

# Factored out for fipy.variables.surfactantConvectionVariable._ConvectionCoeff
# for some reason
def _putAdd(vector, ids, additionVector, mask=False):
    """This is a temporary replacement for Numeric.put as it was not doing
    what we thought it was doing.
    """
    additionVector = numerix.array(additionVector)

    if numerix.sometrue(mask):
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
        """ This is a temporary replacement for Numeric.put as it was not doing
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
        """ This is a temporary replacement for Numeric.put as it was not doing
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
