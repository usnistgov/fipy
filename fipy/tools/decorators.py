#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "decorators.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: James O'Beirne <james.obeirne@gmail.com>
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
 # Portions of this code are copied and/or derived from numpy.lib.utils
 #
 #   http://numpy.scipy.org/
 #   https://github.com/numpy/numpy/blob/master/numpy/lib/utils.py
 #
 # Copyright (c) 2005-2010, NumPy Developers.
 # All rights reserved.
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions are
 # met:
 #
 #     * Redistributions of source code must retain the above copyright
 #        notice, this list of conditions and the following disclaimer.
 #
 #     * Redistributions in binary form must reproduce the above
 #        copyright notice, this list of conditions and the following
 #        disclaimer in the documentation and/or other materials provided
 #        with the distribution.
 #
 #     * Neither the name of the NumPy Developers nor the names of any
 #        contributors may be used to endorse or promote products derived
 #        from this software without specific prior written permission.
 #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 # "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 # LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 # A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 # OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 # LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 # DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 # THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 # OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 # ###################################################################
 ##

import re
import sys
import warnings

__all__ = ["deprecate"]

# Stolen from `numpy.lib.utils`
if sys.version_info < (2, 4):
    # Can't set __name__ in 2.3
    import new
    def _set_function_name(func, name):
        func = new.function(func.__code__, func.__globals__,
                            name, func.__defaults__, func.__closure__)
        return func
else:
    def _set_function_name(func, name):
        func.__name__ = name
        return func

class _Deprecate(object):
    """
    Decorator class to deprecate old functions.

    Refer to `deprecate` for details.

    Stolen from `numpy.lib.utils`

    See Also
    --------
    numerix.deprecate

    """
    def __init__(self, old_name=None, new_name=None, message=None,
                 old_string=":func:`%s` is deprecated",
                 new_string="use :func:`%s` instead",
                 version="UNKNOWN"):
        self.old_name = old_name
        self.new_name = new_name
        self.old_string = old_string
        self.new_string = new_string
        self.message = message
        self.version = version

    def old_name_from_func(self, func):
        old_name = self.old_name
        if old_name is None:
            try:
                old_name = func.__name__
            except AttributeError:
                old_name = func.__name__
        return old_name

    def new_name_old_name(self, old_name):
        return self.new_name

    def __call__(self, func, *args, **kwargs):
        """
        Decorator call.  Refer to ``decorate``.

        """
        new_name = self.new_name
        message = self.message

        old_name = self.old_name_from_func(func=func)
        new_name = self.new_name_old_name(old_name=old_name)

        if new_name is None:
            depwarn = (self.old_string + "!") % old_name
            depdoc = ""
        else:
            depwarn = (self.old_string + ", " + self.new_string + "!") % \
                       (old_name, new_name)
            depdoc = self.new_string % new_name

        if message is not None:
            depwarn += "\n" + message
            depdoc += "\n" + message

        def newfunc(*args,**kwds):
            """:func:`arrayrange` is deprecated, use :func:`arange` instead!"""
            warnings.warn(depwarn, DeprecationWarning, stacklevel=2)
            return func(*args, **kwds)

        newfunc = _set_function_name(newfunc, old_name)

        depdoc = (["", "", ".. deprecated:: %s" % self.version]
                  + ["   " + s for s in depdoc.split('\n')]
                  + ["", ""])
        doc = func.__doc__
        if doc is None:
            doc = depdoc
        else:
            from textwrap import dedent
            doc = dedent(doc).split('\n')
            doc[1:1] = depdoc
        newfunc.__doc__ =  '\n'.join(doc)
        try:
            d = func.__dict__
        except AttributeError:
            pass
        else:
            newfunc.__dict__.update(d)
        return newfunc

def deprecate(*args, **kwargs):
    """Issues a generic `DeprecationWarning`.

    This function may also be used as a decorator.

    Parameters
    ----------
    func : function
        The function to be deprecated.
    old_name : str, optional
        The name of the function to be deprecated. Default is None, in which
        case the name of `func` is used.
    new_name : str, optional
        The new name for the function. Default is None, in which case
        the deprecation message is that `old_name` is deprecated. If given,
        the deprecation message is that `old_name` is deprecated and `new_name`
        should be used instead.
    message : str, optional
        Additional explanation of the deprecation.  Displayed in the docstring
        after the warning.

    Returns
    -------
    old_func : function
        The deprecated function.
    """
    if args:
        fn = args[0]
        args = args[1:]

        return _Deprecate(*args, **kwargs)(fn)
    else:
        return _Deprecate(*args, **kwargs)

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
