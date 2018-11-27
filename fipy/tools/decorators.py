import re
import sys
import warnings

__all__ = ["deprecate"]

# Stolen from `numpy.lib.utils`
if sys.version_info < (2, 4):
    # Can't set __name__ in 2.3
    import new
    def _set_function_name(func, name):
        func = new.function(func.func_code, func.func_globals,
                            name, func.func_defaults, func.func_closure)
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
                old_name = func.func_name
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
