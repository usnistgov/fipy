from __future__ import unicode_literals

__docformat__ = 'restructuredtext'

from package.subpackage.base import Base

class Object(Base):
    def __init__(self, arg1, arg2=None, arg3='string'):
        """
        This method, like all those whose names begin and end with
        "``__``" are special.  You won't ever need to call these
        methods directly, but :term:`Python` will invoke them for you under
        certain circumstances, which are described in the
        `Python Reference Manual: Special Method Names`_.

        As an example, the :meth:`__init__` method is invoked when you create an object, as in::

            obj = Object(arg1=something, arg3=somethingElse, ...)

        .. _`Python Reference Manual: Special Method Names`: https://docs.python.org/3/reference/datamodel.html#special-method-names

        Parameters
        ----------
        arg1
            this argument is required. :term:`Python` supports named arguments,
            so you must either list the value for `arg1`  first::

                obj = Object(val1, val2)

            or you can specify the arguments in any order, as long as they are named::

                obj = Object(arg2=val2, arg1=val1)

        arg2
            this argument may be omitted, in which case it will be
            assigned a default value of ``None``.  If you do not use named
            arguments (and we recommend that you do), all required
            arguments must be specified before any optional arguments.
        arg3
            this argument may be omitted, in which case it will be
            assigned a default value of ``'string'``.
        """
        pass

    def method2(self):
        """
        :class:`Object` provides a new definition for the behavior of
        :meth:`method2`, whereas the behavior of
        :meth:`~package.subpackage.base.Base.method1` is defined by
        :class:`~package.subpackage.base.Base`.
        """
        pass
