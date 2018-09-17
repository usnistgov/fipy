## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 #
 #  FILE: "object.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #
 # ========================================================================
 # This software was developed by employees of the National Institute
 # of Standards and Technology, an agency of the Federal Government.
 # Pursuant to title 17 section 105 of the United States Code,
 # works of NIST employees are not subject to copyright
 # protection, and this software is considered to be in the public domain.
 # FiPy is an experimental system.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
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
 #  Description:
 #
 #    This file is purely for illustration of documentation conventions
 #
 # ###################################################################
 ##

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

        :Parameters:
          - `arg1`: this argument is required. :term:`Python` supports named arguments,
            so you must either list the value for `arg1`  first::

                obj = Object(val1, val2)

            or you can specify the arguments in any order, as long as they are named::

                obj = Object(arg2=val2, arg1=val1)

          - `arg2`: this argument may be omitted, in which case it will be
            assigned a default value of ``None``.  If you do not use named
            arguments (and we recommend that you do), all required
            arguments must be specified before any optional arguments.
          - `arg3`: this argument may be omitted, in which case it will be
            assigned a default value of ``'string'``.

        .. _`Python Reference Manual: Special Method Names`: http://www.python.org/doc/ref/specialnames.html
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
