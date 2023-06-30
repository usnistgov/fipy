
"""
This module can be found in the file :file:`package/subpackage/base.py`.  You make it
available to your script by either::

    import package.subpackage.base

in which case you refer to it by its full name of :mod:`package.subpackage.base`, or::

    from package.subpackage import base

in which case you can refer simply to :mod:`~package.subpackage.base`.
"""
from __future__ import unicode_literals
from builtins import object
__docformat__ = 'restructuredtext'

class Base(object):
    """
    With very few exceptions, the name of a class will be the capitalized
    form of the module it resides in.  Depending on how you imported the
    module above, you will refer to either :class:`package.subpackage.object.Object` or
    :class:`object.Object`.  Alternatively, you can use::

        from package.subpackage.object import Object

    and then refer simply to :class:`~package.subpackage.object.Object`. For
    many classes, there is a shorthand notation::

        from fipy import Object


    :term:`Python` is an object-oriented language and the FiPy framework
    is composed of objects or classes.  Knowledge of object-oriented
    programming (OOP) is not necessary to use either Python or
    FiPy, but a few concepts are useful.  OOP involves two main
    ideas:

    **encapsulation**
      an object binds data with actions or "methods".  In most cases, you will
      not work with an object's data directly; instead, you will set, retrieve,
      or manipulate the data using the object's methods.

      Methods are functions that are attached to objects and that have
      direct access to the data of those objects.  Rather than
      passing the object data as an argument to a function::

          fn(data, arg1, arg2, ...)

      you instruct an object to invoke an appropriate method::

          object.meth(arg1, arg2, ...)

      If you are unfamiliar with object-oriented practices, there probably seems
      little advantage in this reordering.  You will have to trust us that
      the latter is a much more powerful way to do things.

    **inheritance**
      specialized objects are derived or inherited from more general objects.
      Common behaviors or data are defined in base objects and specific
      behaviors or data are either added or modified in derived objects.
      Objects that declare the existence of certain methods, without actually
      defining what those methods do, are called "abstract".  These objects
      exist to define the behavior of a family of objects, but rely on their
      descendants to actually provide that behavior.

      Unlike many object-oriented languages, :term:`Python` does not prevent the
      creation of abstract objects, but we will include a notice like

      .. attention:: This class is abstract. Always create one of its subclasses.

      for abstract classes which should be used for documentation but
      never actually created in a :term:`FiPy` script.
    """
    def __init__(self):
        pass

    def method1(self):
        """
        This is one thing that you can instruct any object that derives from
        :class:`~package.subpackage.base.Base` to do,
        by calling ``myObjectDerivedFromBase.``:meth:`method1`

        Parameters
        ----------
        self : object
            This special argument refers to the object that is being created.

            .. attention::

               *self* is supplied automatically by the :term:`Python` interpreter to all
               methods.  You don't need to (and should not) specify it yourself.
        """
        pass

    def method2(self):
        """
        This is another thing that you can instruct any object that derives from
        :class:`~package.subpackage.base.Base` to do.
        """
        pass
