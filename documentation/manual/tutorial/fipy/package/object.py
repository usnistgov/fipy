## -*-Pyth-*-
 # ###################################################################
 #  FiPy - a finite volume PDE solver in Python
 # 
 #  FILE: "object.py"
 #                                    created: 8/16/05 {10:41:26 AM} 
 #                                last update: 8/18/05 {9:38:05 AM} 
 #  Author: Jonathan Guyer <guyer@nist.gov>
 #  Author: Daniel Wheeler <daniel.wheeler@nist.gov>
 #  Author: James Warren   <jwarren@nist.gov>
 #    mail: NIST
 #     www: http://www.ctcms.nist.gov/fipy/
 #  
 # ========================================================================
 # This document was prepared at the National Institute of Standards
 # and Technology by employees of the Federal Government in the course
 # of their official duties.  Pursuant to title 17 Section 105 of the
 # United States Code this document is not subject to copyright
 # protection and is in the public domain.  object.py
 # is an experimental work.  NIST assumes no responsibility whatsoever
 # for its use by other parties, and makes no guarantees, expressed
 # or implied, about its quality, reliability, or any other characteristic.
 # We would appreciate acknowledgement if the document is used.
 # 
 # This document can be redistributed and/or modified freely
 # provided that any derivative works bear some notice that they are
 # derived from it, and any modified versions bear some notice that
 # they have been modified.
 # ========================================================================
 #  See the file "license.terms" for information on usage and  redistribution
 #  of this file, and for a DISCLAIMER OF ALL WARRANTIES.
 #  
 #  Description: 
 #   
 #    This file is purely for illustration of documentation conventions
 #
 #  History
 # 
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  1939-08-16 JEG 1.0 original
 # ###################################################################
 ##

"""
This module can be found in the file `fipy/package/object.py`.  You make it
available to your script by either:
    
    import fipy.package.object
    
in which case you refer to it by its full name of `fipy.package.object`, or:
    
    from fipy.package import object
    
in which case you can refer simply to `object`.
"""
__docformat__ = 'restructuredtext'

from fipy.package.base import Base

class Object(Base):
    """
    With very few exceptions, the name of a class will be the capitalized
    form of the module it resides in.  Depending on how you imported the
    module above, you will refer to either `fipy.package.object.Object` or
    `object.Object`.  Alternatively, you can use:
        
        from fipy.package.object import Object
        
    and then refer simply to `Object`. There is a shorthand notation::
        
        from fipy import Object
        
    but it is still experimental and does not work for all of the objects
    in FiPy.
    
    Python_ is an object-oriented language and the FiPy framework
    is composed of objects or classes.  Knowledge of object-oriented
    programming (OOP) is not necessary to use either Python or
    FiPy, but a few concepts are useful.  OOP involves two main
    ideas

    **encapsulation** 
      an object binds data with actions or 
      "methods". In most cases, you will not work with an 
      object's data directly; instead, you will set, retrieve, or 
      manipulate the data using the object's methods.

      Methods are functions that are attached to objects and that have
      direct access to the data of those objects.  Rather than 
      passing the object data as an argument to a function::
          
          fn(data, arg1, arg2, ...)
          
      you instruct an object to invoke an appropriate method::
          
          object.meth(arg1, arg2, ...)
          
    **inheritance** 
      objects are derived or inherited from more 
      abstract objects. Common behaviors or data are defined in 
      base objects and specific behaviors or data are either added 
      or modified in derived objects.
      
    .. _Python:               http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://www.python.org/

    """
    def __init__(self, arg1, arg2 = None, arg3 = 'string'):
        """
        This method, like all those whose names begin and end with
        "`__`" are special.  You won't ever need to call these
        methods directly, but Python_ will invoke them for you under
        certain circumstances, which are described in the 
        `Python Reference Manual: Special Method Names`_ |citePythonSpecialMethods|.

        As an example, the `__init__` method is invoked when you create an object, as in::
        
            obj = Object(arg1 = something, arg3 = somethingElse, ...)
            
        :Parameters:
          - `self`: this special argument
            refers to the object that is being created and is
            supplied automatically by the Python interpreter
            to all methods.  You don't need to (and should
            not) specify it yourself.
          - `arg1`: this argument is required. Python_ supports named arguments, 
            so you must either list the value for *arg1*  first::
                      
                obj = Object(val1, val2)
               
            or you can specify the arguments in any order, as long as they are named::
            
                obj = Object(arg2 = val2, arg1 = val1)
                
          - `arg2`: this argument may be omitted, in which case it will be 
            assigned a default value of `None`.  If you do not use named
            arguments (and we recommend that you do), all required
            arguments must be specified before any optional arguments.
          - `arg3`: this argument may be omitted, in which case it will be
            assigned a default value of `'string'`.
        
        .. _Python:               http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://www.python.org/
        .. _`Python Reference Manual: Special Method Names`: http://www.nist.gov/cgi-bin/exit_nist.cgi?url=http://www.python.org/doc/ref/specialnames.html
        .. |citePythonSpecialMethods| raw:: latex

           \cite[\S 3.3]{PythonReference}

        """
        pass

##     def __abs__(self):
##         """
##         This method is invoked whenever you take the absolute value of an
##         `Object` by applying the `abs()` function.
##         """
##         pass
## 
##     def __add__(self, other):
##         """
##         This method is invoked whenever you add `other` to an `Object`::
##             
##             a = Object(...) + other
##         """
##         pass
##         
##     def __radd__(self, other):
##         """
##         This method is invoked whenever you add an `Object` to `other`::
##             
##             a = other + Object(...)
##         """
##         pass
##         
##     def __and__(self, other):
##         """
##         This method is invoked whenever you logically and an `Object` with `other`::
##             
##             a = Object(...) and other
##         """
##         pass
##         
##     def __array__(self):
##         """ 
##         This method is invoked whenever you attempt to coerce an `Object` to a Numeric array::
##             
##             a = Numeric.array(Object(...))
##         """
##         pass
##         
##     def __call__(self, arg):
##         """
##         This method is invoked whenever you attempt to call an `Object` like a function::
##             
##             o = Object(...)
##             a = o(arg = val)
##         """
##         pass
##         
##     def __div__(self, other):
##         """
##         This method is invoked whenever you divide an `Object` by `other`::
##             
##             a = Object(...) / other
##         """
##         pass
## 
##     def __rdiv__(self, other):
##         """
##         This method is invoked whenever you divide `other` by an `Object`::
##             
##             a = other / Object(...)
##         """
##         pass
## 
##     def __eq__(self, other):
##         """
##         This method is invoked whenever you test if an `Object` is equal to
##         `other`::
##             
##             a = (Object(...) == other)
##         """
##         pass
##         
##     def __float__(self):
##         """
##         This method is invoked whenever you convert an `Object` to a
##         floating point number by applying the `float()` function.
##         """
##         pass
##         
##     def __ge__(self, other):
##         """
##         This method is invoked whenever you test if an `Object` is greater than or equal to
##         `other`::
##             
##             a = (Object(...) >= other)
##         """
##         pass
##         
##     def __getitem__(self, index):
##         """
##         This method is invoked whenever you treat an `Object` as an array
##         and "slice" it by requesting the contents at a particular
##         index (or indices)::
##             
##             a = Object(...)[index]
##         """
##         pass
##         
##     def __getstate__(self):
##         """
##         This method is invoked whenever you "pickle" an `Object` to
##         persistant storage.
##         """
##         pass
##         
##     def __gt__(self, other):
##         """
##         This method is invoked whenever you test if an `Object` is greater than
##         `other`::
##             
##             a = (Object(...) > other)
##         """
##         pass
## 
##     def __le__(self, other):
##         """
##         This method is invoked whenever you test if an `Object` is less than
##         or equal to `other`::
##             
##             a = (Object(...) <= other)
##         """
##         pass
## 
##     def __len__(self, index):
##         """
##         This method is invoked whenever you treat an `Object` as a list
##         and requesting its length with the `len()` function.
##         """
##         pass
## 
## 
##     def __lt__(self, other):
##         """
##         This method is invoked whenever you test if an `Object` is less than
##         `other`::
##             
##             a = (Object(...) < other)
##         """
##         pass
## 
##     def __mod__(self, other):
##         """
##         This method is invoked whenever you request the integer remainder
##         of dividing an `Object` by `other`::
##             
##             a = Object(...) % other
##         """
##         pass
## 
##     def __rmod__(self, other):
##         """
##         This method is invoked whenever you request the integer remainder
##         of dividing `other` by an `Object`::
##             
##             a = other % Object(...)
##         """
##         pass
## 
##     def __mul__(self, other):
##         """
##         This method is invoked whenever you multiply an `Object` by
##         `other`::
##             
##             a = Object(...) * other
##         """
##         pass
## 
##     def __rmul__(self, other):
##         """
##         This method is invoked whenever you multiply a `other` by an
##         `Object`::
##             
##             a = other * Object(...)
##         """
##         pass
## 
##     def __ne__(self, other):
##         """
##         This method is invoked whenever you test if an `Object` is not
##         equal to `other`::
##             
##             a = (Object(...) != other)
##         """
##         pass
##         
##     def __neg__(self):
##         """
##         This method is invoked whenever you negate an `Object`::
##             
##             a = -Object(...)
##         """
##         pass
##         
##     def __pos__(self):
##         """
##         This method is invoked whenever you posivate(?) an `Object`::
##             
##             a = +Object(...)
##         """
##         pass
##         
##     def __pow__(self, other):
##         """
##         This method is invoked whenever you exponentiate `Object` to a power
##         `other`::
##             
##             a = Object(...) ** other
##         """
##         pass
## 
##     def __rpow__(self, other):
##         """
##         This method is invoked whenever you exponentiate `other` to a power
##         of an `Object`::
##             
##             a = other ** Object(...)
##         """
##         pass
## 
##     def __repr__(self):
##         """
##         This method is invoked whenever you request a printable
##         representation of an `Object` with the `repr()` function (see also
##         `__str__`), which should be an acceptable argument to `eval()`.
##         """
##         pass
## 
##     def __str__(self):
##         """
##         This method is invoked whenever you request a nicely printable
##         representation of an `Object` with the `str()` function (see also
##         `__repr__`).
##         """
##         pass
## 
##     def __sub__(self, other):
##         """
##         This method is invoked whenever you subtract `other` from an
##         `Object`::
##             
##             a = Object(...) - other
##         """
##         pass
## 
##     def __rsub__(self, other):
##         """
##         This method is invoked whenever you subtract an `Object` from
##         `other`::
##             
##             a = other - Object(...)
##         """
##         pass
## 
##     def __setitem__(self, index, value):
##         """
##         This method is invoked whenever you treat an `Object` as an array
##         and assign `value` to a "slice" at a particular index (or
##         indices)::
##             
##             Object(...)[index] = value
##         """
##         pass
## 
##     def __setstate__(self, dict):
##         """
##         This method is invoked whenever you create an `Object` from
##         "pickled" persistant storage.
##         """
##         pass
        
