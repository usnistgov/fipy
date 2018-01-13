#!/usr/bin/env python

## -*-Pyth-*-
 # ###################################################################
 #  FiPy - Python-based finite volume PDE solver
 #
 #  FILE: "physicalField.py"
 #
 #  Author: Jonathan Guyer <guyer@nist.gov>
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
 #  Copyright 1997-2006 by Konrad Hinsen, except as noted below.
 #
 #  Permission to use, copy, modify, and distribute this software and its
 #  documentation for any purpose and without fee is hereby granted,
 #  provided that the above copyright notice appear in all copies and that
 #  both that copyright notice and this permission notice appear in
 #  supporting documentation.
 #
 #  THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE,
 #  INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO
 #  EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, INDIRECT OR
 #  CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF
 #  USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 #  OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 #  PERFORMANCE OF THIS SOFTWARE.
 #
 #  Description:
 #
 # Physical fields or quantities with units
 #
 # Based on PhysicalQuantities of the Scientific package, written by Konrad
 # Hinsen <hinsen@cnrs-orleans.fr>
 #
 #  History
 #
 #  modified   by  rev reason
 #  ---------- --- --- -----------
 #  2003-12-28 JEG 1.0 support for dimensional fields
 #  1998/09/29 GPW     now supports conversions with offset
 #                     (for temperature units)
 #  1998/09/28 GPW     now removes __args__ from local dict after eval
 # ###################################################################
 ##


"""
Physical quantities with units.

This module derives from `Konrad Hinsen`_'s :class:`PhysicalQuantity`
<http://dirac.cnrs-orleans.fr/ScientificPython/ScientificPythonManual/Scientific.Physics.PhysicalQuantities-module.html>.

This module provides a data type that represents a physical
quantity together with its unit. It is possible to add and
subtract these quantities if the units are compatible, and
a quantity can be converted to another compatible unit.
Multiplication, subtraction, and raising to integer powers
is allowed without restriction, and the result will have
the correct unit. A quantity can be raised to a non-integer
power only if the result can be represented by integer powers
of the base units.

The values of physical constants are taken from the 2002
recommended values from CODATA_. Other conversion factors
(e.g. for British units) come from `Appendix B of NIST Special Publication 811`_.

.. warning::

   We can't guarantee for the correctness of all entries in the unit table,
   so use this at your own risk!

.. _Konrad Hinsen:                              mailto:hinsen@cnrs-orleans.fr
.. _PhysicalQuantity:                           http://starship.python.net/~hinsen/ScientificPython/ScientificPythonManual/Scientific_31.html
.. _CODATA:                                     http://www.codata.org/
.. _Appendix B of NIST Special Publication 811: http://physics.nist.gov/Pubs/SP811/appenB9.html
"""
__docformat__ = 'restructuredtext'

import re

from fipy.tools import numerix
from fipy.tools.numerix import MA
from fipy.tools.numerix import umath

from fipy.tools.dimensions.NumberDict import _NumberDict

__all__ = ["PhysicalField", "PhysicalUnit"]

# Class definitions

class PhysicalField(object):
    """
    Physical field or quantity with units
    """


    def __init__(self, value, unit = None, array = None):
        """
        Physical Fields can be constructed in one of two ways:

          - `PhysicalField(*value*, *unit*)`, where `*value*` is a number of
            arbitrary type and `*unit*` is a string containing the unit name

                >>> print PhysicalField(value = 10., unit = 'm')
                10.0 m

          - `PhysicalField(*string*)`, where `*string*` contains both the value
            and the unit. This form is provided to make interactive use more
            convenient

                >>> print PhysicalField(value = "10. m")
                10.0 m


        Dimensionless quantities, with a `unit` of 1, can be specified
        in several ways

            >>> print PhysicalField(value = "1")
            1.0 1
            >>> print PhysicalField(value = 2., unit = " ")
            2.0 1
            >>> print PhysicalField(value = 2.)
            2.0 1

        Physical arrays are also possible (and are the reason this code
        was adapted from `Konrad Hinsen`_'s original PhysicalQuantity_).
        The `value` can be a Numeric_ `array`:

            >>> a = numerix.array(((3.,4.),(5.,6.)))
            >>> print PhysicalField(value = a, unit = "m")
            [[ 3.  4.]
             [ 5.  6.]] m

        or a `tuple`:

            >>> print PhysicalField(value = ((3.,4.),(5.,6.)), unit = "m")
            [[ 3.  4.]
             [ 5.  6.]] m

        or as a single value to be applied to every element of a supplied array:

            >>> print PhysicalField(value = 2., unit = "m", array = a)
            [[ 2.  2.]
             [ 2.  2.]] m

        Every element in an array has the same unit, which is stored only
        once for the whole array.

        .. _Konrad Hinsen: mailto:hinsen@cnrs-orleans.fr
        .. _PhysicalQuantity: http://starship.python.net/~hinsen/ScientificPython/ScientificPythonManual/Scientific_31.html
        .. _Numeric: http://www.numpy.org
        """
        if hasattr(value, "value") and not isinstance(value, PhysicalField):
            value = value.value
        elif hasattr(value, "getValue") and callable(value.getValue):
            value = value.getValue()

        if isinstance(value, PhysicalField):
            if unit is not None:
                if not value.unit.isDimensionless():
                    value = value.inUnitsOf(unit)
                unit = _findUnit(unit)
            else:
                unit = value.unit

            if hasattr(value.value, 'copy'):
                value = value.value.copy()
            else:
                value = value.value
        elif unit is not None:
            unit = _findUnit(unit)
        elif type(value) is type(''):
            s = value.strip()
            match = PhysicalField._number.match(s)
            if match is None:
                value = 1
                unit = _findUnit(s)
#                   raise TypeError, 'No number found'
            else:
                value = float(match.group(0))
                unit = _findUnit(s[len(match.group(0)):])

        if type(value) in [type([]),type(())]:
            value = [PhysicalField(item,unit) for item in value]
            if unit is None:
                unit = value[0].unit
            normalized = []
            for item in value:
                if item.unit == unit:
                    normalized += [item.value]
                else:
                    normalized += [item.inUnitsOf(unit).value]
            value = numerix.array(normalized)

        if unit is None:
            unit = _unity
##             unit = _findUnit("")

        self.value = numerix.array(value)
        self._unit = unit
        if array is not None:
            array[:] = self.value
            self.value = array

    _number = re.compile('[+-]?[0-9]+(\\.[0-9]*)?([eE][+-]?[0-9]+)?')

    def ravel(self):
        return self.value.ravel()

    def copy(self):
        """
        Make a duplicate.

            >>> a = PhysicalField(1, unit = 'inch')
            >>> b = a.copy()

        The duplicate will not reflect changes made to the original

            >>> a.convertToUnit('cm')
            >>> print a
            2.54 cm
            >>> print b
            1 inch

        Likewise for arrays

            >>> a = PhysicalField(numerix.array((0,1,2)), unit  = 'm')
            >>> b = a.copy()
            >>> a[0] = 3
            >>> print a
            [3 1 2] m
            >>> print b
            [0 1 2] m

        """
        if hasattr(self.value, 'copy'):
            return PhysicalField(value = self.value.copy(), unit = self.unit)
        else:
            return PhysicalField(value = self.value, unit = self.unit)

    def __str__(self):
        """
        Return human-readable form of a physical quantity

            >>> print PhysicalField(value = 3., unit = "eV")
            3.0 eV
        """
        return str(self.value) + ' ' + self.unit.name()

    def __repr__(self):
        """
        Return representation of a physical quantity suitable for re-use

            >>> PhysicalField(value = 3., unit = "eV")
            PhysicalField(3.0,'eV')
        """
        value = self.value
        if value.shape == ():
            value = str(value)
        else:
            value = repr(value)

        if self.unit is _unity:
            return value
        else:
            return (self.__class__.__name__ + '(' + value + ',' +
                    repr(self.unit.name()) + ')')

    def tostring(self, max_line_width=75, precision=8, suppress_small=False, separator=' '):
        """
        Return human-readable form of a physical quantity

            >>> p = PhysicalField(value = (3., 3.14159), unit = "eV")
            >>> print p.tostring(precision = 3, separator = '|')
            [ 3.   | 3.142] eV
        """
        from fipy.tools import numerix
        return numerix.tostring(self.value, max_line_width = max_line_width,
                                precision = precision,
                                suppress_small = suppress_small,
                                separator = separator) + ' ' + self.unit.name()


    def _sum(self, other, sign1 = lambda a: a, sign2 = lambda b: b):
        selfValue = self.value

        if _isVariable(other):
            return sign2(other) + self.__class__(value = sign1(selfValue), unit = self.unit)

        if type(other) is type(''):
            other = PhysicalField(value = other)

        if not isinstance(other,PhysicalField):
            if numerix.alltrue(other == 0):
                new_value = sign1(selfValue)
            elif self.unit.isDimensionlessOrAngle() or self.unit.isInverseAngle():
                new_value = sign1(selfValue) + sign2(other)
            else:
                raise TypeError, str(self) + ' and ' + str(other) + ' are incompatible.'
        else:
            new_value = (sign1(selfValue)
                         + sign2(other.value) * other.unit.conversionFactorTo(self.unit))
        return self.__class__(value=new_value, unit=self.unit)

    def __add__(self, other):
        """
        Add two physical quantities, so long as their units are compatible.
        The unit of the result is the unit of the first operand.

            >>> print PhysicalField(10., 'km') + PhysicalField(10., 'm')
            10.01 km
            >>> print PhysicalField(10., 'km') + PhysicalField(10., 'J')
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return self._sum(other)

    __radd__ = __add__
    add = __add__

    def __sub__(self, other):
        """
        Subtract two physical quantities, so long as their units are compatible.
        The unit of the result is the unit of the first operand.

            >>> print PhysicalField(10., 'km') - PhysicalField(10., 'm')
            9.99 km
            >>> print PhysicalField(10., 'km') - PhysicalField(10., 'J')
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return self._sum(other, sign2 = lambda b: -b)

    subtract = __sub__

    def __rsub__(self, other):
        return self._sum(other, sign1 = lambda a: -a)

    def __mul__(self, other):
        """
        Multiply two physical quantities.  The unit of the result is the
        product of the units of the operands.

            >>> print PhysicalField(10., 'N') * PhysicalField(10., 'm') == PhysicalField(100., 'N*m')
            True

        As a special case, if the result is dimensionless, the value
        is returned without units, rather than with a dimensionless unit
        of `1`. This facilitates passing physical quantities to packages
        such as Numeric that cannot use units, while ensuring the quantities
        have the desired units.

            >>> print (PhysicalField(10., 's') * PhysicalField(2., 'Hz'))
            20.0
        """
        if _isVariable(other):
            return other.__mul__(self)
        if type(other) is type(''):
            other = PhysicalField(value = other)
        if not isinstance(other,PhysicalField):
            return self.__class__(value = self.value*other, unit = self.unit)
        value = self.value*other.value
        unit = self.unit*other.unit
        if unit.isDimensionless():
            if unit.factor != 1:
                return value * unit.factor
            else:
                return value
        else:
            return self.__class__(value = value, unit = unit)

    __rmul__ = __mul__
    multiply = __mul__

    def __truediv__(self, other):
        """
        Divide two physical quantities.  The unit of the result is the
        unit of the first operand divided by the unit of the second.

            >>> print PhysicalField(10., 'm') / PhysicalField(2., 's')
            5.0 m/s

        As a special case, if the result is dimensionless, the value
        is returned without units, rather than with a dimensionless unit
        of `1`. This facilitates passing physical quantities to packages
        such as Numeric_ that cannot use units, while ensuring the quantities
        have the desired units

            >>> print (PhysicalField(1., 'inch')
            ...        / PhysicalField(1., 'mm'))
            25.4

        .. _Numeric: http://www.numpy.org
        """
        if _isVariable(other):
            return other.__rtruediv__(self)
        if type(other) is type(''):
            other = self.__class__(value = other)
        if not isinstance(other,PhysicalField):
            value = self.value/other
            unit = self.unit
        else:
            value = self.value/other.value
            unit = self.unit/other.unit
        if unit.isDimensionless():
            return value*unit.factor
        else:
            return self.__class__(value = value, unit = unit)

    __div__ = __truediv__
    divide = __truediv__

    def __rtruediv__(self, other):
        if _isVariable(other):
            return other.__truediv__(self)
        if type(other) is type(''):
            other = PhysicalField(value = other)
        if not isinstance(other,PhysicalField):
            value = other/self.value
            unit = pow(self.unit, -1)
        else:
            value = other.value/self.value
            unit = other.unit/self.unit
        if unit.isDimensionless():
            return value*unit.factor
        else:
            return self.__class__(value = value, unit = unit)

    __rdiv__ = __rtruediv__

    def __mod__(self, other):
        """
        Return the remainder of dividing two physical quantities.  The unit of the result is the
        unit of the first operand divided by the unit of the second.

            >>> print PhysicalField(11., 'm') % PhysicalField(2., 's')
            1.0 m/s
        """
        if _isVariable(other):
            return other.__rmod__(self)
        if type(other) is type(''):
            other = self.__class__(value = other)
        if not isinstance(other,PhysicalField):
            value = self.value % other
            unit = self.unit
        else:
            value = self.value % other.value
            unit = self.unit/other.unit
        if unit.isDimensionless():
            return value*unit.factor
        else:
            return self.__class__(value = value, unit = unit)

    def __pow__(self, other):
        """
        Raise a `PhysicalField` to a power. The unit is raised to the same power.

            >>> print PhysicalField(10., 'm')**2
            100.0 m**2
        """
        if type(other) is type(''):
            other = PhysicalField(value = other)
        return self.__class__(value = pow(self.value, float(other)), unit = pow(self.unit, other))

    def __rpow__(self, other):
        return pow(other, float(self))

    def __abs__(self):
        """
        Return the absolute value of the quantity. The `unit` is unchanged.

            >>> print abs(PhysicalField(((3.,-2.),(-1.,4.)), 'm'))
            [[ 3.  2.]
             [ 1.  4.]] m
        """
        return self.__class__(value = abs(self.value), unit = self.unit)

    def __pos__(self):
        return self

    def __neg__(self):
        """
        Return the negative of the quantity. The `unit` is unchanged.

            >>> print -PhysicalField(((3.,-2.),(-1.,4.)), 'm')
            [[-3.  2.]
             [ 1. -4.]] m
        """
        return self.__class__(value = -self.value, unit = self.unit)

    def sign(self):
        """
        Return the sign of the quantity. The `unit` is unchanged.

            >>> from fipy.tools.numerix import sign
            >>> print sign(PhysicalField(((3.,-2.),(-1.,4.)), 'm'))
            [[ 1. -1.]
             [-1.  1.]]
        """
        return numerix.sign(self.value)

    def __nonzero__(self):
        """
        Test if the quantity is zero.

        Should this only pass if the unit offset is zero?
        """
        return self.value != 0

    def _inMyUnits(self, other):
        if _isVariable(other):
            other = other.value

        if not isinstance(other, PhysicalField):
            if type(other) is type(''):
                other = PhysicalField(other)
            elif numerix.alltrue(other == 0) or self.unit.isDimensionlessOrAngle():
                other = PhysicalField(value = other, unit = self.unit)
            else:
                raise TypeError, 'Incompatible units'
        return other.inUnitsOf(self.unit)

    def __getitem__(self, index):
        """
        Return the specified element of the array. The unit of the result
        will be the unit of the array.

            >>> a = PhysicalField(((3.,4.),(5.,6.)),"m")
            >>> print a[1,1]
            6.0 m
        """
        return PhysicalField(self.value[index],self.unit)

    def __setitem__(self, index, value):
        """
        Assign the specified element of the array, performing appropriate
        conversions.

            >>> a = PhysicalField(((3.,4.),(5.,6.)),"m")
            >>> a[0,1] = PhysicalField("6 ft")
            >>> print a
            [[ 3.      1.8288]
             [ 5.      6.    ]] m
            >>> a[1,0] = PhysicalField("2 min")
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        if type(value) is type(''):
            value = PhysicalField(value)
        if isinstance(value, PhysicalField) or _isVariable(value):
            value = self._inMyUnits(value).value
        self.value[index] = value

    def itemset(self, value):
        """
        Assign the value of a scalar array, performing appropriate
        conversions.

            >>> a = PhysicalField(4.,"m")
            >>> a.itemset(PhysicalField("6 ft"))
            >>> print a.allclose("1.8288 m")
            1
            >>> a = PhysicalField(((3.,4.),(5.,6.)),"m")
            >>> try:
            ...     a.itemset(PhysicalField("6 ft"))
            ... except IndexError:
            ...     # NumPy 1.7 has changed the exception type
            ...     raise ValueError("can only place a scalar for an  array of size 1")
            Traceback (most recent call last):
                ...
            ValueError: can only place a scalar for an  array of size 1
            >>> a.itemset(PhysicalField("2 min"))
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        if type(value) is type(''):
            value = PhysicalField(value)
        if isinstance(value, PhysicalField) or _isVariable(value):
            value = self._inMyUnits(value).value
        self.value.itemset(value)


##    __array_priority__ and __array_wrap__ are required to override
##    the default behavior of numpy. If a numpy array and a Variable
##    are in a binary operation and numpy is first, then numpy will,
##    by default, try and do everything it can to get a a raw numpy
##    array out of Variable. __array_wrap__ seems to have been
##    introduced into masked array to fix this issue. __array_wrap__ is
##    called after the operation is done so it could hurt efficiency badly.
##    Something else needs to be done to stop the initial evaluation.

    __array_priority__ = 100.0

    def __array_wrap__(self, arr, context=None):
        """
        Required to prevent numpy not calling the reverse binary operations.
        Both the following tests are examples ufuncs.

        >>> from fipy.tools.dimensions.physicalField import PhysicalField
        >>> print type(numerix.array([1.0, 2.0]) * PhysicalField([1.0, 2.0], unit="m"))
        <class 'fipy.tools.dimensions.physicalField.PhysicalField'>

        >>> print type(numerix.array([1.0, 2.0]) * PhysicalField([1.0, 2.0]))
        <class 'fipy.tools.dimensions.physicalField.PhysicalField'>

        >>> from scipy.special import gamma as Gamma # doctest: +SCIPY
        >>> print type(Gamma(PhysicalField([1.0, 2.0]))) is type(numerix.array(1)) # doctest: +SCIPY
        1
        """
        result = arr

        if context is not None:
            (func, args, _) = context
            def __makePhysical(v):
                if not (isinstance(v, PhysicalField) or _isVariable(v)):
                    v = PhysicalField(v)
                return v
            args = [__makePhysical(arg) for arg in args]

            meth = getattr(args[0], func.__name__, None)
            if meth is not None and callable(meth):
                result = meth(*args[1:])

        return result

    def __array__(self, t = None):
        """
        Return a dimensionless `PhysicalField` as a Numeric_ ``array``.

            >>> print numerix.array(PhysicalField(((2.,3.),(4.,5.)),"m/m"))
            [[ 2.  3.]
             [ 4.  5.]]

        As a special case, fields with angular units are converted to base
        units (radians) and then assumed dimensionless.

            >>> print numerix.array(PhysicalField(((2.,3.),(4.,5.)),"deg"))
            [[ 0.03490659  0.05235988]
             [ 0.06981317  0.08726646]]


        If the array is not dimensionless, the numerical value in its base
        units is returned.

            >>> numerix.array(PhysicalField(((2.,3.),(4.,5.)),"mm"))
            array([[ 0.002,  0.003],
                   [ 0.004,  0.005]])

        .. _Numeric: http://www.numpy.org
        """
        if self.unit.isAngle():
            value = self.inRadians()
        else:
            value = self.inSIUnits().value

        return numerix.array(value, t)

#         if self.unit.isDimensionlessOrAngle():
#             value = self.numericValue
# ##            if type(value) is type(numerix.array((0))) and (t is None or t == value.typecode()):
#             if type(value) is type(numerix.array((0))) and (t is None or t == value.dtype.char):
#                 return value
#             else:
#                 return numerix.array(self.numericValue, t)
#         else:
#             raise TypeError, 'Numeric array value must be dimensionless'

    @property
    def _array(self):
        if self.unit.isDimensionlessOrAngle():
            return self.value
        else:
            raise TypeError, 'Numeric array value must be dimensionless'

    def __float__(self):
        """
        Return a dimensionless PhysicalField quantity as a float.

            >>> float(PhysicalField("2. m/m"))
            2.0

        As a special case, quantities with angular units are converted to
        base units (radians) and then assumed dimensionless.

            >>> print numerix.round_(float(PhysicalField("2. deg")), 6)
            0.034907

        If the quantity is not dimensionless, the conversion fails.

            >>> float(PhysicalField("2. m"))
            Traceback (most recent call last):
                ...
            TypeError: Not possible to convert a PhysicalField with dimensions to float

        Just as a Numeric_ `array` cannot be cast to float, neither can PhysicalField arrays

            >>> float(PhysicalField(((2.,3.),(4.,5.)),"m/m")) # doctest: +ELLIPSIS
            Traceback (most recent call last):
                ...
            TypeError: only ...-1 arrays can be converted to Python scalars

        .. _Numeric: http://www.numpy.org
        """
        if self.unit.isDimensionlessOrAngle():
            return float(self.numericValue)
        else:
            raise TypeError, 'Not possible to convert a PhysicalField with dimensions to float'

    def __gt__(self,other):
        """
        Compare `self` to `other`, returning an array of boolean values
        corresponding to the test against each element.

            >>> a = PhysicalField(((3.,4.),(5.,6.)),"m")
            >>> print numerix.allclose(a > PhysicalField("13 ft"),
            ...                        [[False, True], [ True, True]])
            True

        Appropriately formatted dimensional quantity strings can also be
        compared.

            >>> print numerix.allclose(a > "13 ft",
            ...                        [[False, True], [ True, True]])
            True

        Arrays are compared element to element

            >>> print numerix.allclose(a > PhysicalField(((3.,13.),(17.,6.)),"ft"),
            ...                        [[ True, True], [False, True]])
            True

        Units must be compatible

            >>> print a > PhysicalField("1 lb")
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units

        And so must array dimensions

            >>> print a > PhysicalField(((3.,13.,4.),(17.,6.,2.)),"ft") #doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                ...
            ValueError: shape mismatch: objects cannot be broadcast to a single shape
        """
        other = self._inMyUnits(other)
        return self.value > other.value

    def __lt__(self,other):
        other = self._inMyUnits(other)
        return self.value < other.value

    def __le__(self,other):
        other = self._inMyUnits(other)
        return self.value <= other.value

    def __eq__(self,other):
        other = self._inMyUnits(other)
        return self.value == other.value

    __hash__ = object.__hash__

    def __ne__(self,other):
        other = self._inMyUnits(other)
        return self.value != other.value

    def __ge__(self,other):
        other = self._inMyUnits(other)
        return self.value >= other.value

    def __len__(self):
        if type(self.value) in [type(1), type(1.)]:
            return 1
        else:
            return len(self.value)

    def convertToUnit(self, unit):
        """
        Changes the unit to `unit` and adjusts the value such that the
        combination is equivalent.  The new unit is by a string containing
        its name.  The new unit must be compatible with the previous unit
        of the object.

            >>> e = PhysicalField('2.7 Hartree*Nav')
            >>> e.convertToUnit('kcal/mol')
            >>> print e
            1694.27557621 kcal/mol
        """
        unit = _findUnit(unit)
        self.value = _convertValue (self.value, self.unit, unit)
        self._unit = unit

    def inRadians(self):
        """
        Converts an angular quantity to radians and returns the numerical value.

            >>> print PhysicalField(((2.,3.),(4.,5.)),"rad").inRadians()
            [[ 2.  3.]
             [ 4.  5.]]
            >>> print PhysicalField(((2.,3.),(4.,5.)),"deg").inRadians()
            [[ 0.03490659  0.05235988]
             [ 0.06981317  0.08726646]]

         As a special case, assumes a dimensionless quantity is already in
         radians.

             >>> print PhysicalField(((2.,3.),(4.,5.))).inRadians()
             [[ 2.  3.]
              [ 4.  5.]]

         It's an error to convert a quantity with non-angular units

            >>> print PhysicalField(((2.,3.),(4.,5.)),"m").inRadians()
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units

        """
        if self.unit.isDimensionless():
            return self.value
        else:
            return self.inUnitsOf("rad").value

    def inDimensionless(self):
        """
        Returns the numerical value of a dimensionless quantity.

             >>> print PhysicalField(((2.,3.),(4.,5.))).inDimensionless()
             [[ 2.  3.]
              [ 4.  5.]]

         It's an error to convert a quantity with units

            >>> print PhysicalField(((2.,3.),(4.,5.)),"m").inDimensionless()
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units

        """
        if self.unit.isDimensionless():
            return self.value
        else:
            raise TypeError, 'Incompatible units'

    def inUnitsOf(self, *units):
        """
	Returns one or more `PhysicalField` objects that express the same
        physical quantity in different units.  The units are specified by
        strings containing their names.  The units must be compatible with
        the unit of the object.  If one unit is specified, the return value
        is a single `PhysicalField`.

        >>> freeze = PhysicalField('0 degC')
        >>> print freeze.inUnitsOf('degF').allclose("32.0 degF")
        1

        If several units are specified, the return value is a tuple of
	`PhysicalField` instances with with one element per unit such that
        the sum of all quantities in the tuple equals the the original
        quantity and all the values except for the last one are integers.
        This is used to convert to irregular unit systems like
        hour/minute/second.  The original object will not be changed.

        >>> t = PhysicalField(314159., 's')
        >>> print numerix.allclose([e.allclose(v) for (e, v) in zip(t.inUnitsOf('d','h','min','s'),
        ...                                                         ['3.0 d', '15.0 h', '15.0 min', '59.0 s'])],
        ...                        True)
        1
        """
        units = map(_findUnit, units)
        if len(units) == 1:
            unit = units[0]
            value = _convertValue (self.value, self.unit, unit)
            return self.__class__(value = value, unit = unit)
        else:
            units.sort()
            result = []
            value = self.value
            unit = self.unit
            for i in range(len(units)-1,-1,-1):
                value = value*unit.conversionFactorTo(units[i])
                if i == 0:
                    rounded = value
                else:
                    rounded = _round(value)
                result.append(self.__class__(value = rounded, unit = units[i]))
                value = value - rounded
                unit = units[i]
            return tuple(result)

    def getsctype(self, default=None):
        """

        Returns the Numpy sctype of the underlying array.

            >>> PhysicalField(1, 'm').getsctype() == numerix.NUMERIX.obj2sctype(numerix.array(1))
            True
            >>> PhysicalField(1., 'm').getsctype() == numerix.NUMERIX.obj2sctype(numerix.array(1.))
            True
            >>> PhysicalField((1,1.), 'm').getsctype() == numerix.NUMERIX.obj2sctype(numerix.array((1., 1.)))
            True

        """

        if not hasattr(self, 'typecode'):
            self.typecode = numerix.obj2sctype(rep=numerix.array(self.numericValue),
                                               default=default)

        return self.typecode

    def _getUnit(self):
        """
        Return the unit object of `self`.

            >>> PhysicalField("1 m").unit
            <PhysicalUnit m>
        """
        return self._unit

    def _setUnit(self, unit):
        """
        Change the unit object of `self` to `unit`

            >>> a = PhysicalField(value="1 m")
            >>> a.unit = ("m**2/s")
            >>> print a
            1.0 m**2/s
        """
        self._unit = _findUnit(unit)

    unit = property(_getUnit, _setUnit)

    @property
    def numericValue(self):
        """
        Return the `PhysicalField` without units, after conversion to base SI units.

            >>> print numerix.round_(PhysicalField("1 inch").numericValue, 6)
            0.0254
        """
        return self.inSIUnits().value

    # Contributed by Berthold Hoellmann
    def inBaseUnits(self):
        """
        Return the quantity with all units reduced to their base SI elements.

            >>> e = PhysicalField('2.7 Hartree*Nav')
            >>> print e.inBaseUnits().allclose("7088849.01085 kg*m**2/s**2/mol")
            1
        """
        if self.unit.factor != 1:
            new_value = self.value * self.unit.factor
        else:
            new_value = self.value
        num = ''
        denom = ''
        for i in xrange(9):
            unit = _base_names[i]
            power = self.unit.powers[i]
            if power < 0:
                denom = denom + '/' + unit
                if power < -1:
                    denom = denom + '**' + str(-power)
            elif power > 0:
                num = num + '*' + unit
                if power > 1:
                    num = num + '**' + str(power)
        if len(num) == 0:
            num = '1'
        else:
            num = num[1:]
        return self.__class__(value = new_value, unit = num + denom)

    def inSIUnits(self):
        """
        Return the quantity with all units reduced to SI-compatible elements.

            >>> e = PhysicalField('2.7 Hartree*Nav')
            >>> print e.inSIUnits().allclose("7088849.01085 kg*m**2/s**2/mol")
            1
        """
        if self.unit.factor != 1:
            return self.inBaseUnits()
        else:
            return self

    def isCompatible (self, unit):
        unit = _findUnit (unit)
        return self.unit.isCompatible (unit)

    def arccos(self):
        """
        Return the inverse cosine of the `PhysicalField` in radians

            >>> print PhysicalField(0).arccos().allclose("1.57079632679 rad")
            1

        The input `PhysicalField` must be dimensionless

            >>> print numerix.round_(PhysicalField("1 m").arccos(), 6)
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return PhysicalField(value=umath.arccos(self.inDimensionless()), unit = "rad")

    def arccosh(self):
        """
        Return the inverse hyperbolic cosine of the `PhysicalField`

            >>> print numerix.allclose(PhysicalField(2).arccosh(),
            ...                        1.31695789692)
            1

        The input `PhysicalField` must be dimensionless

            >>> print numerix.round_(PhysicalField("1. m").arccosh(), 6)
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return umath.arccosh(self.inDimensionless())

    def arcsin(self):
        """
        Return the inverse sine of the `PhysicalField` in radians

            >>> print PhysicalField(1).arcsin().allclose("1.57079632679 rad")
            1

        The input `PhysicalField` must be dimensionless

            >>> print numerix.round_(PhysicalField("1 m").arcsin(), 6)
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return PhysicalField(value = umath.arcsin(self.inDimensionless()), unit = "rad")

    def sqrt(self):
        """
        Return the square root of the `PhysicalField`

            >>> print PhysicalField("100. m**2").sqrt()
            10.0 m

        The resulting unit must be integral

            >>> print PhysicalField("100. m").sqrt()
            Traceback (most recent call last):
                ...
            TypeError: Illegal exponent
        """
        return pow(self, 0.5)

    def sin(self):
        """
        Return the sine of the `PhysicalField`

            >>> print PhysicalField(numerix.pi/6,"rad").sin()
            0.5
            >>> print PhysicalField(30.,"deg").sin()
            0.5

        The units of the PhysicalField must be an angle

            >>> PhysicalField(30.,"m").sin()
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return umath.sin(self.inRadians())

    def sinh(self):
        """
        Return the hyperbolic sine of the `PhysicalField`

            >>> PhysicalField(0.).sinh()
            0.0

        The units of the `PhysicalField` must be dimensionless

            >>> PhysicalField(60.,"m").sinh()
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return umath.sinh(self.inDimensionless())

    def cos(self):
        """
        Return the cosine of the `PhysicalField`

            >>> print numerix.round_(PhysicalField(2*numerix.pi/6,"rad").cos(), 6)
            0.5
            >>> print numerix.round_(PhysicalField(60.,"deg").cos(), 6)
            0.5

        The units of the `PhysicalField` must be an angle

            >>> PhysicalField(60.,"m").cos()
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return umath.cos(self.inRadians())

    def cosh(self):
        """
        Return the hyperbolic cosine of the `PhysicalField`

            >>> PhysicalField(0.).cosh()
            1.0

        The units of the `PhysicalField` must be dimensionless

            >>> PhysicalField(60.,"m").cosh()
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return umath.cosh(self.inDimensionless())

    def tan(self):
        """
        Return the tangent of the `PhysicalField`

            >>> numerix.round_(PhysicalField(numerix.pi/4,"rad").tan(), 6)
            1.0
            >>> numerix.round_(PhysicalField(45,"deg").tan(), 6)
            1.0

        The units of the `PhysicalField` must be an angle

            >>> PhysicalField(45.,"m").tan()
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return umath.tan(self.inRadians())

    def tanh(self):
        """
        Return the hyperbolic tangent of the `PhysicalField`

            >>> print numerix.allclose(PhysicalField(1.).tanh(), 0.761594155956)
            True

        The units of the `PhysicalField` must be dimensionless

            >>> PhysicalField(60.,"m").tanh()
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return umath.tanh(self.inDimensionless())

    def arctan2(self,other):
        """
        Return the arctangent of `self` divided by `other` in radians

            >>> print numerix.round_(PhysicalField(2.).arctan2(PhysicalField(5.)), 6)
            0.380506

        The input `PhysicalField` objects must be in the same dimensions

            >>> print numerix.round_(PhysicalField(2.54, "cm").arctan2(PhysicalField(1., "inch")), 6)
            0.785398

            >>> print numerix.round_(PhysicalField(2.).arctan2(PhysicalField("5. m")), 6)
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return PhysicalField(value=umath.arctan2(self.value,
                                                 self._inMyUnits(other).value), unit="rad")

    def arctan(self):
        """
        Return the arctangent of the `PhysicalField` in radians

            >>> print numerix.round_(PhysicalField(1).arctan(), 6)
            0.785398

        The input `PhysicalField` must be dimensionless

            >>> print numerix.round_(PhysicalField("1 m").arctan(), 6)
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return PhysicalField(value= umath.arctan(self.inDimensionless()), unit = "rad")

    def arctanh(self):
        """
        Return the inverse hyperbolic tangent of the `PhysicalField`

            >>> print PhysicalField(0.5).arctanh()
            0.549306144334

        The input `PhysicalField` must be dimensionless

            >>> print numerix.round_(PhysicalField("1 m").arctanh(), 6)
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return umath.arctanh(self.inDimensionless())


    def log(self):
        """
        Return the natural logarithm of the `PhysicalField`

            >>> print numerix.round_(PhysicalField(10).log(), 6)
            2.302585

        The input `PhysicalField` must be dimensionless

            >>> print numerix.round_(PhysicalField("1. m").log(), 6)
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return umath.log(self.inDimensionless())

    def log10(self):
        """
        Return the base-10 logarithm of the `PhysicalField`

            >>> print numerix.round_(PhysicalField(10.).log10(), 6)
            1.0

        The input `PhysicalField` must be dimensionless

            >>> print numerix.round_(PhysicalField("1. m").log10(), 6)
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        return umath.log10(self.inDimensionless())


    def floor(self):
        """
        Return the largest integer less than or equal to the `PhysicalField`.

            >>> print PhysicalField(2.2,"m").floor()
            2.0 m
        """
        return self.__class__(value = umath.floor(self.value), unit = self.unit)

    def ceil(self):
        """
        Return the smallest integer greater than or equal to the `PhysicalField`.

            >>> print PhysicalField(2.2,"m").ceil()
            3.0 m
        """
        return self.__class__(value = umath.ceil(self.value), unit = self.unit)

    def conjugate(self):
        """
        Return the complex conjugate of the `PhysicalField`.

            >>> print PhysicalField(2.2 - 3j,"ohm").conjugate() == PhysicalField(2.2 + 3j,"ohm")
            True
        """
        return self.__class__(value = umath.conjugate(self.value), unit = self.unit)

    def dot(self, other):
        """
        Return the dot product of `self` with `other`. The resulting
        unit is the product of the units of `self` and `other`.

            >>> v = PhysicalField(((5.,6.),(7.,8.)), "m")
            >>> print PhysicalField(((1.,2.),(3.,4.)), "m").dot(v)
            [ 26.  44.] m**2

        """
        if not isinstance(other,PhysicalField):
            return self.__class__(value = numerix.dot(self.value,other), unit = self.unit)
        value = numerix.dot(self.value,other.value)
        unit = self.unit*other.unit
        if unit.isDimensionless():
            return value*unit.factor
        else:
            return self.__class__(value = value, unit = unit)

    def take(self, indices, axis = 0):
        """
        Return the elements of `self` specified by the elements of `indices`.
        The resulting `PhysicalField` array has the same units as the original.

            >>> print PhysicalField((1.,2.,3.),"m").take((2,0))
            [ 3.  1.] m

        The optional third argument specifies the axis along which the selection
        occurs, and the default value (as in the example above) is 0, the first axis.

            >>> print PhysicalField(((1.,2.,3.),(4.,5.,6.)),"m").take((2,0), axis = 1)
            [[ 3.  1.]
             [ 6.  4.]] m

        """
        return self.__class__(value = numerix.take(self.value, indices, axis), unit = self.unit)

    def put(self, indices, values):
        """
        `put` is the opposite of `take`. The values of `self` at the locations
        specified in `indices` are set to the corresponding value of `values`.

        The `indices` can be any integer sequence object with values suitable for
        indexing into the flat form of `self`. The `values` must be any sequence of
        values that can be converted to the typecode of `self`.

            >>> f = PhysicalField((1.,2.,3.),"m")
            >>> f.put((2,0), PhysicalField((2.,3.),"inch"))
            >>> print f
            [ 0.0762  2.      0.0508] m

        The units of `values` must be compatible with `self`.

            >>> f.put(1, PhysicalField(3,"kg"))
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units
        """
        numerix.put(self.value, indices, self._inMyUnits(values).value)

    @property
    def shape(self):
        """Tuple of array dimensions."""
        from fipy.tools import numerix
        return numerix.getShape(self.value)

    @property
    def itemsize(self):
        return self.value.itemsize

    def reshape(self, shape):
        """
        Changes the shape of `self` to that specified in `shape`

            >>> print PhysicalField((1.,2.,3.,4.),"m").reshape((2,2))
            [[ 1.  2.]
             [ 3.  4.]] m

        The new shape must have the same size as the existing one.

            >>> print PhysicalField((1.,2.,3.,4.),"m").reshape((2,3)) # doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                ...
            ValueError: total size of new array must be unchanged
        """
        return self.__class__(value = numerix.reshape(self.value, shape), unit = self.unit)

    def sum(self, index = 0):
        """
        Returns the sum of all of the elements in `self` along the
        specified axis (first axis by default).

            >>> print PhysicalField(((1.,2.),(3.,4.)), "m").sum()
            [ 4.  6.] m
            >>> print PhysicalField(((1.,2.),(3.,4.)), "m").sum(1)
            [ 3.  7.] m
        """
        return self.__class__(value = numerix.sum(self.value, index), unit = self.unit)

    def allclose(self, other, atol = None, rtol = 1.e-8):
        """
        This function tests whether or not `self` and `other` are equal subject
        to the given relative and absolute tolerances.  The formula used is::

            | self - other | < atol + rtol * | other |

        This means essentially that both elements are small compared to `atol` or
        their difference divided by `other`'s value is small compared to `rtol`.
        """
        other = self._inMyUnits(other)
        if atol is None:
            atol = PhysicalField(1.e-5, self.unit)
        else:
            atol = self._inMyUnits(atol)

        return MA.allclose(self.value, other.value, atol = atol.value, rtol = rtol)

    def allequal(self, other):
        """
        This function tests whether or not `self` and `other` are exactly equal.
        """
        other = self._inMyUnits(other)
        return MA.allequal(self.value, other.value)

class PhysicalUnit:
    """
    A `PhysicalUnit` represents the units of a `PhysicalField`.
    """

    def __init__(self, names, factor, powers, offset=0):
        """
        This class is not generally not instantiated by users of this module,
        but rather it is created in the process of constructing a `PhysicalField`.

        :Parameters:
          - `names`: the name of the unit
          - `factor`: the multiplier between the unit and the fundamental SI unit
          - `powers`: a nine-element `list`, `tuple`, or Numeric_ `array` representing
            the fundamental SI units of ["m", "kg", "s", "A", "K", "mol", "cd", "rad", "sr"]
          - `offset`: the displacement between the zero-point of the unit and the zero-point
            of the corresponding fundamental SI unit.

        .. _Numeric: http://www.numpy.org
        """
        if type(names) == type(''):
            self.names = _NumberDict()
            self.names[names] = 1
        else:
            self.names = names
        self.factor = factor
        self.offset = offset
        self.powers = numerix.array(powers)

    def __repr__(self):
        """
        Return representation of a physical unit

            >>> PhysicalUnit('m',   1.,    [1,0,0,0,0,0,0,0,0])
            <PhysicalUnit m>
        """
        return '<PhysicalUnit ' + self.name() + '>'

    __str__ = __repr__

    def _checkSame(self, other):
        if not isinstance(other, PhysicalUnit):
            if other == 1:
                return self.isDimensionless()
            else:
                raise TypeError, 'PhysicalUnits can only be compared with other PhysicalUnits'
        if not numerix.alltrue(self.powers == other.powers):
            raise TypeError, 'Incompatible units'

    def __eq__(self, other):
        """
        Determine if units are identical

        >>> a = PhysicalField("1. m")
        >>> b = PhysicalField("3. ft")
        >>> a.unit == b.unit
        0
        >>> a.unit == b.inBaseUnits().unit
        1

        Units can only be compared with other units

        >>> a.unit == 3
        Traceback (most recent call last):
            ...
        TypeError: PhysicalUnits can only be compared with other PhysicalUnits
        """
        self._checkSame(other)
        return self.factor == other.factor

    def __lt__(self, other):
        self._checkSame(other)
        return self.factor < other.factor

    def __le__(self, other):
        self._checkSame(other)
        return self.factor <= other.factor

    def __ne__(self, other):
        self._checkSame(other)
        return self.factor != other.factor

    def __gt__(self, other):
        self._checkSame(other)
        return self.factor > other.factor

    def __ge__(self, other):
        self._checkSame(other)
        return self.factor >= other.factor

    def __mul__(self, other):
        """
        Multiply units together

            >>> a = PhysicalField("1. m")
            >>> b = PhysicalField("3. ft")
            >>> a.unit * b.unit == _findUnit('ft*m')
            True
            >>> a.unit * b.inBaseUnits().unit
            <PhysicalUnit m**2>
            >>> c = PhysicalField("1. s")
            >>> d = PhysicalField("3. Hz")
            >>> c.unit * d.unit == _findUnit('Hz*s')
            True
            >>> c.unit * d.inBaseUnits().unit
            <PhysicalUnit 1>

        or multiply units by numbers

            >>> a.unit * 3.
            <PhysicalUnit m*3.0>

        Units must have zero offset to be multiplied

            >>> e = PhysicalField("1. kB")
            >>> f = PhysicalField("25. degC")
            >>> e.unit * f.unit
            Traceback (most recent call last):
                ...
            TypeError: cannot multiply units with non-zero offset
            >>> e.unit * f.inBaseUnits().unit
            <PhysicalUnit kB*K>
        """
        if self.offset != 0 or (isinstance(other,PhysicalUnit) and other.offset != 0):
            raise TypeError, "cannot multiply units with non-zero offset"
        if isinstance(other,PhysicalUnit):
            return PhysicalUnit(self.names+other.names,
                                self.factor*other.factor,
                                self.powers + other.powers)
        else:
            return PhysicalUnit(self.names+{str(other): 1},
                                self.factor*other,
                                self.powers,
                                self.offset * other)

    __rmul__ = __mul__

    def __truediv__(self, other):
        """
        Divide one unit by another

            >>> a = PhysicalField("1. m")
            >>> b = PhysicalField("3. ft")
            >>> a.unit / b.unit
            <PhysicalUnit m/ft>
            >>> a.unit / b.inBaseUnits().unit
            <PhysicalUnit 1>
            >>> c = PhysicalField("1. s")
            >>> d = PhysicalField("3. Hz")
            >>> c.unit / d.unit
            <PhysicalUnit s/Hz>
            >>> c.unit / d.inBaseUnits().unit
            <PhysicalUnit s**2/1>

        or divide units by numbers

            >>> a.unit / 3.
            <PhysicalUnit m/3.0>

        Units must have zero offset to be divided

            >>> e = PhysicalField("1. J")
            >>> f = PhysicalField("25. degC")
            >>> e.unit / f.unit
            Traceback (most recent call last):
                ...
            TypeError: cannot divide units with non-zero offset
            >>> e.unit / f.inBaseUnits().unit
            <PhysicalUnit J/K>
        """
        if self.offset != 0 or (isinstance(other,PhysicalUnit) and other.offset != 0):
            raise TypeError, "cannot divide units with non-zero offset"
        if isinstance(other,PhysicalUnit):
            return PhysicalUnit(self.names-other.names,
                                self.factor/other.factor,
                                self.powers - other.powers)
        else:
            return PhysicalUnit(self.names+{str(other): -1},
                                self.factor/other, self.powers)

    __div__ = __truediv__

    def __rtruediv__(self, other):
        """
        Divide something by a unit

            >>> a = PhysicalField("1. m")
            >>> 3. / a.unit
            <PhysicalUnit 3.0/m>

        Units must have zero offset to be divided

            >>> b = PhysicalField("25. degC")
            >>> 3. / b.unit
            Traceback (most recent call last):
                ...
            TypeError: cannot divide units with non-zero offset
            >>> 3. / b.inBaseUnits().unit
            <PhysicalUnit 3.0/K>
        """
        if self.offset != 0 or (isinstance(other,PhysicalUnit) and other.offset != 0):
            raise TypeError, "cannot divide units with non-zero offset"
        if isinstance(other,PhysicalUnit):
            return PhysicalUnit(other.names-self.names,
                                other.factor/self.factor,
                                map(lambda a,b: a-b,
                                    other.powers, self.powers))
        else:
            return PhysicalUnit({str(other): 1}-self.names,
                                other/self.factor,
                                -self.powers)

    __rdiv__ = __rtruediv__

    def __pow__(self, other):
        """
        Raise a unit to an integer power

            >>> a = PhysicalField("1. m")
            >>> a.unit**2
            <PhysicalUnit m**2>
            >>> a.unit**-2
            <PhysicalUnit 1/m**2>

        Non-integer powers are not supported

            >>> a.unit**0.5
            Traceback (most recent call last):
                ...
            TypeError: Illegal exponent

        Units must have zero offset to be exponentiated

            >>> b = PhysicalField("25. degC")
            >>> b.unit**2
            Traceback (most recent call last):
                ...
            TypeError: cannot exponentiate units with non-zero offset
            >>> b.inBaseUnits().unit**2
            <PhysicalUnit K**2>
        """
        if self.offset != 0:
            raise TypeError, "cannot exponentiate units with non-zero offset"
        if type(other) == type(0):
            return PhysicalUnit(other*self.names, pow(self.factor, other),
                                self.powers*other)

        other = float(other)

        rounded = int(umath.floor(other+0.5))
        if abs(other-rounded) < 1.e-10:
            other = int(other)
            return PhysicalUnit(other*self.names, pow(self.factor, other),
                                self.powers*other)

        inv_exp = 1./other
        rounded = int(umath.floor(inv_exp+0.5))
        if abs(inv_exp-rounded) < 1.e-10:
            if numerix.logical_and.reduce(self.powers % rounded == 0):
                f = pow(self.factor, other)
                p = self.powers // rounded
                if reduce(lambda a, b: a and b,
                          map(lambda x, e=rounded: x%e == 0,
                              self.names.values())):
                    names = self.names // rounded
                else:
                    names = _NumberDict()
                    if f != 1.:
                        names[str(f)] = 1
                    for i in range(len(p)):
                        names[_base_names[i]] = p[i]
                return PhysicalUnit(names, f, p)
            else:
                raise TypeError, 'Illegal exponent'
        raise TypeError, 'Only integer and inverse integer exponents allowed'

    def conversionFactorTo(self, other):
        """
        Return the multiplication factor between two physical units

            >>> a = PhysicalField("1. mm")
            >>> b = PhysicalField("1. inch")
            >>> print numerix.round_(b.unit.conversionFactorTo(a.unit), 6)
            25.4

        Units must have the same fundamental SI units

            >>> c = PhysicalField("1. K")
            >>> c.unit.conversionFactorTo(a.unit)
            Traceback (most recent call last):
                ...
            TypeError: Incompatible units

        If units have different offsets, they must have the same factor

            >>> d = PhysicalField("1. degC")
            >>> c.unit.conversionFactorTo(d.unit)
            1.0
            >>> e = PhysicalField("1. degF")
            >>> c.unit.conversionFactorTo(e.unit)
            Traceback (most recent call last):
                ...
            TypeError: Unit conversion (K to degF) cannot be expressed as a simple multiplicative factor
        """
        if not numerix.alltrue(self.powers == other.powers):
            if self.isDimensionlessOrAngle() and other.isDimensionlessOrAngle():
                return self.factor/other.factor
            else:
                raise TypeError, 'Incompatible units'
        if self.offset != other.offset and self.factor != other.factor:
            raise TypeError, \
                  ('Unit conversion (%s to %s) cannot be expressed ' +
                   'as a simple multiplicative factor') % \
                  (self.name(), other.name())
        return self.factor/other.factor

    def conversionTupleTo(self, other): # added 1998/09/29 GPW
        """
        Return a `tuple` of the multiplication factor and offset between two physical units

            >>> a = PhysicalField("1. K").unit
            >>> b = PhysicalField("1. degF").unit
            >>> [str(numerix.round_(element,6)) for element in b.conversionTupleTo(a)]
            ['0.555556', '459.67']
        """
        if not numerix.alltrue(self.powers == other.powers):
            raise TypeError, 'Incompatible units'

        # let (s1,d1) be the conversion tuple from 'self' to base units
        #   (ie. (x+d1)*s1 converts a value x from 'self' to base units,
        #   and (x/s1)-d1 converts x from base to 'self' units)
        # and (s2,d2) be the conversion tuple from 'other' to base units
        # then we want to compute the conversion tuple (S,D) from
        #   'self' to 'other' such that (x+D)*S converts x from 'self'
        #   units to 'other' units
        # the formula to convert x from 'self' to 'other' units via the
        #   base units is (by definition of the conversion tuples):
        #     ( ((x+d1)*s1) / s2 ) - d2
        #   = ( (x+d1) * s1/s2) - d2
        #   = ( (x+d1) * s1/s2 ) - (d2*s2/s1) * s1/s2
        #   = ( (x+d1) - (d1*s2/s1) ) * s1/s2
        #   = (x + d1 - d2*s2/s1) * s1/s2
        # thus, D = d1 - d2*s2/s1 and S = s1/s2
        factor = self.factor / other.factor
        offset = self.offset - (other.offset * other.factor / self.factor)
        return (factor, offset)

    def isCompatible (self, other):     # added 1998/10/01 GPW
        """
        Returns a list of which fundamental SI units are compatible between
        `self` and `other`

            >>> a = PhysicalField("1. mm")
            >>> b = PhysicalField("1. inch")
            >>> print numerix.allclose(a.unit.isCompatible(b.unit),
            ...                        [True, True, True, True, True, True, True, True, True])
            True
            >>> c = PhysicalField("1. K")
            >>> print numerix.allclose(a.unit.isCompatible(c.unit),
            ...                        [False, True, True, True, False, True, True, True, True])
            True

        """
        return self.powers == other.powers

    def isDimensionless(self):
        """
        Returns `True` if the unit is dimensionless

            >>> PhysicalField("1. m/m").unit.isDimensionless()
            1
            >>> PhysicalField("1. inch").unit.isDimensionless()
            0
        """
        return not numerix.logical_or.reduce(self.powers)

    def isAngle(self):
        """
        Returns `True` if the unit is an angle

            >>> PhysicalField("1. deg").unit.isAngle()
            1
            >>> PhysicalField("1. rad").unit.isAngle()
            1
            >>> PhysicalField("1. inch").unit.isAngle()
            0
        """
        return self.powers[7] == 1 and \
               numerix.add.reduce(self.powers) == 1

    def isInverseAngle(self):
        """
        Returns `True` if the 1 divided by the unit is an angle

            >>> PhysicalField("1. deg**-1").unit.isInverseAngle()
            1
            >>> PhysicalField("1. 1/rad").unit.isInverseAngle()
            1
            >>> PhysicalField("1. inch").unit.isInverseAngle()
            0
        """
        return self.powers[7] == -1 and \
               numerix.add.reduce(self.powers) == -1


    def isDimensionlessOrAngle(self):
        """
        Returns `True` if the unit is dimensionless or an angle

            >>> PhysicalField("1. m/m").unit.isDimensionlessOrAngle()
            1
            >>> PhysicalField("1. deg").unit.isDimensionlessOrAngle()
            1
            >>> PhysicalField("1. rad").unit.isDimensionlessOrAngle()
            1
            >>> PhysicalField("1. inch").unit.isDimensionlessOrAngle()
            0
        """
        return self.isDimensionless() or self.isAngle()

    def setName(self, name):
        """
        Set the name of the unit to `name`

            >>> a = PhysicalField("1. m/s").unit
            >>> a
            <PhysicalUnit m/s>
            >>> a.setName('meterpersecond')
            >>> a
            <PhysicalUnit meterpersecond>
        """
        self.names = _NumberDict()
        self.names[name] = 1

    def name(self):
        """
        Return the name of the unit

            >>> PhysicalField("1. m").unit.name()
            'm'
            >>> (PhysicalField("1. m") / PhysicalField("1. s")
            ...  / PhysicalField("1. s")).unit.name()
            'm/s**2'
        """
        num = ''
        denom = ''
        for unit in self.names.keys():
            power = self.names[unit]
            if power < 0:
                denom = denom + '/' + unit
                if power < -1:
                    denom = denom + '**' + str(-power)
            elif power > 0:
                num = num + '*' + unit
                if power > 1:
                    num = num + '**' + str(power)
        if len(num) == 0:
            num = '1'
        else:
            num = num[1:]
        return num + denom

# Helper functions

def _findUnit(unit):
    """
    Return the `PhysicalUnit` corresponding to `unit`

        >>> _findUnit('m')
        <PhysicalUnit m>
        >>> _findUnit('')
        <PhysicalUnit 1>
        >>> _findUnit(1.)
        <PhysicalUnit 1>
        >>> _findUnit(PhysicalField("4 N*m").unit) == _findUnit('m*N')
        True
        >>> _findUnit(2.)
        Traceback (most recent call last):
            ...
        TypeError: 2.0 is not a unit
    """
##     print unit, type(unit)

    if type(unit) == type(''):
        name = unit.strip()
        if len(name) == 0 or unit == '1':
            unit = _unity
        else:
            unit = eval(name, _unit_table)
            for cruft in ['__builtins__', '__args__']:
                try: del _unit_table[cruft]
                except: pass

    if not isinstance(unit,PhysicalUnit):
        if unit == 1:
            unit = _unity
        else:
            raise TypeError, str(unit) + ' is not a unit'
    return unit

def _round(x):
    if umath.greater(x, 0.):
        return umath.floor(x)
    else:
        return umath.ceil(x)


def _convertValue (value, src_unit, target_unit):
    (factor, offset) = src_unit.conversionTupleTo(target_unit)
    return (value + offset) * factor

def _Scale(quantity, scaling):
    """
    Normalize `quantity` by `scaling`.

    `quantity` can be a `PhysicalField`

        >>> print numerix.round_(_Scale(PhysicalField("1. inch"), PhysicalField("1. mm")), 6)
        25.4

    or a value-unit string convertable to a `PhysicalField`

        >>> print numerix.round_(_Scale("1. inch", PhysicalField("1. mm")), 6)
        25.4

    or a dimensionless number. A dimensionless number is left alone.

        >>> print numerix.round_(_Scale(PhysicalField(2.), PhysicalField("1. mm")), 6)
        2.0

    It is an error for the result to have dimensions.

        >>> print _Scale(PhysicalField("1. s"), PhysicalField("1. mm"))
        Traceback (most recent call last):
            ...
        TypeError: <PhysicalUnit s> and <PhysicalUnit m> are incompatible
    """

    quantity = PhysicalField(quantity)

    if not quantity.unit.isDimensionless():
        scaling = PhysicalField(scaling)
        # normalize quantity to scaling
        # error will be thrown if incompatible
        dimensionless = quantity / scaling
    else:
        # Assume quantity is a dimensionless number and return it.
        # Automatically throws an error if it's not a number.
        dimensionless = quantity

    if isinstance(dimensionless,PhysicalField) and not dimensionless.unit.isDimensionless():
        raise TypeError, repr(quantity.inBaseUnits().unit) + ' and ' \
        + repr(scaling.inBaseUnits().unit) \
        + ' are incompatible'

    return dimensionless

def _isVariable(var):
    from fipy.variables import variable
    return isinstance(var, variable.Variable)

# SI unit definitions

_base_names = ['m', 'kg', 's', 'A', 'K', 'mol', 'cd', 'rad', 'sr']

_base_units = [('m',   PhysicalUnit('m',   1.,    [1,0,0,0,0,0,0,0,0])),
               ('g',   PhysicalUnit('g',   0.001, [0,1,0,0,0,0,0,0,0])),
               ('s',   PhysicalUnit('s',   1.,    [0,0,1,0,0,0,0,0,0])),
               ('A',   PhysicalUnit('A',   1.,    [0,0,0,1,0,0,0,0,0])),
               ('K',   PhysicalUnit('K',   1.,    [0,0,0,0,1,0,0,0,0])),
               ('mol', PhysicalUnit('mol', 1.,    [0,0,0,0,0,1,0,0,0])),
               ('cd',  PhysicalUnit('cd',  1.,    [0,0,0,0,0,0,1,0,0])),
               ('rad', PhysicalUnit('rad', 1.,    [0,0,0,0,0,0,0,1,0])),
               ('sr',  PhysicalUnit('sr',  1.,    [0,0,0,0,0,0,0,0,1])),
               ]

_prefixes = [('Y',  1.e24),
             ('Z',  1.e21),
             ('E',  1.e18),
             ('P',  1.e15),
             ('T',  1.e12),
             ('G',  1.e9),
             ('M',  1.e6),
             ('k',  1.e3),
             ('h',  1.e2),
             ('da', 1.e1),
             ('d',  1.e-1),
             ('c',  1.e-2),
             ('m',  1.e-3),
             ('mu', 1.e-6),
             ('n',  1.e-9),
             ('p',  1.e-12),
             ('f',  1.e-15),
             ('a',  1.e-18),
             ('z',  1.e-21),
             ('y',  1.e-24),
             ]

_unit_table = {}

for unit in _base_units:
    _unit_table[unit[0]] = unit[1]

def _addUnit(name, unit):
    if name in _unit_table:
        raise KeyError, 'Unit ' + name + ' already defined'
    if type(unit) == type(''):
        unit = eval(unit, _unit_table)
        for cruft in ['__builtins__', '__args__']:
            try: del _unit_table[cruft]
            except: pass
    unit.setName(name)
    _unit_table[name] = unit

def _addPrefixed(unit):
    for prefix in _prefixes:
        name = prefix[0] + unit
        _addUnit(name, prefix[1]*_unit_table[unit])


# SI derived units; these automatically get prefixes

_unit_table['kg'] = PhysicalUnit('kg',   1., [0,1,0,0,0,0,0,0,0])

_addUnit('Hz', '1/s')                # Hertz
_addUnit('N', 'm*kg/s**2')           # Newton
_addUnit('Pa', 'N/m**2')             # Pascal
_addUnit('J', 'N*m')                 # Joule
_addUnit('W', 'J/s')                 # Watt
_addUnit('C', 's*A')                 # Coulomb
_addUnit('V', 'W/A')                 # Volt
_addUnit('F', 'C/V')                 # Farad
_addUnit('ohm', 'V/A')               # Ohm
_addUnit('S', 'A/V')                 # Siemens
_addUnit('Wb', 'V*s')                # Weber
_addUnit('T', 'Wb/m**2')             # Tesla
_addUnit('H', 'Wb/A')                # Henry
_addUnit('lm', 'cd*sr')              # Lumen
_addUnit('lx', 'lm/m**2')            # Lux
_addUnit('Bq', '1/s')                # Becquerel
_addUnit('Gy', 'J/kg')               # Gray
_addUnit('Sv', 'J/kg')               # Sievert

del _unit_table['kg']

for unit in _unit_table.keys():
    _addPrefixed(unit)

# Fundamental constants

_unit_table['pi'] = umath.pi
_addUnit('c', '299792458.*m/s')      # speed of light
_addUnit('mu0', '4.e-7*pi*N/A**2')   # permeability of vacuum
_addUnit('eps0', '1/mu0/c**2')       # permittivity of vacuum
_addUnit('Grav', '6.6742e-11*m**3/kg/s**2') # gravitational constant
_addUnit('hplanck', '6.6260693e-34*J*s')     # Planck constant
_addUnit('hbar', 'hplanck/(2*pi)')   # Planck constant / 2pi
_addUnit('e', '1.60217653e-19*C')    # elementary charge
_addUnit('me', '9.1093826e-31*kg')   # electron mass
_addUnit('mp', '1.67262171e-27*kg')  # proton mass
_addUnit('Nav', '6.0221415e23/mol')  # Avogadro number
_addUnit('kB', '1.3806505e-23*J/K')  # Boltzmann constant

_addUnit('gn', '9.80665*m/s**2')     # standard gravitational acceleration
# Time units

_addUnit('min', '60*s')              # minute
_addUnit('h', '60*min')              # hour
_addUnit('d', '24*h')                # day
_addUnit('wk', '7*d')                # week
_addUnit('yr', '365*d')              # year
_addUnit('yrJul', '365.25*d')        # Julian year
_addUnit('yrSid', '365.2564*d')      # sidereal year

# Length units

_addUnit('inch', '2.54*cm')          # inch
_addUnit('ft', '12*inch')            # foot
_addUnit('yd', '3*ft')               # yard
_addUnit('mi', '5280.*ft')           # (British) mile
_addUnit('nmi', '1852.*m')           # Nautical mile
_addUnit('Ang', '1.e-10*m')          # Angstrom
_addUnit('lyr', 'c*yrJul')           # light year
_addUnit('Bohr', '4*pi*eps0*hbar**2/me/e**2')  # Bohr radius

# Area units

_addUnit('ha', '10000*m**2')         # hectare
_addUnit('acres', 'mi**2/640')       # acre
_addUnit('b', '1.e-28*m')            # barn

# Volume units

_addUnit('l', 'dm**3')               # liter
_addUnit('dl', '0.1*l')
_addUnit('cl', '0.01*l')
_addUnit('ml', '0.001*l')
_addUnit('tsp', '4.928922*ml')  # teaspoon
_addUnit('tbsp', '3*tsp')            # tablespoon
_addUnit('floz', '2*tbsp')           # fluid ounce
_addUnit('cup', '8*floz')            # cup
_addUnit('pt', '16*floz')            # pint
_addUnit('qt', '2*pt')               # quart
_addUnit('galUS', '4*qt')            # US gallon
_addUnit('galUK', '4.54609*l')       # British gallon

# Mass units

_addUnit('amu', '1.6605402e-27*kg')  # atomic mass units
_addUnit('lb', '4.5359237e-1*kg')    # pound
_addUnit('oz', 'lb/16')              # ounce
_addUnit('ton', '2000*lb')           # ton

# Force units

_addUnit('dyn', '1.e-5*N')           # dyne (cgs unit)

# Energy units

_addUnit('erg', '1.e-7*J')           # erg (cgs unit)
_addUnit('eV', 'e*V')                # electron volt
_addPrefixed('eV')
_addUnit('Hartree', 'me*e**4/16/pi**2/eps0**2/hbar**2')
_addUnit('invcm', 'hplanck*c/cm')    # Wavenumbers/inverse cm
_addUnit('Ken', 'kB*K')              # Kelvin as energy unit
_addUnit('cal', '4.184*J')           # thermochemical calorie
_addUnit('kcal', '1000*cal')         # thermochemical kilocalorie
_addUnit('cali', '4.1868*J')         # international calorie
_addUnit('kcali', '1000*cali')       # international kilocalorie
_addUnit('Btui', '1055.05585262*J')  # international British thermal unit

# Power units

_addUnit('hpEl', '746*W')            # electric horsepower
_addUnit('hpUK', '745.7*W')            # horsepower

# Pressure units

_addUnit('bar', '1.e5*Pa')           # bar (cgs unit)
_addUnit('atm', '101325.*Pa')        # standard atmosphere
_addUnit('Torr', 'atm/760')          # torr = mm of mercury
_addUnit('psi', 'lb*gn/inch**2')     # pounds per square inch

# Angle units

_addUnit('deg', 'pi*rad/180')        # degrees

# Temperature units -- can't use the 'eval' trick that _addUnit provides
# for degC and degF because you can't add units
kelvin = _findUnit ('K')
_addUnit ('degR', '(5./9.)*K')       # degrees Rankine
_addUnit ('degC', PhysicalUnit (None, 1.0, kelvin.powers, 273.15))
_addUnit ('degF', PhysicalUnit (None, 5./9., kelvin.powers, 459.67))
del kelvin

_unity = eval("m/m", _unit_table)

def _getUnitStrings():

    working_table = _unit_table.copy()

    def _getSortedUnitStrings(unitDict):
        strings = []
        keys = unitDict.keys()
        keys.sort(key=str.lower)
        for key in keys:
            if key in unitDict:
                unit = unitDict[key]
                if isinstance(unit, PhysicalUnit):
                    tmp = PhysicalField(value = 1, unit = unit)
                    strings.append("%10s = %s" % (str(tmp), str(tmp.inBaseUnits())))

                    del working_table[key]

        return strings

    def _deleteFactors(unit):
        for prefix, factor in _prefixes:
            if (prefix + unit.name()) in working_table:
                del working_table[prefix + unit.name()]


    units = []

    units.append("\nBase SI units::\n")
    units.append("\t%s" % ", ".join(_base_names))
    for name in _base_names:
        if name in working_table:
            del working_table[name]
    for name, unit in _base_units:
        _deleteFactors(unit)

    units.append("\nSI prefixes::\n")
    for prefix, factor in _prefixes:
        units.append("%10s = %g" % (prefix, factor))

    units.append("\nUnits derived from SI (accepting SI prefixes)::\n")
    derived = {}
    for key in working_table.keys():
        if key in working_table:
            unit = working_table[key]
            if isinstance(unit, PhysicalUnit) and unit.factor == 1:
                derived[unit.name()] = unit
                _deleteFactors(unit)

    units.extend(_getSortedUnitStrings(derived))

    units.append("\nOther units that accept SI prefixes::\n")
    prefixed = {}
    for key in working_table.keys():
        if key in working_table:
            unit = working_table[key]
            isPrefixed = 1
            if isinstance(unit, PhysicalUnit):
                for prefix, factor in _prefixes:
                    if (prefix + key) not in working_table:
                        isPrefixed = 0
                        break
                if isPrefixed:
                    prefixed[unit.name()] = unit
                    _deleteFactors(unit)

    units.extend(_getSortedUnitStrings(prefixed))

    units.append("\nAdditional units and constants::\n")
    units.extend(_getSortedUnitStrings(working_table))

    return "\n".join(units)

__doc__ += _getUnitStrings()

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
##     print _getUnitStrings()
    _test()
