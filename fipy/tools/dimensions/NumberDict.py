# Dictionary containing numbers
#
# These objects are meant to be used like arrays with generalized
# indices. Non-existent elements default to zero. Global operations
# are addition, subtraction, and multiplication/division by a scalar.
#
# Written by Konrad Hinsen <hinsen@cnrs-orleans.fr>
# last revision: 1999-7-23
#

__all__ = []

from fipy.tools.dimensions import DictWithDefault

class _NumberDict(DictWithDefault._DictWithDefault):

    """Dictionary storing numerical values

    Constructor: _NumberDict()

    An instance of this class acts like an array of number with
    generalized (non-integer) indices. A value of zero is assumed
    for undefined entries. _NumberDict instances support addition,
    and subtraction with other _NumberDict instances, and multiplication
    and division by scalars.
    """

    def __init__(self):
        DictWithDefault._DictWithDefault.__init__(self, 0)

    def __str__(self):
        return str(self.data)

    def __repr__(self):
        return repr(self.data)

    def __coerce__(self, other):
        if type(other) == type({}):
            new = _NumberDict()
            new.data = other
            other = new
        return self, other

    def __add__(self, other):
        sum = _NumberDict()
        for key in self.keys():
            sum[key] = self[key]
        for key in other.keys():
            sum[key] = sum[key] + other[key]
        return sum

    __radd__ = __add__

    def __sub__(self, other):
        sum = _NumberDict()
        for key in self.keys():
            sum[key] = self[key]
        for key in other.keys():
            sum[key] = sum[key] - other[key]
        return sum

    def __neg__(self):
        neg = _NumberDict()
        for key in list(self.keys()):
            neg[key] = -self[key]
        return neg

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        new = _NumberDict()
        for key in self.keys():
            new[key] = other*self[key]
        return new

    __rmul__ = __mul__

    def __floordiv__(self, other):
        new = _NumberDict()
        for key in self.keys():
            new[key] = self[key]//other
        return new

    __div__ = __floordiv__
