from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

from fipy.tools.dimensions.physicalField import PhysicalField

from fipy.tools import numerix

class _ModPhysicalField(PhysicalField):

    def mod(argument):
        return numerix.fmod(argument + 3. * numerix.pi, 2. * numerix.pi) - numerix.pi
    mod = staticmethod(mod)

    def __sub__(self, other):
        if isinstance(other, _ModPhysicalField):
            return self.__class__(value= self.mod(self.inRadians() - other.inRadians()), unit="rad")
        else:
            return self._sum(other, sign2 = lambda b: -b)

    def __rsub__(self, other):
        if isinstance(other, _ModPhysicalField):
            return self.__class__(value = self.mod(argument=other.inRadians() - self.inRadians()), unit="rad")
        else:
            return self._sum(other, sign1 = lambda a: -a)

    def ravel(self):
        return self.value.ravel()
