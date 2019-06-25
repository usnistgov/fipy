from __future__ import unicode_literals
__docformat__ = 'restructuredtext'

__all__ = []

import os

from fipy.terms.term import Term
from fipy.terms.explicitSourceTerm import _ExplicitSourceTerm
from fipy.terms import ExplicitVariableError

class _AbstractBinaryTerm(Term):
    def __init__(self, term, other):

        if not isinstance(other, Term):
            other = _ExplicitSourceTerm(coeff=other, var=term.var)

        self.term = term
        self.other = other

        if term.var is None:
            if other.var is None:
                pass
            else:
                raise ExplicitVariableError
        else:
            if other.var is None:
                if isinstance(other, _ExplicitSourceTerm):
                    other.var = term.var
                else:
                    raise ExplicitVariableError

        Term.__init__(self, var=self._vars[0])

    def _addNone(self, arg0, arg1):
        if arg0 is None and arg1 is None:
            return None
        elif arg0 is None:
            return arg1
        elif arg1 is None:
            return arg0
        else:
            return arg0 + arg1

    def __neg__(self):
        r"""
         Negate a `_BinaryTerm`.

           >>> -(__NonDiffusionTerm(coeff=1.) - __NonDiffusionTerm(coeff=2.))
           (__NonDiffusionTerm(coeff=-1.0) + __NonDiffusionTerm(coeff=2.0))

        """

        return (-self.term) + (-self.other)

    def _calcVars(self):
        """Collect (non-redundant) list of all `CellVariables`
        this binary term solves for.

        note: cannot use a set because its order can be different
        on different processors
        """
        ids = [id(v) for v in self.term._vars]
        othervars = [v for v in self.other._vars if id(v) not in ids]
        return self.term._vars + othervars

    @property
    def _vars(self):
        if not hasattr(self, '_internalVars'):
            self._internalVars = self._calcVars()
        return self._internalVars

    @property
    def _transientVars(self):
        return self.term._transientVars + self.other._transientVars

    @property
    def _diffusionVars(self):
        return self.term._diffusionVars + self.other._diffusionVars

    def _checkVar(self, var):
        self.term._checkVar(var)
        self.other._checkVar(var)


from fipy.terms.nonDiffusionTerm import _NonDiffusionTerm
class __NonDiffusionTerm(_NonDiffusionTerm):
    """
    Dummy subclass for tests
    """
    pass

def _test():
    import fipy.tests.doctestPlus
    return fipy.tests.doctestPlus.testmod()

if __name__ == "__main__":
    _test()
