#!/usr/bin/env python


__docformat__ = 'restructuredtext'

__all__ = []

class _AbstractRepresentation(object):
    def __init__(self, mesh):
        self.mesh = mesh

    def getstate(self):
        """Collect the necessary information to ``pickle`` the `Mesh` to persistent storage.
        """
        raise NotImplemented

    def setstate(self, state):
        """Populate a new `Mesh` from ``pickled`` persistent storage.
        """
        raise NotImplemented

    def repr(self):
        raise NotImplemented
