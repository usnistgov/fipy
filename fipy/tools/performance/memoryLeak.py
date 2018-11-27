#!/usr/bin/env python


"""

This python script is ripped from
http://www.nightmare.com/medusa/memory-leaks.html

It outputs the top 100 number of outstanding references for each
object.

"""

__all__ = []

import sys
import types

def _get_refcounts(theClass = None):
    d = {}
    sys.modules
    # collect all classes
    for m in sys.modules.values():
        for sym in dir(m):
            o = getattr (m, sym)
            if type(o) is types.ClassType:
                if theClass is not None and o is not theClass:
                    continue
                d[o] = sys.getrefcount (o)
    # sort by refcount
    pairs = map (lambda x: (x[1],x[0]), d.items())
    pairs.sort()
    pairs.reverse()
    return pairs

def _print_top_N(n = 100, theClass = None):
    for n, c in _get_refcounts(theClass)[:n]:
        print '%10d %s' % (n, c.__name__)

if __name__ == '__main__':
    print_top_N()
