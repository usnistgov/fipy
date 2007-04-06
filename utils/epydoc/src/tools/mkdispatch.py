#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-
"""A tool to allow creation of both single page and multi-page manual.

For each file name in argv, detect title names and print a directive
referring to it as anchor in an html file.
"""

# $Id: mkdispatch.py 1547 2007-02-21 17:34:54Z dvarrazzo $
__version__ = "$Revision: 1547 $"[11:-2]
__author__ = "Daniele Varrazzo"
__copyright__ = "Copyright (C) 2007 by Daniele Varrazzo"

import sys, os

def parse_pairs(fn):
    """Parse a file and return a list of directives to create links."""
    outfile = os.path.splitext(os.path.split(fn)[1])[0] + '.html'
    rv = []
    prev = None
    for curr in open(fn):
        curr = curr.rstrip()
        if prev is not None: 
            if curr and curr[0] in "'^-=~":
                if curr == curr[0] * len(curr):
                    rv.append(".. _%s: %s#%s" % 
                        (prev, outfile, get_anchor(prev)))
        prev = curr

    return rv

import string
charmap = {}
charmap.update(zip(string.ascii_lowercase, string.ascii_lowercase))
charmap.update(zip(string.ascii_uppercase, string.ascii_lowercase))
charmap[' '] = '-'
for k in '()':
    charmap[k] = ''

def get_anchor(s):
    # IndexErrors are expected to test for what else include in the map
    return "".join(map(charmap.__getitem__, s))

if __name__ == '__main__':
    for fn in sys.argv[1:]:
        for dir in parse_pairs(fn):
            print dir

