#!/Users/guyer/anaconda/envs/fipy/bin/python2.7

"""futurize doctests

[futurize doesn't run on doctests](https://github.com/PythonCharmers/python-future/issues/103)
This script attempts to make it do so.

Co-opted from 2to3
"""
from __future__ import unicode_literals

import sys
from lib2to3.main import main

fix = sys.argv[1]
fix = fix.split(".")
fixer_pkg, fix = fix[:-1], fix[-1]
fixer_pkg = ".".join(fixer_pkg)
fix = fix.replace("fix_", "", 1)

sys.exit(main(fixer_pkg, ["--doctests_only", "--fix", fix] + sys.argv[2:]))
