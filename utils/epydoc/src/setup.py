#!/usr/local/bin/python
#
# Distutils setup script for the Natural Language
# Processing Toolkit
#
# Created [05/27/01 09:04 PM]
# Edward Loper
#

from distutils.core import setup
import re, sys, epydoc

VERSION = str(epydoc.__version__)
(AUTHOR, EMAIL) = re.match('^(.*?)\s*<(.*)>$', epydoc.__author__).groups()
URL = epydoc.__url__
LICENSE = epydoc.__license__

if '--format=wininst' in sys.argv:
    SCRIPTS = ['scripts/epydoc.pyw', 'scripts/epydoc.py']
else:
    SCRIPTS = ['scripts/epydoc', 'scripts/epydocgui']

setup(name="epydoc",
      description="Edward Loper's API Documentation Generation Tool",
      version=VERSION,
      author=AUTHOR,
      author_email=EMAIL,
      license=LICENSE,
      url=URL,
      scripts=SCRIPTS,
      packages=['epydoc', 'epydoc.markup', 'epydoc.test'])

