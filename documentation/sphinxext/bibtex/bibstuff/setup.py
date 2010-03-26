#!/usr/bin/env python

# bibstuff setup script
# Date: 2009-02-13

from distutils.core import setup
setup(name="bibstuff",
      version="1.0.0",  
      description="Bibtex database utilities",
      author="Dylan W. Schwilk and Alan G. Isaac",
      author_email="point@pricklysoft.org",
      url = "http://www.pricklysoft.org/software/bibstuff.html",
      license = "MIT",
      packages = ["bibstyles"],
      py_modules=["bibfile","bibgrammar", "bibname", "ebnf_sp"],
      scripts=["add2bib.py", "bib4txt.py", "bibsearch.py", "biblabel.py", "jabbrev.py", "reflist.py"] )
