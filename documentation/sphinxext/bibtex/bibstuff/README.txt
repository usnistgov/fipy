================================
     README: BibStuff (TM)
================================

:authors: Dylan W. Schwilk and Alan G. Isaac
:web site: http://www.pricklysoft.org
:source code: http://code.google.com/p/bibstuff/
:date: 2009-02-13
:version: 1.0

BibStuff License
================

See `license.txt`_, which must be included when this software is
distributed.

Installation
============

Simply type 'python setup.py install' in the unpacked directory.


Command-line tools (scripts)
============================
 
These tools are installed in you python scripts directory or they can
be run directly from where they were unpacked. Each of these tools has
a command line interface and provides the -h option to describe usage.

   * biblabel.py 
      Creates unique keys for entries bibtex database(s).  default keys
      look like Schwilk+Isaac:2006 or Smith+Johnson+etal:1999 Command
      line options allow you to change the default behavior.


   * bibsearch.py
      Search through a bibtex database for entries by key or by
      regular expression.  Results can be output as a (minimally)
      formatted reference, a full bibtex entry, or by key.  Note that
      bibsearch always takes a database by name (-f option or first
      argument) standard input is used for search terms.


   * bib4txt.py
      Creates formatted references for a text dodument.  (Useful for
      reStructuredText documents.) Interacts with a Bibtex style
      database file (without using LaTeX or bibtex).  The source text
      file should include citation references in reStructuredText
      format: a citation key enclosed in brackets, followed by an
      underscore.  Citation keys cannot be all digits.  The source
      document can be output with formatted citation references
      substituted.  In this case, the reference list is added to the
      end of the file.


   * bibname.py
      Create list of author/editor names for a bibtex database.
      Options allow you to specify a name template.  See the module
      documentation for details.

   * jabbrev.py
      Replaces all journal names in a bibtex file with alternative
      names (abbreviations).  The abbreviation file should be in the
      format: <ABBREVIATION> = <LONG_NAME> (see
      /examples/journal_names.txt).  I've also provides a short script
      in the /examples directory that will take the list of journal
      abbreviations at
      http://www.csa.com/htbin/sjldisp.cgi?filename=/wais/data/srcjnl/biologset
      and produce a format readable by jabbrev.py


   * reflist.py
      Creates a list of keys from a latex .bbl file.  This tool simply
      extracts reference keys from the bbl file.  This is useful for
      creating a bibtex database limited to those references which
      occur only in a single latex file.
      
      example: reflist.py my_doc.bbl | bibsearch.py -l my_db.bib > new_db.bib

Modules
=======

Package: bibstuff

Package: bibstuff/bibstyles


Testing
=======

No tests available yet


Related Projects
================

For more information on BibTeX, see the excellent discussion in
chapter 13 section 4 of `The LaTeX Companion`_.

..	_license.txt: ./license.txt

..	_`The LaTeX Companion`: http://www.awprofessional.com/bookstore/product.asp?isbn=0201362996&rl=1
