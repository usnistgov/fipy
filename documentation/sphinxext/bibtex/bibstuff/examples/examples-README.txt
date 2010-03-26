================================
     README: BibStuff examples
================================

:authors: Dylan W. Schwilk and Alan G. Isaac
:web site: http://www.pricklysoft.org
:source code: http://code.google.com/p/bibstuff/
:date: 2009-02-13
:version: 1.0

Some odds and ends files and scripts
====================================

rst_input.txt
	An example reStructuredText document using bibstuff style
	citations for processing with bib4txt.py.

testout.txt
	Example output bibliography from running bib4txt.py on
	rst_input.txt. 
	`> python bib4txt.py -i examples/rst_input.txt -no examples/new-testout.txt examples/example.bib`

jmaker.py
	A little script to translate the list of journal abbreviations
	from `Cambridge Scientific Abstracts`_ to the format used by
	jabbrev.py.  The CSA format file is journals_from_csa.txt and the
	the output produced by running jmaker.py on this file is
	journal_abbreviations.txt. 

schwilk.bib
	An example bibliography database


.. _`Cambridge Scientific Abstracts` : http://www.csa.com/htbin/sjldisp.cgi?filename=/wais/data/srcjnl/biologset

