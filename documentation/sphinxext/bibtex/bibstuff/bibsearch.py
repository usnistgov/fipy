#! /usr/bin/env python
# File: bibsearch.py
"""
Utility for extracting references from bibtex database file.
Extract formatted reference, citekey, or entry from a bibtex file.
Search by key or by regular expression.

bibsearchy.py -h gives usage options.

The script allows style based formatting.  The default style
produces a reference for pasting into a plain text file.

Example::

	python bibsearch.py my_database.bib Smith:1998
	   -> produces a formated reference if citekey Smith:1998 is found
	
	cat ref_list.txt | python bibsearch.py -l my_database.bib
		-> produces a bibtex-format file of all references in list.


:author: Dylan Schwilk
:contact: http://www.schwilk.org
:author: Alan G Isaac
:contact: http://www.american.edu/cas/econ/faculty/isaac/isaac1.htm
:license: MIT (see `license.txt`_)
:date: 2006-08-19
:see: reflist.py (useful in conjuction with bibsearch.py)
:TODO: add additional search capabilities
:TODO: add HTML output option
:TODO: add output options (e.g., to file)

.. _license.txt: ./license.txt
"""
__docformat__ = "restructuredtext en"
__authors__  =    ["Dylan W. Schwilk", "Alan G. Isaac"]
__version__ = "1.8.1"
__needs__ = '2.4'


###################  IMPORTS  ##################################################
#imports from standard library
import string, sys, os
import logging
logging.basicConfig(format='\n%(levelname)s:\n%(message)s\n')
bibsearch_logger = logging.getLogger('bibstuff_logger')

#local imports
import bibfile, bibgrammar
import bibstyles
################################################################################

 
def main():
	"""Command-line tool.
	See bibsearch.py -h for help.
	"""

	from optparse import OptionParser
	
	usage = "usage: %prog [options] FILE [search strings]"
	parser = OptionParser(usage=usage, version ="%prog " + __version__)

	parser.add_option("-k", "--key", action="store_true", dest="citekey_output", 
					  default=False, help="Output citekey rather than reference")
	parser.add_option("-l", "--long", action="store_true", dest="long_output", 
					  default=False, help="Output entire bibtex entry")
	parser.add_option("-r", "--regex", action="store_true", dest="search_input", 
					  default=False, help="Search for regular expression rather than key")
	parser.add_option("-s", "--stylefile", action="store", dest="stylefile", default="default.py",
					  help="Specify user-chosen style file",metavar="FILE")
	parser.add_option("-f", "--field", action="store", type="string", dest="field",
					  default=None,
					  help="Search only FIELD; default=%default.",
					  metavar="FIELD")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
					  help="Print INFO messages to stdout, default=%default")
	parser.add_option("-V", "--very_verbose", action="store_true", dest="very_verbose", default=False,
					  help="Print DEBUG messages to stdout, default=%default")
	
	(options, args) = parser.parse_args()
	if options.verbose:
		bibsearch_logger.setLevel(logging.INFO)
	if options.very_verbose:
		bibsearch_logger.setLevel(logging.DEBUG)
	bibsearch_logger.debug("Script running.\nargs=%s\nstyle file=%s"
	             %(args, options.stylefile)
	            )

	try:
		src = open(args[0]).read()
	except :
		print("Error: No bibtex file found.")
		sys.exit(1)
	# If no search string was sepcified was specified, read search strings from stdin
	if len(args) < 2 :
		searches = string.split(sys.stdin.read())
	else :
		searches = args[1:]

	# create object to store parsed .bib file
	parsed_bibfile = bibfile.BibFile()
	# store a parsed .bib file in parsed_bibfile
	bibgrammar.Parse(src, parsed_bibfile)

	# list of entries
	entrylist = []
	if options.field:
		for s in searches:
			entrylist.extend( parsed_bibfile.search_entries(s, field=options.field) )
	elif options.search_input:
		for s in searches:
			entrylist.extend(parsed_bibfile.search_entries(s))
	else:
		entrylist = parsed_bibfile.get_entrylist(searches, discard=True)

	if entrylist:  #found some matches -> output the list in desired format
		result = ""
		if options.citekey_output:
			result = "\n".join(e.citekey for e in entrylist )
		elif options.long_output :
			result = "\n".join(str(e) for e in entrylist)
		else :
			# style based formated references
			style_stmt = "import bibstyles.%s as style"%os.path.splitext(options.stylefile)[0]
			exec style_stmt in globals()
			citation_manager = style.CitationManager([parsed_bibfile],
													 citekeys=[e.citekey for e in entrylist],
													 citation_template=style.CITATION_TEMPLATE)
			cite_processor = bibstyles.shared.CiteRefProcessor(citation_manager)
			result = citation_manager.make_citations()
		print(result)
	else: #did not find any matches
		bibsearch_logger.info("No matches.")


 
if __name__ == '__main__':
	main()
