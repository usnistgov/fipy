#! /usr/bin/env python
#File: bibname.py
"""
Parses bibtex-formatted author/editor raw names and provides
formatting functions (e.g., via bibstyles/shared.NamesFormatter).

:author: Dylan W. Schwilk
:contact: http://www.schwilk.org
:author: Alan G. Isaac
:contact: http://www.american.edu/cas/econ/faculty/isaac/isaac1.htm
:copyright: 2009 by Dylan Schwilk and Alan G Isaac
:license: MIT (see `license.txt`_)
:date: 2009-01-14
:since: 2006-08-29


:note: Major change as of 2008-07-02. Now the ebnf grammar and processor
       handles parsing of a list of names (a bibtex names field such as editor
       or author) and parses the single author name into its fvlj parts. This
       eliminates the need for the original hand-coded parse_raw_names_parts
       function. Moved to using names_dicts rather than names_parts. The
       grammar handles latex accents and ligatures as well as braces strings so
       that a name such as {Barnes and Noble, Inc} is parsed as a single name
       and not split on the " and ".
:todo: The dispatch processor does not currently strip the leading and trailing
       braces from latex/bibtex strings. Not hard to add (see bibfile.py). This
       should be done eventually.
:todo: The grammar does not support quoted strings, only braces strings. Could
       be added fairly simply

.. _license.txt: ./license.txt
"""
__docformat__ = "restructuredtext en"
__authors__  =    ["Dylan W. Schwilk", "Alan G. Isaac"]
__version__ =    '2.0'
__needs__ = '2.4'


################ IMPORTS #############################
# import from standard library
import logging
logging.basicConfig(format='\n%(levelname)s:\n%(message)s\n')
bibname_logger = logging.getLogger('bibstuff_logger')

# import dependencies
import simpleparse
from simpleparse.dispatchprocessor import dispatch, DispatchProcessor, getString, lines
from string import maketrans
# BibStuff imports
import bibstyles, bibfile, bibgrammar
######################################################


################## Global Variables ##################

# constant needed for populating dicts in names_dicts with empty lists for
# missing parts
nameparts = ("first","last","von","jr")

# The EBNF description of a bibtex name field (such as a list of author names).
ebnf_bibname = r"""
namelist := sp*, name, (_and_, name)*
<_and_>  := sp+, "and", sp+
name     := vlf / fvl / fl / vljf / fvlj / l
>l<      :=  last
>vlf<    := (von, sp+)*, last, (sp+, last)*, comma, (sp*, first)+
>fl<     := first, sp+, (first, sp+, ?(capitalized/capstring))*, last
>fvl<    := (first, sp+)+, (von, sp+)+, last, (sp+, last)*
>fvlj<   := fvl, comma, jr
>vljf<   := (von, sp+)*, last, (sp+, last)*, comma, jr,  comma,  first, (sp+ , first)*
von      := lowercase / lowerstring
first    := capitalized / capstring
last     := capitalized / capstring
jr       := "jr" / "Jr" / "JR" /  "Junior" / "junior" /
            "Sr" / "sr" / "II" / "III" / "IV" / "2nd" / "3rd" / "4th"
<comma>           := sp*, ',', sp*
<capitalized>     := capital  , anyc*    
<lowercase>       := ?lowerc, -"and ", anyc*  # Mustn't grab the delimiter _and_ for a part
<ltx_accent>      := '\\`' / "\\'" / '\\^' / '\\"'  /  '\\H' / '\\~' / '\\c' / '\\=' / '\\b' / '\\.' /
                      '\\d' / '\\u' / '\\v' / '\\t'
<ltx_ij_accent>   := '\\^{\\i}' / '\\"{\\i}' / '\\^{\\j}' / '\\"{\\j}'
<ltx_ligature_uc> := '\\AE' / '\\OE' / '\\AA' / '\\O'
<ltx_ligature_lc> := '\\ae' / '\\oe' / '\\aa' / '\\o' / '\\ss'
<capital>         := ('{',capital,'}') / [A-Z] /
                     (ltx_accent, [A-Z]) / (ltx_accent, '{' , [A-Z] , '}') /
                     ltx_ligature_uc
<lowerc>          := ('{',lowerc,'}') / [a-z] / (ltx_accent, [a-z]) /
                     (ltx_accent, '{' , [a-z] , '}') /
                     ltx_ij_accent / ltx_ligature_lc
<anyc>            := [~'-] / capital / lowerc
<string>              :=  '{' , braces_string?, '}'
<capstring>           := '{' , cap_braces_string?, '}'
<lowerstring>         := '{' , lower_braces_string?, '}'
<cap_braces_string>   := ( (capital, -[{}]*) / capstring)+ 
<lower_braces_string> := ( (capital, -[{}]*) / lowerstring)+
<braces_string>       := (-[{}]+ / string)+
<sp>                  := [ \t\n\r.]
"""

bibnamelist_parser = simpleparse.parser.Parser(ebnf_bibname, 'namelist')

######################################################

# ----------- Public Classes and Functions -----------------#


# ----------------------------------------------------------
# BibName
# -------
# Parser processor for bibtex names
# ----------------------------------------------------------
class BibName( simpleparse.dispatchprocessor.DispatchProcessor ):
	"""Processes a bibtex names entry (author, editor, etc) and
	stores the resulting raw_names_parts.
	
	:note: a BibName object should be bibstyle independent.
	"""
	def __init__(self, raw_names=None, from_field=None) :  #:note: 2006-07-25 add initialization based on raw name
		"""initialize a BibName instance
		
		:Parameters:
			`raw_names` : str
				the raw names (e.g., unparsed author field of a BibEntry instance)
			`from_field` : str
				the entry field for the raw name

		:note: 2006-08-02 add `from_field` argument (set by `BibEntry.make_names`)
		"""
		self.from_field = from_field
		self.raw_names = raw_names
		self.names_dicts = []
		#populate self.names_dicts from raw_names
		if raw_names:
			self.parse_raw_names(raw_names)

	###############  PRODUCTION FUNCTIONS  #######################
	# Handle each name by adding new dict to list "names_dicts", then
	# handle each name part by adding to last dict in names_dict list.

	def name(self, (tag,start,stop,subtags), buffer):
		"""Prduction function to process a single name in a nameslist"""
		self.names_dicts.append({}) # add new dict to list
		for part in subtags:
			dispatch(self, part, buffer)
		# Create empty lists for missing parts
		for p in nameparts:
			if not self.names_dicts[-1].has_key(p):
				self.names_dicts[-1][p] = []

	def last(self, (tag,start,stop,subtags), buffer ):
		"""Processes last name part in a single name of a bibtex names field"""
		if self.names_dicts[-1].has_key("last"):
			self.names_dicts[-1]["last"].append(buffer[start:stop])
		else:
			self.names_dicts[-1]["last"] = [buffer[start:stop],]

	def first(self, (tag,start,stop,subtags), buffer ):
		"""Processes first name part in a single name of a bibtex names field"""
		if self.names_dicts[-1].has_key("first"):
			self.names_dicts[-1]["first"].append(buffer[start:stop])
		else:
			self.names_dicts[-1]["first"] = [buffer[start:stop],]

	def von(self, (tag,start,stop,subtags), buffer ): 
		"""Processes von name part in a single name of a bibtex names field"""
		if self.names_dicts[-1].has_key("von"):
			self.names_dicts[-1]["von"].append(buffer[start:stop])
		else:
			self.names_dicts[-1]["von"] = [buffer[start:stop],]

	def jr(self, (tag,start,stop,subtags), buffer ):
		"""Processes jr name part in a single name of a bibtex names field"""
		# Just on jr part so simple add list with one item
		self.names_dicts[-1]["jr"] = [ buffer[start:stop],]
		
	##############  HELPER FUNCTIONS  ######################

	def parse_raw_names(self, raw_name):
		"""This function can be used to populate an empty BibName
		instance or replace all the name values currently contained in
		an instance. It parses the names field with the bibname grammar"""
		self.names_dicts = []  # Replace extant list of  names
		bibnamelist_parser.parse(raw_name,  processor =  self)

	def get_names_dicts(self):  #:note: renamed
		"""
		Return a list of name dicts,
		one dict per name,
		having the fields: first , von, last, jr
		"""
		return self.names_dicts

	
	#ai: method to get last names, which is needed by bibstyle.py and by
	#some style sortkeys
	def get_last_names(self):
		"""Return list of strings, where each string is a last name.
		
		:TODO: graceful handling of missing names parts
		"""
		result = list(' '.join(name_dict['last']) for name_dict in self.names_dicts)
		#bibname_logger.debug("BibName.get_last_names result: "+str(result))
		return result

	def format(self, names_formatter):
		"""
		format a BibName object into a string useful for citations

		:note: called by the BibEntry class in bibfile.py when entry formatting
			is requested
		"""
		return names_formatter.format_names(self)


def getNames(src) :
	"""Returns list of name dicts. Each dict has keys "first", "last",
	"von", "jr". `src` is a string is in bibtex name format.
	"""
	try :
		p = BibName(src)  #:note: 2006-07-25 allow initialization w src
		return p.get_names_dicts()  #:note: 2006-07-25 renamed
	except :
		bibname_logger.error('Error in name %s' % src)
		raise


# command-line version
if __name__ =="__main__":
	import sys
	from optparse import OptionParser
	
	usage = "usage: %prog [options] filenames"

	parser = OptionParser(usage=usage, version ="%prog " + __version__)
	parser.add_option("-t", "--template", action="store", type="string", \
					  dest="template", default = 'f{.}. |v |l| jr', help="Name format template")
	parser.add_option("-i", "--initials", action="store_true", dest="initials", \
					  default = True, help="Initialize first names")
	parser.add_option("-I", "--no-initials", action="store_false", dest="initials", \
					  default = True, help="do not initialize first names")
	parser.add_option("-l", "--last-names", action="store_true", dest="last_names", \
					  default = False, help="Print last names only.")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
					  help="Print INFO messages to stdout, default=%default")

	# get options
	(options, args) = parser.parse_args()
	if options.verbose:
		bibname_logger.setLevel(logging.INFO)
	if options.last_names:
		options.template = 'l'
	if options.initials :
		initials = 'f'  # only first names.  Does any style ever use initials for anything else?
	else :
		initials = ''

	if len(args) == 0 :
		src = sys.stdin.read()
	else :
		flist = list()
		for fname in args:
			try:
				flist.append(open(fname,'r'))
			except IOError :
				bibname_logger.warn('Error in filelist: %s.'%fname)
		src = '\n'.join(f.read() for f in flist)
		map(lambda f: f.close(), flist)

	if not src:
		bibname_logger.error("No bibtex source database found")
		sys.exit(1)
	else:
		bfile = bibfile.BibFile()
		bibgrammar.Parse(src, bfile)

	names_formatter = bibstyles.shared.NamesFormatter(template_list=[options.template]*2,initials=initials)
	for entry in bfile.entries:
		print entry.format_names(names_formatter)

