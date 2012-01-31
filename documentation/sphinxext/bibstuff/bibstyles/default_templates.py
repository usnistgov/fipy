#File: default_templates.py
"""
Provides default templates for style writers.
Used by the default style (default.py).

CITATION TEMPLATE
-----------------

book = '%(year)s. %(title)s.',
article  = '%(year)s. %(title)s. %(journal)s %(volume)s, %(pages)s.',
misc = '%(year)s.  %(title)s.',

`default_type` : str
	provide template for string interpolation (use fields as keys)
	e.g.,  '  %(year)s. %(title)s.'
`name_first` : str
	name template for primary name in citation (see NameFormatter documentation)
`name_other` : str
	name template for remaining names in citation (see NameFormatter documentation)
`name_name_sep` : 2-tuple of str
	first element separates each name from the next,
	second element separates penultimate name from ultimate
`etal` : str
	replacement for name when max_citation_names exceeded (e.g., ', et al.')
`initials` : str
	first letter of first (f), von (v), last (l), jr (j) (e.g., 'f')
`max_citation_names` : int
	maximum number of names to format for a citation definition
`indent_left` : int
	left indent for citation definitions
`citation_sep` : str
	separator between citations (e.g.,  "\n\n")
`names_details_sep` : str
	separator between the names and the details in a citation definition (e.g., '. ')
	
:author: Alan G Isaac
:contact: http://www.american.edu/cas/econ/faculty/isaac/isaac1.htm
:copyright: 2006 by Alan G Isaac
:license: MIT (see `license.txt`_)
:since: 2006-08-19

.. _license.txt: ./license.txt
"""
__docformat__ = "restructuredtext en"
import re


#default post processing of citations
# remove newlines (may be returned by bib parser)
_crs = re.compile(r"\s*$\s*", re.MULTILINE)
# remove accents and braces
_aigu = re.compile(r"\\'")
_specialchars = re.compile(r'{\\([a-zA-Z])}')
_deletechars = re.compile(r'[\\{}]')
_ldquote = re.compile(r'``')

def default_post_processor(citations_as_string):
	result = citations_as_string
	result = _crs.sub(' ', result)
	result = _aigu.sub('', result)
	result = _specialchars.sub(r'\1', result)
	result = _deletechars.sub('', result)
	result = _ldquote.sub('"', result)
	return result

DEFAULT_CITEREF_TEMPLATE = dict(
max_cite_names = 2,
citeref_sep = ", ",
)

"""
initials
	string containing none, any, or all of f,v,l,j
:TODO: add separate editor handling
"""

DEFAULT_CITATION_TEMPLATE = dict(
book = '(%(year)s) *%(title)s*. %(address)s: %(publisher)s.',
article  = '%(year)s. %(title)s. *%(journal)s* %(volume)s, %(pages)s.',
techreport  = '(%(year)s) "%(title)s". %(institution)s %(type)s %(number)s. %(url)s',
inproceedings  = '(%(year)s) "%(title)s". In %(editor)s (Eds.) *%(booktitle)s*, %(address)s: %(publisher)s.',
incollection  = '(%(year)s) "%(title)s". In %(editor)s (Eds.) *%(booktitle)s*, %(address)s: %(publisher)s.',
misc = '%(year)s.  %(title)s.',
default_type = '  %(year)s. %(title)s.',
name_first = 'v |l,| j,| f',
name_other = 'f |v |l|, j',
name_name_sep = (', ',', and '),
etal = ', et al.',
initials = '',
max_citation_names = 3,
indent_left = 3,
citation_sep = "\n\n",
names_details_sep = '. ',
post_processor = default_post_processor
)

