#File: jass_style.py
"""
Style for the Journal of Artificial Societies and Social Simulation
The JASSS page for authors__ claims they use the `Oxford style`__,
but these are clearly different.

__http://jasss.soc.surrey.ac.uk/admin/submit.html
__http://www.usq.edu.au/library/help/ehelp/ref_guides/oxford.htm

Examples (from the journal website) inconsistent in use of commas, quotes, ...

	HASTIE, R (1986) "Experimental evidence on group accuracy".
	In Jablin F M, Putnam L L, Roberts K H and Porter L W (Eds.)
	Handbook of Organizational Communication: An Interdisciplinary Perspective, Beverly Hills, CA: Sage.

	KALAKOTA R and Whinston A B (1996)
	Frontiers of Electronic Commerce. Reading, MA: Addison-Wesley Publishing Company, Inc..

	KARPINSKI R (1997) Extranets emerge as next challenge for marketers. Netmarketing, April 1997. pp. 1-4.

	LEE H L and Billington C (1992) Managing Supply Chain Inventory: Pitfalls and Opportunities.
	Sloan Management Review, Spring 1992. pp. 65-73.

	RICHIARDI, M, Leombruni, R, Sonnessa, M and Saam, N (2006).
	'A Common Protocol for Agent-Based Social Simulation'.
	Journal of Artificial Societies and Social Simulation 9(1) http://jasss.soc.surrey.ac.uk/9/1/15.html. 

:contact: http://www.american.edu/cas/econ/faculty/isaac/isaac1.htm
:author: Alan G Isaac
:license: MIT (see `license.txt`_)
:date: 2006-08-01

.. _license.txt: ../license.txt
"""
__docformat__ = "restructuredtext en"
__author__  =   "Alan G. Isaac"
__version__ = "0.7"
__needs__ = '2.4'

# import everything from a useful style
from default import *

#####################################################################
############  Override the style choices  ###########################
#####################################################################


######## ADJUST CITATION TEMPLATE FOR NEW STYLE  ###########
######## note: see help for bibstyles.shared.NameFormatter for name details
CITATION_TEMPLATE.update(dict(
indent_left=3,
name_first = 'V |L |f',
name_other = 'v |l |f',
initials = 'f',
max_citation_names = 5,
name_name_sep = (', ',' and '),
names_details_sep = ' ',
article  = '(%(year)s) "%(title)s". *%(journal)s* %(volume)s, %(month)s %(year)s. pp. %(pages)s. %(url)s',
inproceedings  = '(%(year)s) "%(title)s". In %(editor)s (Eds.) *%(booktitle)s*, %(address)s: %(publishers)s',
incollection  = '(%(year)s) "%(title)s". In %(editor)s (Eds.) *%(booktitle)s*, %(address)s: %(publishers)s',
book = '(%(year)s) *%(title)s*. %(address)s: %(publisher)s.',
techreport  = '(%(year)s) "%(title)s". %(institution)s %(type)s %(number)s. %(url)s',
))


# Redefine the CitationManager class, even if "unchanged".
# (This is necessary if you want to change any of the global formatting functions,
#  and usually you change 'format_inline_cite')
class CitationManager(shared.CitationManager):
	def format_inline_cite(self, cite_key_list):
		pass
	################### CITATION FORMATTING ########################
	def get_citation_label(self,entry,citation_template=None):
		return '.. ['+entry.citekey+']\n'  #:TODO: ? allow use of key in place of citekey ?

	#sort_key for sorting list of references
	# (choice of field_list is a formatting decision)
	def sortkey(self,bibentry):
		return make_sort_key(bibentry,['Author','Year'])

