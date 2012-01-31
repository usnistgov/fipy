""" style1

Edited from bibstuff - jasss_style - details:

:contact: http://www.american.edu/cas/econ/faculty/isaac/isaac1.htm
:author: Alan G Isaac
:license: MIT (see `license.txt`_)
:date: 2006-08-01

Modifications Matthew Brett MIT license also
"""

# import everything from a useful style
from bibstuff.bibstyles import default as bbd

######## ADJUST CITATION TEMPLATE FOR NEW STYLE  ###########
######## note: see help for bibstyles.shared.NameFormatter for name details
CITATION_TEMPLATE = bbd.CITATION_TEMPLATE.copy()
CITATION_TEMPLATE.update(dict(
    indent_left=3,
    name_first = 'f |v |l',
    name_other = 'f |v |l',
    initials = '',
    max_citation_names = 5,
    name_name_sep = (', ',' and '),
    names_details_sep = ' ',
    article = '(%(year)s) "%(title)s". *%(journal)s* %(volume)s, %(month)s '
    '%(year)s. pp. %(pages)s. %(url)s %(doi)s',
    inproceedings  = '(%(year)s) "%(title)s". In %(editor)s (Eds.) *%(booktitle)s*, %(address)s: %(publishers)s',
    incollection  = '(%(year)s) "%(title)s". In %(editor)s (Eds.) *%(booktitle)s*, %(address)s: %(publishers)s',
    book = '(%(year)s) *%(title)s*. %(address)s: %(publisher)s.',
    techreport  = '(%(year)s) "%(title)s". %(institution)s %(type)s %(number)s. %(url)s',
))

# Redefine the CitationManager class, even if "unchanged".
# (This is necessary if you want to change any of the global formatting functions,
#  and usually you change 'format_inline_cite')
class CitationManager(bbd.shared.CitationManager):
    default_citation_template = CITATION_TEMPLATE
    def get_citation_label(self,entry,citation_template=None):
        return '.. ['+entry.citekey+'] '
