"""
Contains some simpleparse style ebnf declarations
for use in parsing with simpleparse.

:author: Alan G Isaac
:contact: http://www.american.edu/cas/econ/faculty/isaac/isaac1.htm
:copyright: 2006 by Alan G Isaac
:license: MIT (see `license.txt`_)
:date: 2006-08-19
"""

# AI's modifications for reST
cites_rest = r"""
src                    := plain_or_fn_or_cite*
>plain_or_fn_or_cite<  := cite / fn_or_plain
cite                   := '[', -([]#] / [0-9]+),-[]]+, ']_'
>fn_or_plain<          := fn / plain
fn                     := '[', ('#' / '*' / [0-9]+), ']_'
plain                  := noref_brackets / nopunct+ / punct
>nopunct<              := -punct
>punct<                :=  '[' / ']'
>noref_brackets<       := '[', -[]]+, ']', ?-'_'
"""

cites_only_rest = r"""
src                    := plain_or_fn_or_cite*
>plain_or_fn_or_cite<  := cite / fn_or_plain
cite                   := '[', -([]#] / [0-9]+),-[]]+, ']_'
>fn_or_plain<          := fn / plain
<fn>                   := '[', ('#' / '*' / [0-9]+), ']_'
<plain>                := noref_brackets / nopunct+ / punct
>nopunct<              := -punct
>punct<                :=  '[' / ']'
>noref_brackets<       := '[', -[]]+, ']', ?-'_'
"""

#EXPERIMENTAL VERSION (use is currently recommended)
cites_xp = r"""
src                    := plain_or_known*
>plain_or_known<       := known / plain
>known<                := inline_literal / cite / fn
inline_literal         := '``', -'``'+, '``'
cite                   := '[', -([]#] / [0-9]+),-[]]+, ']_'
fn                     := '[', ('#' / '*' / [0-9]+), ']_'
plain                  := noref_brackets / nopunct+ / punct
>nopunct<              := ?-inline_literal,-punct
>punct<                :=  '[' / ']'
>noref_brackets<       := '[', -[]]+, ']', ?-'_'
"""



# Schwilk's original:
# EBNF description of a simple text file with citations in reST format
addrefs = r'''
src                 := plain_or_cite*
>plain_or_cite<     := cite / plain
cite                := '[', -[]]+, ']_'
plain               := nocite_brackets / nopunct+ / punct
>nopunct<           := -punct
>punct<             :=  '[' / ']'
>nocite_brackets<   := '[', -[]]+, ']', ?-'_'
'''
