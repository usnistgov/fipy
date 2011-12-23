#! /usr/bin/env python
# File: biblabel.py
'''
Simple script to generate automatic labels (keys)
for bibtex database entries.
Default format produces citekeys like:
Schwilk+Isaac:2002 and Isaac+Schwilk+etal:2006.

:author: Dylan Schwilk
:contact: http://www.schwilk.org
:author: Alan G Isaac (esp. refactoring)
:contact: http://www.american.edu/cas/econ/faculty/isaac/isaac1.htm
:copyright: 2006 by Dylan Schwilk
:license: MIT (see `license.txt`_)
:date: $Date: 2006/08/29 15:48:05 $

.. _license.txt: ./license.txt
'''
__docformat__ = "restructuredtext en"
__authors__  =    ['Dylan Schwilk','Alan G. Isaac']
__version__ =    '1.9.1'
__needs__ = '2.4'


###################  IMPORTS  ##################################################
# import from standard library
from string import ascii_lowercase
import logging
logging.basicConfig(format='\n%(levelname)s:\n%(message)s\n')
biblabel_logger = logging.getLogger('bibstuff_logger')

# bibstuff imports
import bibfile, bibname
import bibstyles
################################################################################


def make_entry_citekey(entry, used_citekeys, name_template = 'v_|l', max_names = 2, sep = '+' , ysep = ':', etal = 'etal'):
	'''return new entry key (as string)
	'''
	a = entry['author'] or entry['editor']
	if not a:
		a = 'anon'

	y = entry['year'] or '????'

	name_formatter = bibstyles.shared.NameFormatter(template = name_template)
	names_dicts = entry.get_names().get_names_dicts()
	#make list of 'v_|l' last names, which can possibly have multiple tokens (e.g., two piece last names)
	ls = [name_formatter.format_name(name_dict) for name_dict in names_dicts]
	if len(ls) > max_names:
		ls = ls[:max_names] + [etal]
	#for each name, join the tokens with an underscore (i.e., split on whitespace and then join with '_').
	ls = ['_'.join( s.split() )  for s in ls]
	result =  sep.join(ls) + ysep + y

	#make unique result: if needed, append suffix (sfx) b or c or d . . .
	sfx = ''; c = 1
	while result+sfx in used_citekeys:
		sfx = ascii_lowercase[c%26]*(1+c//26)  #:note: lowercase since BibTeX does not distinguish case
		c += 1
	result += sfx
	return result

		  

#-- Command line version of tool
def main():
	'''Command line version of tool'''
	import sys
	import bibgrammar

	from optparse import OptionParser
	
	usage = "%prog [options] filename(s)"

	parser = OptionParser(usage=usage, version ="%prog " + __version__)
	parser.add_option("-y", "--yearsep", action="store", type="string", \
					  dest="yearsep", default = ':', help="char to separate names and year")
	parser.add_option("-s", "--sep", action="store", type="string", \
					  dest="sep",  default = '+', help="char to separate names")
	parser.add_option("-m", "--maxnames", action="store", type="int", \
					  dest="maxnames",  default = 2, help="Max names to add to key")
	parser.add_option("-e", "--etal", action="store", type="string", \
					  dest="etal",  default = 'etal',help="What to add after max names")
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
					  help="Print INFO messages to stdout, default=%default")

	# get options
	(options, args) = parser.parse_args()
	if options.verbose:
		biblabel_logger.setLevel(logging.INFO)

	# get database as text from .bib file(s) or stdin
	if len(args) > 0 :
		try :
		   src = ''.join(open(f).read() for f in args)
		except:
			print 'Error in filelist'
	else :
		src = sys.stdin.read()
	 

	bfile = bibfile.BibFile()
	bibgrammar.Parse(src, bfile)
	used_citekeys = [] # stores created keys
	for entry in bfile.entries:
		label = make_entry_citekey( entry,
		                        used_citekeys,
		                        max_names = options.maxnames,
		                        sep = options.sep,
		                        ysep = options.yearsep,
		                        etal = options.etal)
		#:note: citekey != key, be careful!
		entry.citekey = label
		used_citekeys.insert(0,label) #prepend to take advantage (in make_entry_citekey) of possibly sorted bfile

	for entry in bfile.entries:
		print entry

if __name__ == '__main__':
	main()
