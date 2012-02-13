#! /usr/bin/env python
# -*- coding: latin-1 -*-
# File: isbn2bib.py
"""
:requires: pyaws v.0.3+ (installation is easy; see below)
:requires: free Amazon web services key http://www.amazon.com/gp/browse.html?node=3435361
:license: MIT
:contact: alan dot isaac at gmail dot com 

Installing pyAWS
----------------

You can get a tarball from:http://svn2.assembla.com
If you use SVN, I'll assume you want the latest version.
(Otherwise, get the tagged version rather than the trunk.)

- Decide where you want your pyaws build directory, say in ``mysvn/pyaws``.
- In a command shell, change to the ``mysvn`` directory.
- Issue the command: svn checkout http://svn2.assembla.com/svn/pyaws/trunk/ pyaws
- Change to your new ``pyaws`` directory.
- Use your python to execute: setup.py install

Your Amazon Web Services key
----------------------------

- it is free from Amazon web services
  http://www.amazon.com/gp/browse.html?node=3435361
- right now I only look for it in bibstuff.cfg,
  which must be in the directory from which you
  call your script, and which must contain the lines::

	[isbn2bib]
	aws_key : your_AWS_key_here 

"""
__docformat__ = "restructuredtext en"
__authors__  =    ['Alan G. Isaac']
__version__ =    '0.1'
__needs__ = '2.4'

#prepare a logger
import logging
logging.basicConfig(format='\n%(levelname)s:\n%(message)s\n')
isbn2bib_logger = logging.getLogger('bibstuff_logger')

#we need pyaws to get data from Amazon
from pyaws import ecs


#need an AWS key to proceed (see above)
import ConfigParser as configparser #anticipate name change
cfg = configparser.ConfigParser()
cfg.read('bibstuff.cfg')
aws_key = cfg.get('isbn2bib','aws_key')
try:
	ecs.setLicenseKey(aws_key)
except AWSException:
	print """Failed to set key.
	Do you have a bibstuff.cfg file
	in your current directory?
	"""

#unfortunately, addresses are not available in bookinfo
# hope it's in my list ...
publisher_addresses = dict()
fh = open('data/publisher_addresses.txt','r')
for line in fh:
	if line.startswith('#') or not line.strip():
		continue
	info = tuple(item.strip() for item in line.split('|') )
	try:
		name = info[0].strip()
		address = info[2].strip()
	except:
		continue #TODO: log error
	publisher_addresses[name] = address
fh.close()

def make_entry(isbn):
	"""
	Return a bibfile.BibEntry instance.
	Calls `make_bookdict`;
	called by `main`.

	:date: 2008-08-31
	:todo: this is reusing too much add2bib code
	"""
	import bibfile
	entry = bibfile.BibEntry()
	entry.entry_type = 'book'
	try:
		bkinfo = ecs.ItemLookup(ItemId=isbn, IdType='ISBN',
			SearchIndex="Books", ResponseGroup="Medium")
	except ecs.AWSException:
		print "ItemLookup failed"
		raise
	bkdict = make_bookdict(bkinfo)
	entry.citekey = bkdict['citekey']
	del bkdict['citekey']  #leaving only real field
	#entry.update(bkdict) #TODO: why does this not work?
	for k,v in bkdict.items():
		entry[k] = v
	return entry


def make_bookdict(bkinfo):
	from collections import defaultdict
	import difflib
	bd = defaultdict(str)
	try:
		author = bkinfo.Author.strip()
		author_last = author.split()[-1].lower()
	except AttributeError:
		author = "unknown"
		author_last = "unknown"
	try:
		date = bkinfo.PublicationDate.strip().split()[-1]
		year = date.strip().split('-')[0]
	except AttributeError:
		date = "unknown"
		year = "unknown"
	bd['citekey'] = "%s-%s"%(author_last,date)
	bd['author'] = author
	bd['date'] = date
	bd['year'] = year
	bd['title'] = bkinfo.Title
	bd['isbn'] = bkinfo.ISBN
	publisher = bkinfo.Manufacturer.strip() #?att name??
	bd['publisher'] = publisher
	#thanks to Greg Pinero for nicer address matching:
	best_pub_matches = difflib.get_close_matches(publisher,publisher_addresses.keys(),1)
	if best_pub_matches:
		bd['address'] = publisher_addresses[best_pub_matches[0]]	   
	return bd


# some test ISBNs:
testISBNS = "0-324-23583-6 9780596529321 0231071949"


#-- Command line version of tool
def main():
	"""Command-line tool.
	See bibsearch.py -h for help.
	"""

	import sys
	import add2bib, bibgrammar

	input = sys.stdin
	output = sys.stdout
	
	from optparse import OptionParser
	
	usage = """
	%prog [options]
	example: %prog -f h -bo BIB_DATABASE 0-324-23583-6 9780596529321
	"""


	parser = OptionParser(usage=usage, version ="%prog " + __version__)
	parser.add_option("-f", "--format", action="store",
	                  dest="format", default='b',
					  help="set format(s) of output\nb: BibTeX\nh: HTML\nt: text", metavar="FORMAT")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile",
					  help="Write formatted references to FILE", metavar="FILE")
	parser.add_option("-n", "--nuke", action="store_true", dest="overwrite", default=False,
					  help="CAUTION! silently overwrite outfile, default=%default")
	parser.add_option("-b", "--backup", action="store_true", dest="backup", default=False,
					  help="backup FILE to FILE.bak, default=%default")
	parser.add_option("-v", "--verbose", action="store_true",
	                  dest="verbose", default=False,
					  help="Print INFO messages to stdout, default=%default")
	parser.add_option("-V", "--very_verbose", action="store_true",
	                  dest="very_verbose", default=False,
					  help="Print DEBUG messages to stdout, default=%default")

	#parse options
	(options, args) = parser.parse_args()

	# open output file for writing (default: stdout)
	if options.outfile:
		if options.backup and os.path.exists(options.outfile):
			shutil.copyfile(options.outfile,options.outfile+".bak")
		if options.overwrite or not os.path.exists(options.outfile):
			output = open(options.outfile,'w')
		else:
			isbn2bib_logger.info("Appending to %s.\n(Use -n option to nuke (overwrite) the old output file.)"
			                     %options.outfile)
			output = open(options.outfile,'a')
	print args
	for isbn in args:
		isbn = isbn.replace('-','')
		entry = make_entry(isbn)
		output.write( str(entry) )
		
	if 'h' in options.format:
		output.write( add2bib.html_format(entry) )
	if 't' in options.format:
		output.write( add2bib.text_format(entry) )

if __name__ == '__main__':
	main()

