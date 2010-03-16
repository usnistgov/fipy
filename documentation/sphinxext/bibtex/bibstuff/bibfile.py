"""
Provides two classes, BibFile and BibEntry for
accessing the parts of a bibtex database.
BibFile inherits from ``simpleparse.dispatchprocessor``.
To fill a BibFile instance, bfi, call bibgrammar.Parse(src, bfi).


:author: Dylan Schwilk (esp. BibFile)
:contact: http://www.schwilk.org
:author: Alan G Isaac (esp. BibEntry)
:contact: http://www.american.edu/cas/econ/faculty/isaac/isaac1.htm
:copyright: 2006 by Dylan Schwilk and Alan G Isaac
:license: MIT (see `license.txt`_)
:date: 2006-08-05
:requires: Python 2.4+
:TODO: make this framework more general, perhaps along the lines of the btparse_ library in btOOL_.

.. _btOOL: http://www.gerg.ca/software/btOOL/doc/btparse.html
.. _btparse: http://www.gerg.ca/software/btOOL/doc/btparse.html
.. _license.txt: ./license.txt
"""
__docformat__ = "restructuredtext en"
__authors__  = ["Dylan W. Schwilk", "Alan G. Isaac"]
__version__ = '1.13'
__needs__ = '2.4'

# options:
# __strict__ = False allows empty citekeys
__strict__ = False # should we be strict with bibtex format?


####################### IMPORTS #####################################
# import from standard library
import re
import sys

# import dependencies
from simpleparse.dispatchprocessor import dispatch, DispatchProcessor, getString, lines

#bibstuff imports
import bibgrammar
import bibname #ai: shd move all bibname into here? possibly
from bibstyles.shared import reformat_para
import logging
bibfile_logger = logging.getLogger('bibstuff_logger')
#####################################################################

###############  GLOBAL VARIABLES  ##################################
months_en = ('January','February','March','April','May','June',
             'July','August','September','October','November','December')
monthslower_en = [m.lower() for m in months_en]
monthmacros_en = [m[:3] for m in monthslower_en]
MONTH_DICT = dict( zip(monthmacros_en, months_en) )
#####################################################################


class BibEntry(dict):
	"""
	Stores a single bibliographic entry.
	Provides a dictionary interface to the fields:
	field keys are case-insensitive and fields are stored
	in the order added.
	
	:note: 2006-08-10 use 'citekey' instead of 'key' since BibTeX allows a 'key' field
	:note: 2008-03-29 'entry_type' instead of 'type' since BibTeX allows a 'type' field
	"""
	def __init__(self,*args,**kwargs):
		dict.__init__(self,*args,**kwargs)
		self._fields = []
	def __repr__(self):
		"""return string representation of entry
		
		:note: 2006-08-11:eliminate final comma, handle months-> macro and journal macros
		"""
		stringrep = '@%s{%s,\n' % (self.entry_type.upper() , self.citekey)
		try:
			mlen = max( len(key_str) for key_str in self._fields )  # for pretty format
		except ValueError: #no fields (not a true entry)
			mlen = 0
			bibfile_logger.warn("Entry apparently has no fields.")
		field_list = []
		for key in self._fields:
			addbraces = True
			spacer = ' '*(mlen - len(key) )
			val = self[key]
			#handle crossref
			if key == 'crossref':
				try: val = val['citekey'] #might be an entry
				except TypeError: pass    #->must be a string
			elif key == 'journal':
				if val.isalpha() and val.islower(): #:TODO: allow punctuation!!
					addbraces = False  #i.e., assume it is a macro
			elif key == 'month':
				# always use month macros if possible
				if val.lower() in monthslower_en + monthmacros_en:
					val = val[:3].lower()
					addbraces = False
			elif key in ("year","number","volume","chapter"):
				try:
					addbraces = not int(val)
				except:
					pass
			if addbraces:
				val = "{" + val + "}"
			field_list.append("  %-*s = %s" % (mlen, key, val))
		stringrep += ",\n".join(field_list)
		stringrep += '\n}\n'
		return stringrep
	def __setitem__(self, key, val):
		key = key.lower()
		dict.__setitem__(self, key, val)
		if key == "key":
			bibfile_logger.info(
			"Setting 'key' as an entry *field*. (Recall 'citekey' holds the entry id.)")
		if key not in self._fields and key not in ["citekey","entry_type"] and val:
			self._fields.append(key)
	def __getitem__(self, field):  #field is usually a BibTeX field but can be a citekey
		field = field.lower()
		if field == "key":
			bibfile_logger.info(
			"Seeking 'key' as an entry *field*. (Recall 'citekey' holds the entry id.)")
		try:
			result = dict.__getitem__(self, field)
		#:TODO: rethink this decision (but it is used for formatting)
		#:note: 20080331 changed KeyError to return '' instead of None
		except KeyError:
			crossref = self.get('crossref', '')
			if isinstance(crossref, self.__class__):
				result = crossref[field]
			else:
				result = ''
		#:note: 20080331 add handling of month macros
		if field == 'month' and result in monthmacros_en:
			result = MONTH_DICT[result]
		return result
	def __delitem__(self,key) :
		key = key.lower()
		try:
			dict.__delitem__(self, key)
		except KeyError:
			pass
		try:
			self._fields.remove(key)
		except ValueError:
			pass

	def set_entry_type(self, val):
		self["entry_type"] = val.lower()  #:note: entry_type stored as lowercase
	def get_entry_type(self):
		return self["entry_type"]
	entry_type = property(get_entry_type, set_entry_type, None, "property: 'entry_type'")

	def set_citekey(self, val):
		self["citekey"] = val
	def get_citekey(self):
		return self["citekey"]
	citekey = property(get_citekey,set_citekey,None,"property: 'citekey'")

	def get_fields(self):
		return self._fields
	def set_fields(self, lst):
		self._fields = lst
	fields = property(get_fields, set_fields, None, "property: 'fields'")

	def search_fields(self, string_or_compiled, field='', ignore_case=True):
		"""Return MatchObject if string_or_compiled found in entry else None.
		Find regular expression in entry.
		If field is omitted, search is through all fields.
		
		:note: used by BibFile's find_re method, which is used in turn by bibsearch.py
		:Parameters:
		  `string_or_compiled` : string to compile or compiled regex
		    pattern for searching
		  `field` : string
		    field to search in self (default: search all fields)
		"""
		if isinstance(string_or_compiled, str):
			if ignore_case:
				reo = re.compile(string_or_compiled, re.MULTILINE | re.IGNORECASE)
			else:
				reo = re.compile(string_or_compiled, re.MULTILINE)
		else: #must have a compiled regular expression
			reo = string_or_compiled
		if not field: #->try all fields (but not citekey)
			for f in self.get_fields():
				found = reo.search( self[f] )
				if found: break # no need to check more fields
		#:note: CAN test 'field in self' (even though an entry will not raise KeyError! see TODO above)
		#       BUT do not test 'field in self' bc want test for empty fields below
		elif self[field]:
			found = reo.search( self[field] )
		else:
			if field in self:
				bibfile_logger.info("Empty field %s in entry\n%s.\n."%(self,field))
			found = None
		return found

	def format_names(self, names_formatter):
		"""return formatted BibName-object if possible else raw name

		:type `names_formatter`: NamesFormatter
		:note: called by CitationManager in format_citation
		:note: 2006-08-08 no longer sets a `_names` attribute
		:TODO: add default name_template useful for .bib files?
		"""
		bibfile_logger.debug("BibEntry.format_names: arg is:"+str(names_formatter))
		names = self.get_names()  #get a BibName instance (or possibly, a string)
		#keep string if stuck with it
		if isinstance(names,str):
			result = names
		else: #assume a BibName instance
			#ask BibName instance to format itself (and it asks a NamesFormatter to do it)
			result = names.format(names_formatter)
		bibfile_logger.debug("BibEntry.format_names result = "+str(result))
		return result

	def get_names(self, entry_formatter=None, try_fields=None):
		"""return (BibName-object if possible else string)

		:note: 2006-08-09 matching change to `make_names`, no longer sets `self._names`
		"""
		if entry_formatter is None:
			if not try_fields:
				try_fields = ['author','editor','organization']
		return self.make_names(entry_formatter, try_fields=try_fields)

	def make_names(self, entry_formatter=None, try_fields=None):
		"""return (BibName-object if possible else string)
		(from "raw" names).
		
		:change: 2006-08-02 altered to return BibName instance and not set _names
		:note: self returns None if field missing (-> no KeyError)
		:note: this method introduces the only dependence on simpleparse (via bibname)
		:TODO: return BibName instance for each available name field??
		:Parameters:
		  - `entry_formatter`: EntryFormatter instance to provide style information
		  - `try_fields`: list of field names to try sequentially; none empty filed -> name
		"""
		if entry_formatter is None:
			for field in try_fields:
				raw_names = self[field]
				if raw_names:
					break
		else:
			raw_names, field = entry_formatter.pick_raw_names(self,try_fields)
		return  bibname.BibName(raw_names,from_field=field)  #names are in a BibName object

	def format_with(self, entry_formatter):
		bibfile_logger.debug("BibEntry.format_with: arg is:"+str(entry_formatter))
		#ask the EntryFormatter to do it
		return entry_formatter.format_entry(self)


# ----------------------------------------------------------
# Bibfile
# -------
# Data storage for bibtex file
# ----------------------------------------------------------
class BibFile( DispatchProcessor ):
	"""Stores parsed bibtex file.  Access entries by key.

	:note: a BibFile object should simply *store* .bib file parts
	       (a list of entries and a macro map) and provide access
	       to these parts
	"""
	def __init__(self) :
		self.entries = []
		self._macroMap = {}

	def get_entrylist(self, citekeys, discard=True):
		"""Return list, the BibEntry instances that were found
		(and None for entries not found, unless discarded).
		"""
		if not citekeys:
			bibfile_logger.warning("get_entrylist: No keys provided; returning empty cited-entry list.")
			return []
		temp = [ (key,self.get_entry_by_citekey(key)) for key in citekeys ]
		bad_keys = [pair[0] for pair in temp if not pair[1]]
		if bad_keys and discard:
			bibfile_logger.warning("Database entries not found for the following keys:\n"+"\n".join(bad_keys))
		if discard:
			result = [pair[1] for pair in temp if pair[1]]
		else: #keep None when occurs in entry list
			result =  [pair[1] for pair in temp]
		#attach cross references
		for entry in result:
			if entry:
				crossref = entry.get('crossref', None)
				if isinstance(crossref, str):
					crossref = self.get_entry_by_citekey(crossref)
					if crossref:
						entry['crossref'] = crossref
		return result
		
	def get_entry_by_citekey(self, citekey):
		"""Return entry or None."""
		for entry in self.entries:
			if entry.citekey == citekey:
				return entry

	"""PRODUCTION FUNCTIONS:
	for parsing, must provide a function for each production name.
	"""

	def string(self, (tag,start,stop,subtags), buffer ):
		"""Return a string, stripping leading and trailing markers"""
		return buffer[start+1:stop-1]

	def number(self, (tag,start,stop,subtags), buffer ):
		"""return a number as a string"""
		return buffer[start:stop]

	def entry_type( self, (tag,start,stop,subtags), buffer ):
		"""Return the entry type"""
		return getString((tag,start,stop,subtags), buffer)

	def citekey( self, (tag,start,stop,subtags), buffer ):
		"""Return the entry's citekey"""
		return getString((tag,start,stop,subtags), buffer)

	# macro name
	def name(self, (tag,start,stop,subtags), buffer ):
		"""Return lookup on name or name if not in map."""
		return self._macroMap.get(buffer[start:stop],buffer[start:stop])

	def field(self, (tag,start,stop,subtags), buffer ):
		"""Process a bibentry field and return tuple of name, value."""
		str = ''
		for t in subtags[1][3]:
			if(t) :
				str += dispatch(self, t, buffer) # concatenate hashed together strings
		return (dispatch(self, subtags[0], buffer), str)
				
	def entry( self, (tag,start,stop,subtags), buffer ):
		"""Process the bibentry and its children.
		"""
		entry = BibEntry()
		entry.entry_type = dispatch(self, subtags[0], buffer)
		entry.citekey  = dispatch(self, subtags[1], buffer)
		for field in subtags[2][3] :
			#bibfile_logger.debug("entry: ready to add field: "+str(dispatch(self, field, buffer)))
			k,v = dispatch(self, field, buffer)
			#:note: entry will force k to lowercase
			entry[k] = v
		self.entries.append(entry)
	

	def macro( self, (tag,start,stop,subtags), buffer ):
		"""Process a macro entry and add macros to macro map"""
		name, str = dispatch(self, subtags[0], buffer)
		"""
		the_type = getString(subtags[0], buffer)
		if  the_type.upper() != 'STRING' :
			# it looks like  a macro, but is not: could be a regular entry with no key
			lineno = lines(0, start, buffer)+1
			bibfile_logger.warning("Entry at line %d has macro syntax, but entry_type is %s" % (lineno ,  the_type))
			if not __strict__: # we can add a dummy key and treat this entry as a regular entry
				entry = BibEntry()
				entry.entry_type = dispatch(self, subtags[0], buffer)
				entry.citekey  = 'KEY'  # dummy key -- or should we be strict?
				for field in subtags[1][3] :
					k,v = dispatch(self, field, buffer)
					#:note: entry will force k to lowercase
					entry[k] = v
				self.entries.append(entry)
				bibfile_logger.warning("Dummy key added to entry at line %d" % lineno)
		else :  # otherwise it is really a macro entry
			for field in subtags[1][3]:
				name, str = dispatch(self, field, buffer)
				self._macroMap[name] = str  
		"""
		self._macroMap[name] = str  
		

	def preamble( self, (tag,start,stop,subtags), buffer ):
		"""Process the given production and it's children"""
		the_type = getString(subtags[0], buffer)
		lineno = lines(0,start,buffer)+1
		if  the_type.upper() != 'PREAMBLE' :
			bibfile_logger.warning("Entry at line %d has preamble syntax but entry_type is %s" % (lineno,the_type))
		else :
			bibfile_logger.warning("Preamble entry on line %d:" % lineno + "\n" + buffer[start:stop])

	def comment_entry( self, (tag,start,stop,subtags), buffer ):
		"""Process the given production and it's children"""
		the_type = getString(subtags[0], buffer)
		lineno = lines(0,start,buffer)+1
		if  the_type.upper() != 'COMMENT' :
			bibfile_logger.warning("Entry at line %d has comment syntax but entry_type is %s" % (lineno,the_type))
		else :
			bibfile_logger.info("Comment entry on line %d:" % lineno + " " + getString(subtags[1],buffer))

	def search_entries(self, string_or_compiled, field='', ignore_case=True):
		"""Return list of matching entries.
		Search for regular expression in the fields of each entry.
		If field is omitted, search is through all fields.
		
		:note: used by bibsearch.py
		:Parameters:
		  `string_or_compiled` : string to compile or compiled regex
		    pattern for searching
		  `field` : string
		    field to search in self (default: search all fields)
		"""
		if isinstance(string_or_compiled, str):
			if ignore_case:
				reo = re.compile(string_or_compiled, re.MULTILINE | re.IGNORECASE)
			else:
				reo = re.compile(string_or_compiled, re.MULTILINE)
		else: #->must have a compiled regular expression
			reo = string_or_compiled
		"""
		Find regex in bib_entry.
		If field is omitted, search is through all fields.
		
		:note: used by bibsearch.py
		"""
		ls = [entry for entry in self.entries
			if entry.search_fields(string_or_compiled=reo, field=field, ignore_case=ignore_case)]
		return ls


# self test
# -------------------------
# usage: bibfile.py DATABASE_FILE
if __name__ == "__main__":
	import sys
	if len(sys.argv) > 1 :
		src = open(sys.argv[1]).read()
		bfile = BibFile()
		bibgrammar.Parse(src, bfile)
		for entry in bfile.entries :
			print entry

	else :
		print "self test usage: bibfile.py DATABASE_FILE"

