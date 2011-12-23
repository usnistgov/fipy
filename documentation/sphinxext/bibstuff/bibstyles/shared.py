#File: shared.py
"""
Utilities and formatting classes for BibStuff,
especially for bib4txt.py.

:author: Alan G Isaac
:contact: http://www.american.edu/cas/econ/faculty/isaac/isaac1.htm
:copyright: 2008 by Alan G. Isaac
:license: MIT (see `license.txt`_)
:since: 2006-08-01
:date: 2008-10-14

.. _license.txt: ../license.txt
"""
__docformat__ = "restructuredtext en"
__version__ = "1.3"

###################  IMPORTS  ##################################################
#import from standard library
import logging, re
#import dependencies
import simpleparse
# We need to import this specifically because simpleparse does not import it by
# default
import simpleparse.dispatchprocessor
#create globals
shared_logger = logging.getLogger('bibstuff_logger')
################################################################################

from .default_templates import DEFAULT_CITATION_TEMPLATE

#allow for a single citation reference to have keys for multiple citations
#ordinarily, you do not override this
CITE_SEP = ','

def append_sep(s,sep):
	"""return s+sep after removing duplicate punctuation at the join
	`s`: string
	`sep`: string
	TODO? restrict characters removed
	"""
	if s[-1]==sep[0]:
		sep = sep[1:]
	return s+sep

def reformat_para(para='', left=0, right=72, just='LEFT'):
	"""Simple paragraph reformatter.  Allows specification
	of left and right margins, and of justification style
	(using constants defined in module).
	:note: Adopted by Schwilk from David Mertz's example in TPiP
	:see:  Mertz, David,  *Text Processing in Python* (TPiP)
	"""
	LEFT, RIGHT, CENTER = 'LEFT', 'RIGHT', 'CENTER'
	words = para.split()
	lines = []
	line  = ''
	word = 0
	end_words = 0
	while not end_words:
		if len(words[word]) > right-left: # Handle very long words
			line = words[word]
			word +=1
			if word >= len(words):
				end_words = 1
		else:							 # Compose line of words
			while len(line)+len(words[word]) <= right-left:
				line += words[word]+' '
				word += 1
				if word >= len(words):
					end_words = 1
					break
		lines.append(line)
		line = ''
	if just.upper() == CENTER:
		r, l = right, left
		return '\n'.join([' '*left+ln.center(r-l) for ln in lines])
	elif just.upper() == RIGHT:
		return '\n'.join([line.rjust(right) for line in lines])
	else: # left justify
		return '\n'.join([' '*left+line for line in lines])

class NamesFormatter(object):
	"""Provides a formatter for BibName instances.
	Instances are initialized with formatting information.
	Use the `format_names` method to produce
	a formatted string representing a BibName instance.

	Sample usage::

		#create an author entry
		n = bibname.BibName('One, Test and Test Two','author') 
		#create a formatter
		nf = bibstyles.shared.NamesFormatter(template_list=['f| v| l| j']*2,initials=False)
		#print the formatted names
		print nf.format_names(n)

	:see: documentation for the `NameFormatter` class
	:note: 2006-08-03 add initials keyword to ``__init__``
	"""
	def __init__(self, citation_template=None, template_list=None, initials=''):
		"""Create name formatters for each template."""
		shared_logger.debug("NamesFormatter.__init__ args: "+str((citation_template,template_list,initials)))
		assert (template_list or citation_template), "Must provide formatting templates."
		if citation_template:
			self.citation_template = citation_template
			self.template_list = [citation_template['name_first'], citation_template['name_other']]
			self.initials = citation_template['initials']
			self.etal = citation_template['etal']
			self.max_citation_names = citation_template['max_citation_names']
			self.name_name_sep = citation_template['name_name_sep']
		else: #set defaults
			self.template_list = template_list
			self.initials = initials
			self.etal = "et al."
			self.max_citation_names = 99
			self.name_name_sep = (', ', ', and ')
		self.formatters = [ NameFormatter(template,self.initials) for template in self.template_list ]

	#get all names, formatted as a string
	def format_names(self,names):
		"""Return string,
		which represents the BibName instance `names`
		formatted as determined by the `NamesFormatter` attributes.

		`NAME FORMATTING TEMPLATES`_ are explained in some detail
		in the doc string for the NameFormatter class.  Briefly:

		Template sections are separated by ``|``.
		Name parts are referred to by first letter: (v)on, (l)last, (j)r or (f)irst.
		These letters may be followed by token separator enclosed in curly braces.
		Any other characters are included as is.

		:type `names`: BibName object
		:note: 2006-07-25 radically refactored from bibname.py's FormatName() function

		.. _`NAME FORMATTING TEMPLATES`: bibstyles/shared.py
		"""
		shared_logger.debug("NamesFormatter.format: Type of names data is "+str(type(names)))
		#get the list of name_dicts from the BibName instance
		#   each name_dict in the list has the keys: first , von, last, jr
		names_dicts = names.get_names_dicts()
		num_names = len(names_dicts)

		#now make a list of formatted names
		#the first name formatted with the first formatter no matter what
		formatted_name_list = [ self.formatters[0].format_name(names_dicts[0]) ]
		#any additional names are formatted with the second formatter (unless too many -> etal)
		if num_names > 1 and num_names <= self.max_citation_names:
			for name_dict in names_dicts[1:]:  #for each name ...
				formatted_name_list.append( self.formatters[1].format_name(name_dict) )
		shared_logger.debug("NamesFormatter.format_names: formatted_name_list: "+str(formatted_name_list))

		#formatted_name_list = [' '.join(names_dicts[0]['last'])]

		#now concatenate the formatted names into the desired result
		result = formatted_name_list.pop(0)
		#first concatenate all but the last
		while len(formatted_name_list) > 1:
			result = append_sep(result,self.name_name_sep[0]) + formatted_name_list.pop(0)
		#finally, add on the last (with the different name_name_sep)
		if formatted_name_list:
			final_name = formatted_name_list.pop(0)
			if final_name != "others":
				result = append_sep(result,self.name_name_sep[1]) + final_name
			else:
				result = append_sep(result,self.etal)
		assert (len(formatted_name_list) == 0)  #obviously
		if num_names > self.max_citation_names:
			result = append_sep(result,self.etal)
		return result


class NameFormatter(object):
	"""Create a NameFormatter object based on a template string.
	
	NAME FORMATTING TEMPLATES
	-------------------------

	The name template takes some explanation.

	Name parts are referred to by part-designator, which is just the part's first letter:
	(v)on, (l)last, (j)r or (f)irst.
	The designator may be capitalized for force upper-casing the entire part.

	Each name part may have one associated section in a name formatting template.
	Sections are separated by '|' and *must* include a part-designator (one of 'FVLJfvlj').
	The presumption is that part-designators will be the only alphabetic characters in a name template.

	A section will generate output iff the name part for that section exists.
	Each section may have a partsep
	(in curly braces, immediately following the part-designator)
	and other characters
	(which may not be any of 'fvljFVLJ').
	The partsep indicates what should separate multiple tokens of the same part
	(e.g., two part last names, or 'van der' for the (v)on part).
	A part separator will replace the default space to separate multiple tokens in a part.
	Any other characters are included as is.

	For example::

		   "v{~}~|l,| j,| f{. }." with initials='f' produces:
		   "McFeely, J. W." or "van~der~Stadt, Jr, C. M."

	:note: has a property -> must be new style class, inherit carefully
	:note: 20080331 allow capital part-designators (FVLJ) to force capitalization
	"""
	def __init__(self, template, initials=''):
		shared_logger.debug("NameFormatter.__init__ args: "+str((template,initials)))
		#set a default partsep
		#:note: not planning to parameterize this default (e.g., in the citation template)
		self.default_partsep = ' '
		#self.partdict = {}  #this will be set by set_template
		self.initials = initials
		self.set_template(template)

	#get one name, formatted
	def format_name(self,name_data):
		"""Return one name (stored in `name_data`) as a formatted string.

		Formats `name_data` according to the `NameFormatter` template.

		:param `name_data`: list of name_parts or name as string
		:type `name_data`: list or string
		"""
		shared_logger.debug("NameFormatter.format_name:\nType of name_data is: "+str(type(name_data)))
		if isinstance( name_data, (list,tuple) ):
			shared_logger.debug("Assume list is a name_parts list.")
			result = self.name_parts2formatted(name_data)  #TODO: currently commented out for testing dicts
		elif isinstance(name_data, dict):
			shared_logger.debug("Assume dict is a name_dict.")
			result = self.name_dict2formatted(name_data)
		elif isinstance(name_data, str):
			result = name_data
		else:
			raise ValueError("Unrecognized name_data type.")
		shared_logger.debug("NameFormatter.format_name result: '"+result+"'")
		return result

	'''
	def name_parts2formatted(self,name_parts):
		"""Returns one fully formatted name, based on a name_parts tuple.
		"""
		shared_logger.debug("name_parts2formatted: name_parts is "+str(name_parts))
		partdict = self.partdict
		shared_logger.debug("name_parts2formatted: partdict is "+str(partdict))
		result = ''
		#name_parts have a fixed order, and each part is a list (e.g., of one person's last names)
		map_names_parts = dict(f=0, v=1, l=2, j=3)
		if self.initials:
			f,v,l,j = name_parts
			name_parts = ([s[0] for s in f],v,l,j)
		for partcode in partdict['parts_order']:
			partsep = partdict[partcode]['partsep']
			part = partsep.join(name_parts[map_names_parts[partcode]])
			if part:
				result += partdict[partcode]['pre'] + part + partdict[partcode]['post']
			shared_logger.debug("%s: %s"%(partcode,result))
		return result
	'''

	def name_dict2formatted(self,name_dict):
		"""Returns one fully formatted name, based on a name_dict.
		the name_dict should have the keys: first , von, last, jr
		"""
		assert ( len(name_dict['last'][0]) > 0 )
		if name_dict['last'][0] == "others":
			return "others"
		shared_logger.debug("name_dict2formatted: name_dict is "+str(name_dict))
		#get the partdict (that was produced from the name template)
		#  recall that the partdict has keys: pre, post, partsep, parts_order
		#  the parts_order value is a string with characters from "FVLJfvlj"
		partdict = self.partdict
		shared_logger.debug("name_dict2formatted: partdict is "+str(partdict))
		result = ''
		#name_dict has keys, and each value is a list (e.g., of one person's last names)
		map_names_parts = dict(f='first', v='von', l='last', j='jr')
		#change names to initials where requested
		if self.initials:
			name_dict = name_dict.copy()
			for partcode in self.initials.lower():
				part_key = map_names_parts[partcode]
				name_dict[part_key] = [s[0] for s in name_dict[part_key]]
		for partcode in partdict['parts_order']:  #keep the parts in the template determined order
			partsep = partdict[partcode]['partsep']
			part = partsep.join(name_dict[map_names_parts[partcode.lower()]])
			if part:
				#force upper case if parcode is uppercase
				if partcode.isupper():
					part = part.upper()
				result += partdict[partcode]['pre'] + part + partdict[partcode]['post']
			shared_logger.debug("%s: %s"%(partcode,result))
		return result

	def get_template(self):
		return self._template
	def set_template(self,template):
		"""Return None.

		sets the name formatting template *and* sets the associated partdict used for actual formatting 
		"""
		shared_logger.debug("NameFormatter.set_template args: "+str(template))
		assert isinstance(template,str), "Provide a name-template string to make a NameFormatter object."
		self._template = template
		self.partdict = self.template2dict(template)
	template = property(get_template,set_template,None,"template property")

	def template2dict(self,template):
		"""
		parse the name formatting template into a partdict to be used for the actual formatting

		:note: parsing a name template into a partdict is trivial, so just do it here
		:note: allow capital part id (to force capitalization)
		"""
		#to keep track of the order of the parts...
		parts_order = ''
		#split a name template into parts (each part shd have part-designator)
		template_parts = template.split('|')
		partdict = {}
		for part in template_parts:
			for partid in 'FVLJfvlj':
				if partid in part:
					parts_order += partid
					pre, temp = part.split(partid)
					if temp and temp[0] == '{':   #found a partsep
						partsep,post = temp[1:].split('}')
					else:
						post = temp
						partsep = self.default_partsep
					partdict[partid] = dict(pre=pre,post=post,partsep=partsep)
					break
		shared_logger.debug("template2dict: name formatting template parsed to:\n"+str(partdict))
		partdict['parts_order'] = parts_order
		return partdict


class CitationManager(object):
	"""
	:TODO: possibly useful for bibsearch.py
	"""
	default_citation_template = DEFAULT_CITATION_TEMPLATE.copy()

	def __init__(self, biblist, citekeys=None, citation_template=None, sortkey=None):
		self.biblist = biblist
		#:alert: set_citekeys -> self._entries created!
		self.set_citekeys(citekeys)
		if citation_template is None:
			citation_template = self.default_citation_template
		self.citation_template = citation_template
		self.entry_formatter = EntryFormatter(citation_template)
		if sortkey: #TODO: ?? remove this possibility ??
			self.sortkey = sortkey
		self.citeref_processor = None

	def __str__(self):
		if self.citation_template and "citation_sep" in self.citation_template:
			citation_sep = self.citation_template['citation_sep']
		else:
			citation_sep = "\n\n"
		return citation_sep.join( [str(entry)  for entry in self._entries] )

	def set_citeref_processor(self, processor):
		self.citeref_processor = processor
	def format_inline_cite(self, cite_key_list):
		"""Returns a formatted inline citation reference.
		Usually used by a CiteRefProcessor object during processing. 
		Usually styles need to override this method.
		"""
		#substitute formatted citation reference into document text
		self.result.append( self.citation_manager.format_inline_cite(entry_list,cite_key_list) )
		return '**[' + ','.join(cite_key_list) + ']_'


	def get_citekeys(self):
		return self._citekeys
	def set_citekeys(self, citekeys):
		"""set self._citekeys to keys **and** make associated entries
		"""
		shared_logger.debug("shared.CitationManager.set_citekeys %s."%citekeys)
		self._citekeys = citekeys
		if citekeys:
			#discard keys that do not have an entry
			self._entries = self.find_entries(citekeys, discard=True)
		else:
			self._entries = []
	citekeys = property(get_citekeys, set_citekeys, None, "citekeys property")


	def find_entries(self, citekeys=None, discard=True):
		"""return all entries if citekeys==None else matching entries
		discard=True -> discard keys that do not have a bib entry
		"""
		if citekeys is None:
			citekeys = self.citekeys
		result = []
		#TODO: check for reuse of citekeys in different BibFile objects
		for bib in self.biblist:
			result.extend(bib.get_entrylist(citekeys,discard=discard))
		return result
	def get_entries(self, citekeys=None):
		if not citekeys:
			return self._entries[:]
		else:
			return self.find_entries(citekeys)
	#note: citation_rank uses unit-based indexing!! (so styles don't have to offset it)
	def get_citation_rank(self, entry, citekeys=None):
		if citekeys is None:
			citekeys = self._citekeys
		if citekeys is None:  #chk
			citekeys = self.citeref_processor.all_citekeys
			self._citekeys = citekeys
		shared_logger.debug("shared.CitationManager.get_citation_rank citekeys %s."%citekeys)
		if entry.citekey not in citekeys:
			rank = None
			msg = 'Entry citekey not in citekeys; citation_rank set to None.'
			shared_logger.error(msg)
		else: # found the citekey in the cite-key list
			rank = 1 + self._citekeys.index(entry.citekey)
		return rank

	def make_sort_key(self, bibentry, field_list):
		"""create a string for sorting.
		Function returns tuple: (sort_string, bibentry key)

		:note: this is essentially what was Bibstyle's makeSortKey method
		"""
		shared_logger.debug("Entering make_sort_key.")
		result = []
		for field in field_list:
			# some special cases
			if field.lower() in [ 'author','editor','names']:
				result.append(' '.join(bibentry.get_names().get_last_names()).lower())
			elif field.lower() == "year":
				result.append(bibentry['year'])
			else :
				w = bibentry[field]
				if w :
					result.append(w)
		shared_logger.debug("Exiting make_sort_key.")
		return result

	def sortkey(self, entry):
		"""
		:note: the sort key is a style consideration and so must be provided by the style;
			therefore, you must usually OVERRIDE this default sort key
		"""
		result = entry.get_names().get_last_names()
		result.append(entry['year'])
		return result
	def sort(self, sortkey=None): #TODO: not currently using this!
		if sortkey:
			self.sortkey = sortkey  # NB!
		if self.sortkey:
			self._entries.sort(key=sortkey) #2.4 dependency (implements stable Schwartzian transform or better)
			shared_logger.debug("Entries are sorted.")

	#citation_label handling can make be style dependent
	# e.g., for numbered citations, see example_numbered.py
	def get_citation_label(self,entry,citation_template=None):
		return ''

	def make_citations(self, entries=None, citation_template=None):
		"""return formatted citations based on list of entries

		:note: called by ../bib4txt.py in make_text_output
		:note: citation order based on order of entries (so must sort ahead of time)
		:note: related functionality was in the old CitationFormatter's FormatReferences() method
		"""
		shared_logger.debug("shared.CitationManager.make_citations: args are:"+str((entries,citation_template)))
		if entries is None:
			if not self._entries: #get entries matching cite keys found by citeref_processor
				self._entries = self.find_entries(self.citeref_processor.all_citekeys)
			entries = self._entries
			msg = "make_citations: entries are: %s"%(self._entries)
			shared_logger.debug(msg)
		entries.sort(key=self.sortkey)  #TODO!!! use more sensible approach (also: 2.4 dependency)
		if citation_template is None:
			citation_template = self.citation_template
		citation_sep = citation_template['citation_sep']
		#:note: in 2.4 join will accept generators; why is the list necessary?
		result = citation_sep.join( [self.format_citation(entry)  for entry in entries] )
		shared_logger.debug("Exiting make_citations.")
		return result

	def format_citation(self, entry):
		citation_template = self.citation_template
		formatter = self.entry_formatter
		result = formatter.format_entry(entry)
		citation_label = self.get_citation_label(entry, citation_template)
		#result = citation_label + reformat_para( append_sep(names,sep)+details, left=citation_template['indent_left'] )
		result = citation_label + reformat_para( result, left=citation_template['indent_left'] )
		return result






class CiteRefProcessor( simpleparse.dispatchprocessor.DispatchProcessor ):
	"""Formats inline citations and substitutes them into text.
	Stores all cite keys in `all_citekeys` (a list, to record citation order).
	Can store `result` as original text with substituted citation references.

	:note: based on the defunct 'addrefs.py' CitationFormatter class
	"""
	def __init__(self, citation_manager):
		"""
		param `parsed_bibfile`: a dispatch processor holding parsed .bib file
		"""
		#associate with citation manager
		citation_manager.set_citeref_processor(self)
		self.citation_manager = citation_manager
		#self.bib = parsed_bibfile
		# result holds the entire processed file, reformatted for inline citation
		self.result = []
		self.all_citekeys = []  #order matters! unique citekeys added as encountered: see `cite`

	def __repr__(self):
		return ''.join(self.result)

	#set up debug message logging
	def log_msg(self,msg):
		shared_logger.debug(msg)

	#PRODUCTION FUNCTIONS
	# define method for EACH production (see the help for DispatchProcessor)

	def cite(self, (tag,start,stop,subtags), buffer ):
		"""Return everything.

		Alternative default def:
		self.result.append( buffer[start:stop])
		"""
		self.log_msg("The following is parsed as cite:\n" + buffer[start:stop])
		"Process cites and format in text citation according to current style"
		# list because allow for a single citation reference to have keys for multiple citations
		cite_key_list = [s.strip() for s in buffer[start+1:stop-2].split(CITE_SEP)]
		#include current cite keys in set of all cite keys
		#  keep track of order of citation (used by some styles)
		for cite_key in cite_key_list:
			if cite_key not in self.all_citekeys:
				self.all_citekeys.append(cite_key)
		#make (ordered) list of entries for the current cite key(s)
		#:note: need entry to be None if cite_key not found, so discard=False
		entry_list = self.citation_manager.find_entries(cite_key_list,discard=False)
		#substitute formatted citation reference into document text
		self.result.append( self.citation_manager.format_inline_cite(cite_key_list) )

	def inline_literal(self, (tag,start,stop,subtags), buffer):
		"Return everything."
		self.result.append( buffer[start:stop] )
		self.log_msg("The following is parsed as inline_literal:\n" + buffer[start:stop])

	def fn(self, (tag,start,stop,subtags), buffer):
		"Return everything."
		self.result.append( buffer[start:stop])
		self.log_msg("The following is parsed as fn:\n" + buffer[start:stop])

	def plain(self, (tag,start,stop,subtags), buffer):
		"Return everything."
		self.result.append( buffer[start:stop])
		self.log_msg("The following is parsed as plain:\n" + buffer[start:stop])
	

class EntryFormatter(object):
	def __init__(self, citation_template):
		self.citation_template = citation_template
		self.names_formatter=NamesFormatter(citation_template)

	def format_entry(self, entry, citation_template=None):
		"""Return string.
		Format an entry (e.g., as a citation, i.e., a single bibliography reference).
		Note that a BibEntry object acts like a dict for Bib fields
		*except* no KeyError (returns None instead).
		`citation_template` holds templates for entry types

		:note: something related to this method was formerly Bibstyle's formatRef method
		:note: called by make_citations (and currently nothing else)
		"""
		shared_logger.debug("Entering format_citation.")
		if citation_template is None:
			citation_template = self.citation_template
		#:note: a BibEntry object will return None if field is missing
		#get the other (not name) fields
		names = self.format_citation_names(entry, citation_template)
		details = self.format_citation_details(entry, citation_template)
		sep = citation_template['names_details_sep']
		result = append_sep(names, sep) + details
		#ai 2009-02-11 by request but, good idea? think about it
		post_processor = citation_template.get('post_processor', None)
		if post_processor:
			result = post_processor(result)
		shared_logger.debug("EntryFormatter.format_citation: result = "+result)
		return result
	def format_citation_names(self, entry, citation_template=None):
		if citation_template is None:
			citation_template = self.citation_template
		#get the names from the entry (as a BibName object)
		names = entry.make_names(self)  #use this entry formatter (self) to make the names
		#use own names_formatter (based on citation_template) to format the names
		result = self.names_formatter.format_names(names)
		#shared_logger.debug("name_name_sep: "+str(template['name_name_sep']))
		#shared_logger.debug("format_citation_names: result = "+result)
		return result
	#TODO: this deserves substantial enhancement, at the least for journal handling for articles
	def format_citation_details(self, entry, citation_template=None):
		"""Return string."""
		if citation_template is None:
			citation_template = self.citation_template
		try:
			type_template = citation_template[entry.entry_type]  #:note: recall entry_type was stored as lowercase
		except KeyError:  #no template exists for this entry_type -> use default
			type_template = citation_template['default_type']
			shared_logger.warning("Unknown entry type: "+entry.entry_type+". Using default format.")
		#:note: entry will return None instead of KeyError
		result = type_template % entry
		return result
	def pick_raw_names(self, entry, fields=None):
		"""Return BibName-object if possible else string
		(from "raw" names).
		
		:type `field`: str
		:note: 2006-08-02 altered to return BibName instance and not set _names
		:note: self returns None if field missing (-> no KeyError)
		:TODO: return BibName instance for each available name field??
		"""
		names_source = dict(
		article = ['author','organization'],
		book = ['author','editor','organization']
		)
		if fields:
			for field in fields:
				raw_names = entry['field']
				if raw_names:
					break
			if not raw_names:
				shared_logger.warning("EntryFormatter.make_names: empty field -> empty BibName object.")
		#raw_names = self['author'] or self['editor'] #TODO: distinguish author and editor
		elif entry.entry_type in names_source:
			for field in names_source[entry.entry_type]:
				raw_names = entry[field]
				if raw_names:
					break
		else: # default formatting
			for field in ['author','editor','organization']:
				raw_names = entry[field]
				if raw_names:
					break
		if not raw_names:
			shared_logger.warning("No raw names for bib citekey "+entry.citekey)
			raw_names = "Anonymous"  #TODO: shd be a formatting choice (use None?)
			field = None
		#return  bibname.BibName(raw_names,from_field=field)  #names are in a BibName object
		return  raw_names, field
	
