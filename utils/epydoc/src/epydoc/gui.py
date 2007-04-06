#!/usr/bin/env python
#
# objdoc: epydoc command-line interface
# Edward Loper
#
# Created [03/15/02 10:31 PM]
# $Id$
#

"""
Graphical interface to epydoc.  This interface might be useful for
systems where it's inconvenient to use the command-line interface
(such as Windows).  It supports all of the features that are supported
by the command-line interface.  It also supports loading and saving of
X{project files}, which store a set of related modules, and the
options that should be used to generate the documentation for those
modules.

Usage::
    epydocgui [OPTIONS] [FILE.prj | MODULES...]

    FILE.prj                  An epydoc GUI project file.
    MODULES...                A list of Python modules to document.
    -V, --version             Print the version of epydoc.
    -h, -?, --help, --usage   Display this usage message
    --debug                   Do not suppress error messages
"""
__docformat__ = 'epytext en'

import sys, os.path, re, glob
from Tkinter import *
from tkFileDialog import askopenfilename, asksaveasfilename
from thread import start_new_thread, exit_thread
from pickle import dump, load

# askdirectory is only defined in python 2.2+; fall back on
# asksaveasfilename if it's not available.
try: from tkFileDialog import askdirectory
except: askdirectory = None

# Include support for Zope, if it's available.
try: import ZODB
except: pass

##/////////////////////////////////////////////////////////////////////////
## CONSTANTS
##/////////////////////////////////////////////////////////////////////////

DEBUG = 0

# Colors for tkinter display
BG_COLOR='#e0e0e0'
ACTIVEBG_COLOR='#e0e0e0'
TEXT_COLOR='black'
ENTRYSELECT_COLOR = ACTIVEBG_COLOR
SELECT_COLOR = '#208070'
MESSAGE_COLOR = '#000060'
ERROR_COLOR = '#600000'
GUIERROR_COLOR = '#600000'
WARNING_COLOR = '#604000'
HEADER_COLOR = '#000000'

# Convenience dictionaries for specifying widget colors
COLOR_CONFIG = {'background':BG_COLOR, 'highlightcolor': BG_COLOR,
                'foreground':TEXT_COLOR, 'highlightbackground': BG_COLOR}
ENTRY_CONFIG = {'background':BG_COLOR, 'highlightcolor': BG_COLOR,
                'foreground':TEXT_COLOR, 'highlightbackground': BG_COLOR,
                'selectbackground': ENTRYSELECT_COLOR,
                'selectforeground': TEXT_COLOR}
SB_CONFIG = {'troughcolor':BG_COLOR, 'activebackground':BG_COLOR,
             'background':BG_COLOR, 'highlightbackground':BG_COLOR}
LISTBOX_CONFIG = {'highlightcolor': BG_COLOR, 'highlightbackground': BG_COLOR,
                  'foreground':TEXT_COLOR, 'selectforeground': TEXT_COLOR,
                  'selectbackground': ACTIVEBG_COLOR, 'background':BG_COLOR}
BUTTON_CONFIG = {'background':BG_COLOR, 'highlightthickness':0, 'padx':4, 
                 'highlightbackground': BG_COLOR, 'foreground':TEXT_COLOR,
                 'highlightcolor': BG_COLOR, 'activeforeground': TEXT_COLOR,
                 'activebackground': ACTIVEBG_COLOR, 'pady':0}
CBUTTON_CONFIG = {'background':BG_COLOR, 'highlightthickness':0, 'padx':4, 
                  'highlightbackground': BG_COLOR, 'foreground':TEXT_COLOR,
                  'highlightcolor': BG_COLOR, 'activeforeground': TEXT_COLOR,
                  'activebackground': ACTIVEBG_COLOR, 'pady':0,
                  'selectcolor': SELECT_COLOR}
SHOWMSG_CONFIG = CBUTTON_CONFIG.copy()
SHOWMSG_CONFIG['foreground'] = MESSAGE_COLOR
SHOWWRN_CONFIG = CBUTTON_CONFIG.copy()
SHOWWRN_CONFIG['foreground'] = WARNING_COLOR
SHOWERR_CONFIG = CBUTTON_CONFIG.copy()
SHOWERR_CONFIG['foreground'] = ERROR_COLOR

# Colors for the progress bar
PROGRESS_HEIGHT = 16
PROGRESS_WIDTH = 200
PROGRESS_BG='#305060'
PROGRESS_COLOR1 = '#30c070'
PROGRESS_COLOR2 = '#60ffa0'
PROGRESS_COLOR3 = '#106030'

# On tkinter canvases, where's the zero coordinate?
if sys.platform.lower().startswith('win'):
    DX = 3; DY = 3
    DH = 0; DW = 7
else:
    DX = 1; DY = 1
    DH = 1; DW = 3

# How much of the progress is in each subtask?
IMPORT_PROGRESS = 0.1
BUILD_PROGRESS  = 0.2
WRITE_PROGRESS  = 1.0 - BUILD_PROGRESS - IMPORT_PROGRESS

##/////////////////////////////////////////////////////////////////////////
## IMAGE CONSTANTS
##/////////////////////////////////////////////////////////////////////////

UP_GIF = '''\
R0lGODlhCwAMALMAANnZ2QDMmQCZZgBmZgAAAAAzM////////wAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAACH5BAEAAAAALAAAAAALAAwAAAQjEMhJKxCW4gzCIJxXZIEwFGDlDadqsii1sq1U0nA64+ON
5xEAOw==
'''
DOWN_GIF = '''\
R0lGODlhCwAMALMAANnZ2QDMmQCZZgBmZgAAAAAzM////////wAAAAAAAAAAAAAAAAAAAAAAAAAA
AAAAACH5BAEAAAAALAAAAAALAAwAAAQmEIQxgLVUCsppsVPngVtXEFfIfWk5nBe4xuSL0tKLy/cu
7JffJQIAOw==
'''
LEFT_GIF='''\
R0lGODlhDAALAKIAANnZ2QDMmQCZZgBmZgAAAAAzM////////yH5BAEAAAAALAAAAAAMAAsAAAM4
CLocgaCrESiDoBshOAoAgBEyMzgAEIGCowsiOLoLgEBVOLoIqlSFo4OgC1RYM4Ogq1RYg6DLVJgA
Ow==
'''
RIGHT_GIF='''\
R0lGODlhDAALAKIAANnZ2QDMmQBmZgCZZgAzMwAAAP///////yH5BAEAAAAALAAAAAAMAAsAAAM5
GIGgyzIYgaCrIigTgaALIigyEQiqKLoTgaAoujuDgKJLVAgqIoJEBQAIIkKEhaArRFgIukqFoMsJ
ADs=
'''

##/////////////////////////////////////////////////////////////////////////
## MessageIO
##/////////////////////////////////////////////////////////////////////////

class _PipeIO:
    """
    A file-like class that captures strings and prints them.  For
    thread safety, it monotonically increases _data (is this
    necessary?).
    """
    def __init__(self):
        self.clear()

    def clear(self):
        self._data = []
        self._n = 0
        self.softspace = 0

    def write(self, str):
        # If we detect an exception, then start writing everything to
        # the real stderr.
        global DEBUG
        if str.startswith('Traceback') or str.startswith('Exception'):
            DEBUG = 1
        if DEBUG: sys.__stderr__.write(str)

        # Save the data.
        self._data.append(str)
        
    def read(self):
        if self._n >= len(self._data):
            return None
        else:
            str = self._data[self._n]
            self._n += 1
            return str

##/////////////////////////////////////////////////////////////////////////
## THREADED DOCUMENTER
##/////////////////////////////////////////////////////////////////////////

def document(options, progress, cancel):
    """
    Create the documentation for C{modules}, using the options
    specified by C{options}.  C{document} is designed to be started in
    its own thread by L{EpydocGUI._go}.

    @param options: The options to use for generating documentation.
        This includes keyword options that can be given to
        L{html.HTMLFormatter}, as well as the option C{outdir}, which
        controls where the output is written to.
    @type options: C{dictionary}
    @param progress: This first element of this list will be modified
        by C{document} to reflect its progress.  This first element
        will be a number between 0 and 1 while C{document} is creating
        the documentation; and the string C{"done"} once it finishes
        creating the documentation.
    @type progress: C{list}
    """
    progress[0] = 0.02
    from epydoc.html import HTMLFormatter
    from epydoc.objdoc import DocMap, report_param_mismatches
    from epydoc.imports import import_module, find_modules
    from epydoc.objdoc import set_default_docformat

    # Set the default docformat.
    set_default_docformat(options.get('docformat', 'epytext'))

    try:
        # Import the modules.
        modnames = options['modules']
        # First, expand packages.
        for name in modnames[:]:
            if os.path.isdir(name):
                modnames.remove(name)
                new_modules = find_modules(name)
                if new_modules: modnames += new_modules
                sys.stderr.write('!!!Error: %r is not a pacakge\n!!!' % name)
                
        modules = []
        for (name, num) in zip(modnames, range(len(modnames))):
            if cancel[0]: exit_thread()
            sys.stderr.write('***Importing %s\n***' % name)
            try:
                module = import_module(name)
                if module not in modules: modules.append(module)
            except ImportError, e:
                sys.stderr.write('!!!Error importing %s: %s\n!!!' % (name, e))
            progress[0] += (IMPORT_PROGRESS*0.98)/len(modnames)
        
        # Create the documentation map.
        verbosity = 0
        document_bases = 1
        document_autogen_vars = 1
        inheritance_groups = (options['inheritance'] == 'grouped')
        inherit_groups = (options['inheritance'] != 'grouped')
        d = DocMap(verbosity, document_bases, document_autogen_vars,
                   inheritance_groups, inherit_groups)
    
        # Build the documentation.
        for (module, num) in zip(modules, range(len(modules))):
            if cancel[0]: exit_thread()
            sys.stderr.write('***Building docs for %s\n***' % module.__name__)
            d.add(module)
            progress[0] += (BUILD_PROGRESS*0.98)/len(modules)

        # Report any param inheritance mismatches.
        report_param_mismatches(d)

        # Create the HTML formatter.
        htmldoc = HTMLFormatter(d, **options)
        numfiles = htmldoc.num_files()
    
        # Set up the progress callback.
        n = htmldoc.num_files()
        def progress_callback(path, numfiles=numfiles,
                              progress=progress, cancel=cancel, num=[1]):
            if cancel[0]: exit_thread()

            # Find the short name of the file we're writing.
            (dir, file) = os.path.split(path)
            (root, d) = os.path.split(dir)
            if d in ('public', 'private'): fname = os.path.join(d, file)
            else: fname = file
                
            sys.stderr.write('***Writing %s\n***' % fname)
            progress[0] += (WRITE_PROGRESS*0.98)/numfiles
            num[0] += 1
    
        # Write the documentation.
        htmldoc.write(options.get('outdir', 'html'), progress_callback)
    
        # We're done.
        sys.stderr.write('***Finished!\n***')
        progress[0] = 'done'

    except SystemExit:
        # Cancel.
        sys.stderr.write('!!!Cancel!\n!!!')
        progress[0] = 'cancel'
        raise
    except Exception, e:
        # We failed.
        sys.stderr.write('!!!Internal error: %s\n!!!' % e)
        progress[0] = 'cancel'
        raise
    except:
        # We failed.
        sys.stderr.write('!!!Internal error\n!!!')
        progress[0] = 'cancel'
        raise
    
##/////////////////////////////////////////////////////////////////////////
## GUI
##/////////////////////////////////////////////////////////////////////////

class EpydocGUI:
    """
    A graphical user interace to epydoc.
    """
    def __init__(self):
        self._afterid = 0
        self._progress = [None]
        self._cancel = [0]
        self._filename = None

        # Store a copy of sys.modules, so that we can restore it
        # later.  This is useful for making sure that we reload
        # everything when we re-build its documentation.  This will
        # *not* reload the modules that are present when the EpydocGUI
        # is created, but that should only contain some builtins, some
        # epydoc modules, Tkinter, pickle, and thread..
        self._old_modules = sys.modules.keys()

        # Create the main window.
        self._root = Tk()
        self._root['background']=BG_COLOR
        self._root.bind('<Control-q>', self.destroy)
        self._root.bind('<Alt-q>', self.destroy)
        self._root.bind('<Alt-x>', self.destroy)
        #self._root.bind('<Control-d>', self.destroy)
        self._root.title('Epydoc')
        self._rootframe = Frame(self._root, background=BG_COLOR,
                               border=2, relief='raised')
        self._rootframe.pack(expand=1, fill='both', padx=2, pady=2)

        # Set up the basic frames.  Do not pack the options frame or
        # the messages frame; the GUI has buttons to expand them.
        leftframe = Frame(self._rootframe, background=BG_COLOR)
        leftframe.pack(expand=1, fill='both', side='left')
        optsframe = Frame(self._rootframe, background=BG_COLOR)
        mainframe = Frame(leftframe, background=BG_COLOR)
        mainframe.pack(expand=1, fill='both', side='top')
        ctrlframe = Frame(mainframe, background=BG_COLOR)
        ctrlframe.pack(side="bottom", fill='x', expand=0)
        msgsframe = Frame(leftframe, background=BG_COLOR)

        self._optsframe = optsframe
        self._msgsframe = msgsframe

        # Initialize all the frames, etc.
        self._init_menubar()
        self._init_progress_bar(mainframe)
        self._init_module_list(mainframe)
        self._init_options(optsframe, ctrlframe)
        self._init_messages(msgsframe, ctrlframe)
        self._init_bindings()

        # Replace stderr with a _PipeIO stream, so we can capture
        # messages.
        self._errpipe = _PipeIO()
        sys.stderr = self._errpipe

        # Open the messages pane by default.
        self._messages_toggle()

        ## For testing options:
        #self._options_toggle()
        
    def _init_menubar(self):
        menubar = Menu(self._root, borderwidth=2,
                       background=BG_COLOR,
                       activebackground=BG_COLOR)
        filemenu = Menu(menubar, tearoff=0)
        filemenu.add_command(label='New Project', underline=0,
                             command=self._new,
                             accelerator='Ctrl-n')
        filemenu.add_command(label='Open Project', underline=0,
                             command=self._open,
                             accelerator='Ctrl-o')
        filemenu.add_command(label='Save Project', underline=0,
                             command=self._save,
                             accelerator='Ctrl-s')
        filemenu.add_command(label='Save As..', underline=5,
                             command=self._saveas,
                             accelerator='Ctrl-a')
        filemenu.add_separator()
        filemenu.add_command(label='Exit', underline=1,
                             command=self.destroy,
                             accelerator='Ctrl-x')
        menubar.add_cascade(label='File', underline=0, menu=filemenu)
        gomenu = Menu(menubar, tearoff=0)
        gomenu.add_command(label='Run Epydoc',  command=self._open,
                           underline=0, accelerator='Alt-g')
        menubar.add_cascade(label='Run', menu=gomenu, underline=0)
        self._root.config(menu=menubar)
        
    def _init_module_list(self, mainframe):
        mframe1 = Frame(mainframe, relief='groove', border=2,
                        background=BG_COLOR)
        mframe1.pack(side="top", fill='both', expand=1, padx=4, pady=3)
        l = Label(mframe1, text="Modules to document:",
                  justify='left', **COLOR_CONFIG) 
        l.pack(side='top', fill='none', anchor='nw', expand=0)
        mframe2 = Frame(mframe1, background=BG_COLOR)
        mframe2.pack(side="top", fill='both', expand=1)
        mframe3 = Frame(mframe1, background=BG_COLOR)
        mframe3.pack(side="bottom", fill='x', expand=0)
        self._module_list = Listbox(mframe2, width=80, height=10,
                                    selectmode='multiple',
                                    **LISTBOX_CONFIG)
        self._module_list.pack(side="left", fill='both', expand=1)
        sb = Scrollbar(mframe2, orient='vertical',**SB_CONFIG)
        sb['command']=self._module_list.yview
        sb.pack(side='right', fill='y')
        self._module_list.config(yscrollcommand=sb.set)
        Label(mframe3, text="Add:", **COLOR_CONFIG).pack(side='left')
        self._module_entry = Entry(mframe3, **ENTRY_CONFIG)
        self._module_entry.pack(side='left', fill='x', expand=1)
        self._module_entry.bind('<Return>', self._entry_module)
        self._module_browse = Button(mframe3, text="Browse",
                                     command=self._browse_module,
                                     **BUTTON_CONFIG) 
        self._module_browse.pack(side='right', expand=0, padx=2)
        
    def _init_progress_bar(self, mainframe):
        pframe1 = Frame(mainframe, background=BG_COLOR)
        pframe1.pack(side="bottom", fill='x', expand=0)
        self._go_button = Button(pframe1, width=4, text='Start',
                                 underline=0, command=self._go,
                                 **BUTTON_CONFIG)
        self._go_button.pack(side='left', padx=4)
        pframe2 = Frame(pframe1, relief='groove', border=2,
                        background=BG_COLOR) 
        pframe2.pack(side="top", fill='x', expand=1, padx=4, pady=3)
        Label(pframe2, text='Progress:', **COLOR_CONFIG).pack(side='left')
        H = self._H = PROGRESS_HEIGHT
        W = self._W = PROGRESS_WIDTH
        c = self._canvas = Canvas(pframe2, height=H+DH, width=W+DW, 
                                  background=PROGRESS_BG, border=0,
                                  selectborderwidth=0, relief='sunken',
                                  insertwidth=0, insertborderwidth=0,
                                  highlightbackground=BG_COLOR)
        self._canvas.pack(side='left', fill='x', expand=1, padx=4)
        self._r2 = c.create_rectangle(0,0,0,0, outline=PROGRESS_COLOR2)
        self._r3 = c.create_rectangle(0,0,0,0, outline=PROGRESS_COLOR3)
        self._r1 = c.create_rectangle(0,0,0,0, fill=PROGRESS_COLOR1,
                                      outline='')
        self._canvas.bind('<Configure>', self._configure)

    def _init_messages(self, msgsframe, ctrlframe):
        self._downImage = PhotoImage(master=self._root, data=DOWN_GIF)
        self._upImage = PhotoImage(master=self._root, data=UP_GIF)

        # Set up the messages control frame
        b1 = Button(ctrlframe, text="Messages", justify='center',
                    command=self._messages_toggle, underline=0, 
                    highlightthickness=0, activebackground=BG_COLOR, 
                    border=0, relief='flat', padx=2, pady=0, **COLOR_CONFIG) 
        b2 = Button(ctrlframe, image=self._downImage, relief='flat', 
                    border=0, command=self._messages_toggle,
                    activebackground=BG_COLOR, **COLOR_CONFIG) 
        self._message_button = b2
        self._messages_visible = 0
        b2.pack(side="left")
        b1.pack(side="left")

        f = Frame(msgsframe, background=BG_COLOR)
        f.pack(side='top', expand=1, fill='both')
        messages = Text(f, width=80, height=10, **ENTRY_CONFIG)
        messages['state'] = 'disabled'
        messages.pack(fill='both', expand=1, side='left')
        self._messages = messages

        # Add a scrollbar
        sb = Scrollbar(f, orient='vertical', **SB_CONFIG)
        sb.pack(fill='y', side='right')
        sb['command'] = messages.yview
        messages['yscrollcommand'] = sb.set

        # Set up some colorization tags
        messages.tag_config('error', foreground=ERROR_COLOR)
        messages.tag_config('warning', foreground=WARNING_COLOR)
        messages.tag_config('guierror', foreground=GUIERROR_COLOR)
        messages.tag_config('message', foreground=MESSAGE_COLOR)
        messages.tag_config('header', foreground=HEADER_COLOR)
        messages.tag_config('uline', underline=1)

        # Keep track of tag state..
        self._in_header = 0
        self._last_tag = 'error'

        # Add some buttons
        buttons = Frame(msgsframe, background=BG_COLOR)
        buttons.pack(side='bottom', fill='x')
        self._show_errors = IntVar(self._root)
        self._show_errors.set(1)
        self._show_warnings = IntVar(self._root)
        self._show_warnings.set(1)
        self._show_messages = IntVar(self._root)
        self._show_messages.set(0)
        Checkbutton(buttons, text='Show Messages', var=self._show_messages,
                    command=self._update_msg_tags, 
                    **SHOWMSG_CONFIG).pack(side='left')
        Checkbutton(buttons, text='Show Warnings', var=self._show_warnings,
                    command=self._update_msg_tags, 
                    **SHOWWRN_CONFIG).pack(side='left')
        Checkbutton(buttons, text='Show Errors', var=self._show_errors,
                    command=self._update_msg_tags, 
                    **SHOWERR_CONFIG).pack(side='left')
        self._update_msg_tags()

    def _update_msg_tags(self, *e):
        elide_errors = not self._show_errors.get()
        elide_warnings = not self._show_warnings.get()
        elide_messages = not self._show_messages.get()
        elide_headers = elide_errors and elide_warnings
        self._messages.tag_config('error', elide=elide_errors)
        self._messages.tag_config('guierror', elide=elide_errors)
        self._messages.tag_config('warning', elide=elide_warnings)
        self._messages.tag_config('message', elide=elide_messages)
        self._messages.tag_config('header', elide=elide_headers)

    def _init_options(self, optsframe, ctrlframe):
        self._leftImage=PhotoImage(master=self._root, data=LEFT_GIF)
        self._rightImage=PhotoImage(master=self._root, data=RIGHT_GIF)

        # Set up the options control frame
        b1 = Button(ctrlframe, text="Options", justify='center',
                    border=0, relief='flat',
                    command=self._options_toggle, padx=2,
                    underline=0, pady=0, highlightthickness=0,
                    activebackground=BG_COLOR, **COLOR_CONFIG) 
        b2 = Button(ctrlframe, image=self._rightImage, relief='flat', 
                    border=0, command=self._options_toggle,
                    activebackground=BG_COLOR, **COLOR_CONFIG) 
        self._option_button = b2
        self._options_visible = 0
        b2.pack(side="right")
        b1.pack(side="right")

        oframe2 = Frame(optsframe, relief='groove', border=2,
                        background=BG_COLOR)
        oframe2.pack(side="right", fill='both',
                     expand=0, padx=4, pady=3, ipadx=4)
        
        Label(oframe2, text="Project Options", font='helvetica -16',
              **COLOR_CONFIG).pack(anchor='w')
        oframe3 = Frame(oframe2, background=BG_COLOR)
        oframe3.pack(fill='x')
        oframe4 = Frame(oframe2, background=BG_COLOR)
        oframe4.pack(fill='x')
        oframe7 = Frame(oframe2, background=BG_COLOR)
        oframe7.pack(fill='x')
        div = Frame(oframe2, background=BG_COLOR, border=1, relief='sunk')
        div.pack(ipady=1, fill='x', padx=4, pady=2)

        Label(oframe2, text="Help File", font='helvetica -16',
              **COLOR_CONFIG).pack(anchor='w')
        oframe5 = Frame(oframe2, background=BG_COLOR)
        oframe5.pack(fill='x')
        div = Frame(oframe2, background=BG_COLOR, border=1, relief='sunk')
        div.pack(ipady=1, fill='x', padx=4, pady=2)

        Label(oframe2, text="CSS Stylesheets", font='helvetica -16',
              **COLOR_CONFIG).pack(anchor='w')
        oframe6 = Frame(oframe2, background=BG_COLOR)
        oframe6.pack(fill='x')

        #==================== oframe3 ====================
        # -n NAME, --name NAME
        row = 0
        l = Label(oframe3, text="Project Name:", **COLOR_CONFIG)
        l.grid(row=row, column=0, sticky='e')
        self._name_entry = Entry(oframe3, **ENTRY_CONFIG)
        self._name_entry.grid(row=row, column=1, sticky='ew', columnspan=3)

        # -u URL, --url URL
        row += 1
        l = Label(oframe3, text="Project URL:", **COLOR_CONFIG)
        l.grid(row=row, column=0, sticky='e')
        self._url_entry = Entry(oframe3, **ENTRY_CONFIG)
        self._url_entry.grid(row=row, column=1, sticky='ew', columnspan=3)

        # -o DIR, --output DIR
        row += 1
        l = Label(oframe3, text="Output Directory:", **COLOR_CONFIG)
        l.grid(row=row, column=0, sticky='e')
        self._out_entry = Entry(oframe3, **ENTRY_CONFIG)
        self._out_entry.grid(row=row, column=1, sticky='ew', columnspan=2)
        self._out_browse = Button(oframe3, text="Browse",
                                  command=self._browse_out,
                                  **BUTTON_CONFIG) 
        self._out_browse.grid(row=row, column=3, sticky='ew', padx=2)

        #==================== oframe4 ====================
        # --no-frames
        row = 0
        self._frames_var = IntVar(self._root)
        self._frames_var.set(1)
        l = Label(oframe4, text="Generate a frame-based table of contents",
                  **COLOR_CONFIG)
        l.grid(row=row, column=1, sticky='w')
        cb = Checkbutton(oframe4, var=self._frames_var, **CBUTTON_CONFIG)
        cb.grid(row=row, column=0, sticky='e')

        # --no-private
        row += 1
        self._private_var = IntVar(self._root)
        self._private_var.set(1)
        l = Label(oframe4, text="Generate documentation for private objects",
                  **COLOR_CONFIG)
        l.grid(row=row, column=1, sticky='w')
        cb = Checkbutton(oframe4, var=self._private_var, **CBUTTON_CONFIG)
        cb.grid(row=row, column=0, sticky='e')

        # --show-imports
        row += 1
        self._imports_var = IntVar(self._root)
        self._imports_var.set(0)
        l = Label(oframe4, text="List imported classes and functions",
                  **COLOR_CONFIG)
        l.grid(row=row, column=1, sticky='w')
        cb = Checkbutton(oframe4, var=self._imports_var, **CBUTTON_CONFIG)
        cb.grid(row=row, column=0, sticky='e')

        #==================== oframe7 ====================
        # --docformat
        row += 1
        l = Label(oframe7, text="Default Docformat:", **COLOR_CONFIG)
        l.grid(row=row, column=0, sticky='e')
        df_var = self._docformat_var = StringVar(self._root)
        self._docformat_var.set('epytext')
        b = Radiobutton(oframe7, var=df_var, text='Epytext',
                        value='epytext', **CBUTTON_CONFIG)
        b.grid(row=row, column=1, sticky='w')
        b = Radiobutton(oframe7, var=df_var, text='ReStructuredText',
                        value='restructuredtext', **CBUTTON_CONFIG)
        b.grid(row=row, column=2, columnspan=2, sticky='w')
        row += 1
        b = Radiobutton(oframe7, var=df_var, text='Plaintext',
                        value='plaintext', **CBUTTON_CONFIG)
        b.grid(row=row, column=1, sticky='w')
        b = Radiobutton(oframe7, var=df_var, text='Javadoc',
                        value='javadoc', **CBUTTON_CONFIG)
        b.grid(row=row, column=2, columnspan=2, sticky='w')
        row += 1

        # Separater
        Frame(oframe7, background=BG_COLOR).grid(row=row, column=1, pady=3)
        row += 1

        # --inheritance
        l = Label(oframe7, text="Inheritance Style:", **COLOR_CONFIG)
        l.grid(row=row, column=0, sticky='e')
        inh_var = self._inheritance_var = StringVar(self._root)
        self._inheritance_var.set('grouped')
        b = Radiobutton(oframe7, var=inh_var, text='Grouped',
                        value='grouped', **CBUTTON_CONFIG)
        b.grid(row=row, column=1, sticky='w')
        b = Radiobutton(oframe7, var=inh_var, text='Listed',
                        value='listed', **CBUTTON_CONFIG)
        b.grid(row=row, column=2, sticky='w')
        b = Radiobutton(oframe7, var=inh_var, text='Included',
                        value='included', **CBUTTON_CONFIG)
        b.grid(row=row, column=3, sticky='w')
        row += 1

        #==================== oframe5 ====================
        # --help-file FILE
        row = 0
        self._help_var = StringVar(self._root)
        self._help_var.set('default')
        b = Radiobutton(oframe5, var=self._help_var,
                        text='Default',
                        value='default', **CBUTTON_CONFIG)
        b.grid(row=row, column=1, sticky='w')
        row += 1
        b = Radiobutton(oframe5, var=self._help_var,
                        text='Select File',
                        value='-other-', **CBUTTON_CONFIG)
        b.grid(row=row, column=1, sticky='w')
        self._help_entry = Entry(oframe5, **ENTRY_CONFIG)
        self._help_entry.grid(row=row, column=2, sticky='ew')
        self._help_browse = Button(oframe5, text='Browse',
                                   command=self._browse_help,
                                   **BUTTON_CONFIG)
        self._help_browse.grid(row=row, column=3, sticky='ew', padx=2)
        
        from epydoc.css import STYLESHEETS
        items = STYLESHEETS.items()
        def _css_sort(css1, css2):
            if css1[0] == 'default': return -1
            elif css2[0] == 'default': return 1
            else: return cmp(css1[0], css2[0])
        items.sort(_css_sort)

        #==================== oframe6 ====================
        # -c CSS, --css CSS
        # --private-css CSS
        row = 0
        l = Label(oframe6, text="Public", **COLOR_CONFIG)
        l.grid(row=row, column=0, sticky='e')
        l = Label(oframe6, text="Private", **COLOR_CONFIG)
        l.grid(row=row, column=1, sticky='w')
        row += 1
        css_var = self._css_var = StringVar(self._root)
        css_var.set('default')
        private_css_var = self._private_css_var = StringVar(self._root)
        private_css_var.set('default')
        for (name, (sheet, descr)) in items:
            b = Radiobutton(oframe6, var=css_var, value=name, **CBUTTON_CONFIG)
            b.grid(row=row, column=0, sticky='e')
            b = Radiobutton(oframe6, var=private_css_var, value=name, 
                            text=name, **CBUTTON_CONFIG)
            b.grid(row=row, column=1, sticky='w')
            l = Label(oframe6, text=descr, **COLOR_CONFIG)
            l.grid(row=row, column=2, sticky='w')
            row += 1
        b = Radiobutton(oframe6, var=css_var, value='-other-', 
                        **CBUTTON_CONFIG)
        b.grid(row=row, column=0, sticky='e')
        b = Radiobutton(oframe6, text='Select File', var=private_css_var, 
                        value='-other-', **CBUTTON_CONFIG)
        b.grid(row=row, column=1, sticky='w')
        self._css_entry = Entry(oframe6, **ENTRY_CONFIG)
        self._css_entry.grid(row=row, column=2, sticky='ew')
        self._css_browse = Button(oframe6, text="Browse",
                                  command=self._browse_css,
                                  **BUTTON_CONFIG) 
        self._css_browse.grid(row=row, column=3, sticky='ew', padx=2)

    def _init_bindings(self):
        self._root.bind('<Delete>', self._delete_module)
        self._root.bind('<Alt-o>', self._options_toggle)
        self._root.bind('<Alt-m>', self._messages_toggle)
        self._root.bind('<F5>', self._go)
        self._root.bind('<Alt-s>', self._go)
        
        self._root.bind('<Control-n>', self._new)
        self._root.bind('<Control-o>', self._open)
        self._root.bind('<Control-s>', self._save)
        self._root.bind('<Control-a>', self._saveas)

    def _options_toggle(self, *e):
        if self._options_visible:
            self._optsframe.forget()
            self._option_button['image'] = self._rightImage
            self._options_visible = 0
        else:
            self._optsframe.pack(fill='both', side='right')
            self._option_button['image'] = self._leftImage
            self._options_visible = 1

    def _messages_toggle(self, *e):
        if self._messages_visible:
            self._msgsframe.forget()
            self._message_button['image'] = self._downImage
            self._messages_visible = 0
        else:
            self._msgsframe.pack(fill='both', side='bottom', expand=1)
            self._message_button['image'] = self._upImage
            self._messages_visible = 1

    def _configure(self, event):
        self._W = event.width-DW

    def _delete_module(self, *e):
        selection = self._module_list.curselection()
        if len(selection) != 1: return
        self._module_list.delete(selection[0])

    def _entry_module(self, *e):
        modules = [self._module_entry.get()]
        if glob.has_magic(modules[0]):
            modules = glob.glob(modules[0])
        for name in modules:
            self.add_module(name, check=1)
        self._module_entry.delete(0, 'end')

    def _browse_module(self, *e):
        title = 'Select a module for documentation'
        ftypes = [('Python module', '.py'),
                  ('Python extension', '.so'),
                  ('All files', '*')]
        filename = askopenfilename(filetypes=ftypes, title=title,
                                   defaultextension='.py')
        if not filename: return
        self.add_module(filename, check=1)
        
    def _browse_css(self, *e):
        title = 'Select a CSS stylesheet'
        ftypes = [('CSS Stylesheet', '.css'), ('All files', '*')]
        filename = askopenfilename(filetypes=ftypes, title=title,
                                   defaultextension='.css')
        if not filename: return
        self._css_entry.delete(0, 'end')
        self._css_entry.insert(0, filename)

    def _browse_help(self, *e):
        title = 'Select a help file'
        self._help_var.set('-other-')
        ftypes = [('HTML file', '.html'), ('All files', '*')]
        filename = askopenfilename(filetypes=ftypes, title=title,
                                   defaultextension='.html')
        if not filename: return
        self._help_entry.delete(0, 'end')
        self._help_entry.insert(0, filename)

    def _browse_out(self, *e):
        ftypes = [('All files', '*')]
        title = 'Choose the output directory'
        if askdirectory is not None:
            filename = askdirectory(mustexist=0, title=title)
            if not filename: return
        else:
            # Hack for Python 2.1 or earlier:
            filename = asksaveasfilename(filetypes=ftypes, title=title,
                                         initialfile='--this directory--')
            if not filename: return
            (f1, f2) = os.path.split(filename)
            if f2 == '--this directory--': filename = f1
        self._out_entry.delete(0, 'end')
        self._out_entry.insert(0, filename)

    def destroy(self, *e):
        if self._root is None: return

        # Unload any modules that we've imported
        for m in sys.modules.keys():
            if m not in self._old_modules: del sys.modules[m]
        self._root.destroy()
        self._root = None

    def add_module(self, name, check=0):
        from epydoc.imports import import_module, find_modules

        # First, expand packages.
        if os.path.isdir(name):
            module_names = find_modules(name)
            if not module_names:
                sys.stderr.write('!!!Error: %r is not a pacakge\n!!!' % name)
                self._update_messages()
                self._root.bell()
        else:
            module_names = [name]

        for name in module_names:
            # Check that it's a good module, if requested.
            if check:
                try: import_module(name)
                except ImportError, e:
                    print >>sys.stderr, e
                    self._update_messages()
                    self._root.bell()
                    continue

            # Add the module to the list of modules.
            self._module_list.insert('end', name)
            self._module_list.yview('end')
        
    def mainloop(self, *args, **kwargs):
        self._root.mainloop(*args, **kwargs)

    def _getopts(self):
        options = {}
        options['modules'] = self._module_list.get(0, 'end')
        options['prj_name'] = self._name_entry.get() or ''
        options['prj_url'] = self._url_entry.get() or None
        options['docformat'] = self._docformat_var.get()
        options['inheritance'] = self._inheritance_var.get()
        options['outdir'] = self._out_entry.get() or 'html'
        options['frames'] = self._frames_var.get()
        options['private'] = self._private_var.get()
        options['show_imports'] = self._imports_var.get()
        if self._help_var.get() == '-other-':
            options['help'] = self._help_entry.get() or None
        else:
            options['help'] = None
        if self._css_var.get() == '-other-':
            options['css'] = self._css_entry.get() or 'default'
        else:
            options['css'] = self._css_var.get() or 'default'
        if self._private_css_var.get() == '-other-':
            options['private_css'] = self._css_entry.get() or 'default'
        else:
            options['private_css'] = self._private_css_var.get() or 'default'
        return options
    
    def _go(self, *e):
        if len(self._module_list.get(0,'end')) == 0:
            self._root.bell()
            return

        if self._progress[0] != None:
            self._cancel[0] = 1
            return

        # Construct the argument list for document().
        opts = self._getopts()
        self._progress[0] = 0.0
        self._cancel[0] = 0
        args = (opts, self._progress, self._cancel)

        # Clear the messages window.
        self._messages['state'] = 'normal'
        self._messages.delete('0.0', 'end')
        self._messages['state'] = 'disabled'
        self._errpipe.clear()

        # Restore the module list.  This will force re-loading of
        # anything that we're documenting.
        for m in sys.modules.keys():
            if m not in self._old_modules:
                del sys.modules[m]

        # Reset the uid cache (otherwise, we would get UID conflicts,
        # since we just re-loaded everything).
        from epydoc.uid import reset_uid_cache
        reset_uid_cache()
    
        # Start documenting
        start_new_thread(document, args)

        # Start the progress bar.
        self._go_button['text'] = 'Stop'
        self._afterid += 1
        dt = 300 # How often to update, in milliseconds
        self._update(dt, self._afterid)

    def _update_messages(self):
        # Update messages (using the _PipeIO that we installed on
        # sys.stderr).
        while 1:
            data = self._errpipe.read()
            if data is None: break
            self._messages['state'] = 'normal'
            if data.startswith('***') and data.endswith('***'):
                self._messages.insert('end', data[3:-3], 'message')
            elif data.startswith('!!!') and data.endswith('!!!'):
                self._messages.insert('end', data[3:-3], 'guierror')
            else:
                if data == '\n':
                    if self._last_tag != 'header2':
                        self._messages.insert('end', '\n', self._last_tag)
                elif data == '='*75:
                    if self._messages.get('end-3c', 'end') == '\n\n\n':
                        self._messages.delete('end-1c')
                    self._in_header = 1
                    self._messages.insert('end', ' '*75, 'uline header')
                    self._last_tag = 'header'
                elif data == '-'*75:
                    self._in_header = 0
                    self._last_tag = 'header2'
                elif self._in_header:
                    self._messages.insert('end', data, 'header')
                    self._last_tag = 'header'
                elif re.match(r'\s*(L\d+:|-)?\s*Warning: ', data):
                    self._messages.insert('end', data, 'warning')
                    self._last_tag = 'warning'
                else:
                    self._messages.insert('end', data, 'error')
                    self._last_tag = 'error'

            self._messages['state'] = 'disabled'
            self._messages.yview('end')

    def _update(self, dt, id):
        if self._root is None: return
        if self._progress[0] is None: return
        if id != self._afterid: return

        # Update the messages box
        self._update_messages()

        # Update the progress bar.
        if self._progress[0] == 'done': p = self._W + DX
        elif self._progress[0] == 'cancel': p = -5
        else: p = DX + self._W * self._progress[0]
        self._canvas.coords(self._r1, DX+1, DY+1, p, self._H+1)
        self._canvas.coords(self._r2, DX, DY, p-1, self._H)
        self._canvas.coords(self._r3, DX+1, DY+1, p, self._H+1)

        # Are we done?
        if self._progress[0] in ('done', 'cancel'):
            if self._progress[0] == 'cancel': self._root.bell()
            self._go_button['text'] = 'Start'
            self._progress[0] = None
            return

        self._root.after(dt, self._update, dt, id)

    def _new(self, *e):
        self._module_list.delete(0, 'end')
        self._name_entry.delete(0, 'end')
        self._url_entry.delete(0, 'end')
        self._docformat_var.set('epytext')
        self._inheritance_var.set('grouped')
        self._out_entry.delete(0, 'end')
        self._module_entry.delete(0, 'end')
        self._css_entry.delete(0, 'end')
        self._help_entry.delete(0, 'end')
        self._frames_var.set(1)
        self._private_var.set(1)
        self._imports_var.set(0)
        self._css_var.set('default')
        self._private_css_var.set('default')
        self._help_var.set('default')
        self._filename = None

    def _open(self, *e):
        title = 'Open project'
        ftypes = [('Project file', '.prj'),
                  ('All files', '*')]
        filename = askopenfilename(filetypes=ftypes, title=title,
                                   defaultextension='.css')
        if not filename: return
        self.open(filename)

    def open(self, prjfile):
        from epydoc.css import STYLESHEETS
        self._filename = prjfile
        try:
            opts = load(open(prjfile, 'r'))
            
            modnames = list(opts.get('modules', []))
            modnames.sort()
            self._module_list.delete(0, 'end')
            for name in modnames:
                self.add_module(name)
            self._module_entry.delete(0, 'end')
                
            self._name_entry.delete(0, 'end')
            if opts.get('prj_name'):
                self._name_entry.insert(0, opts['prj_name'])
                
            self._url_entry.delete(0, 'end')
            if opts.get('prj_url'):
                self._url_entry.insert(0, opts['prj_url'])

            self._docformat_var.set(opts.get('docformat', 'epytext'))
            self._inheritance_var.set(opts.get('inheritance', 'grouped'))

            self._help_entry.delete(0, 'end')
            if opts.get('help') is None:
                self._help_var.set('default')
            else:
                self._help_var.set('-other-')
                self._help_entry.insert(0, opts.get('help'))
                
            self._out_entry.delete(0, 'end')
            self._out_entry.insert(0, opts.get('outdir', 'html'))

            self._frames_var.set(opts.get('frames', 1))
            self._private_var.set(opts.get('private', 1))
            self._imports_var.set(opts.get('show_imports', 0))
            
            self._css_entry.delete(0, 'end')
            if opts.get('css', 'default') in STYLESHEETS.keys():
                self._css_var.set(opts.get('css', 'default'))
            else:
                self._css_var.set('-other-')
                self._css_entry.insert(0, opts.get('css', 'default'))

            if opts.get('private_css', 'default') in STYLESHEETS.keys():
                self._private_css_var.set(opts.get('private_css', 'default'))
            else:
                self._private_css_var.set('-other-')
                self._css_entry.insert(0, opts.get('private_css', 'default'))
                                                   
        except Exception, e:
            sys.stderr.write('!!!Error opening %s: %s\n!!!' % (prjfile, e))
            self._root.bell()
        
    def _save(self, *e):
        if self._filename is None: return self._saveas()
        try:
            opts = self._getopts()
            dump(opts, open(self._filename, 'w'))
        except Exception, e:
            if self._filename is None:
                sys.stderr.write('!!!Error saving: %s\n!!!' %  e)
            else:
                sys.stderr.write('!!!Error saving %s: %s\n!!!' %
                                 (self._filename, e))
            self._root.bell()
             
    def _saveas(self, *e):
        title = 'Save project as'
        ftypes = [('Project file', '.prj'), ('All files', '*')]
        filename = asksaveasfilename(filetypes=ftypes, title=title,
                                     defaultextension='.prj')
        if not filename: return
        self._filename = filename
        self._save()

def _version():
    """
    Display the version information, and exit.
    @rtype: C{None}
    """
    import epydoc
    print "Epydoc version %s" % epydoc.__version__
    sys.exit(0)

# At some point I could add:
#   --show-messages, --hide-messages
#   --show-options, --hide-options
def _usage():
    print
    print 'Usage: epydocgui [OPTIONS] [FILE.prj | MODULES...]'
    print
    print '    FILE.prj                  An epydoc GUI project file.'
    print '    MODULES...                A list of Python modules to document.'
    print '    -V, --version             Print the version of epydoc.'
    print '    -h, -?, --help, --usage   Display this usage message'
    print '    --debug                   Do not suppress error messages'
    print
    sys.exit(0)

def _error(str):
    str = '%s; run "%s -h" for usage' % (str, os.path.basename(sys.argv[0]))
    if len(str) > 80:
        i = str.rfind(' ', 0, 80)
        if i>0: str = str[:i]+'\n'+str[i+1:]
    print >>sys.stderr, str
    sys.exit(1)
    
def gui():
    global DEBUG
    sys.stderr = sys.__stderr__
    projects = []
    modules = []
    for arg in sys.argv[1:]:
        if arg[0] == '-':
            arg = arg.lower()
            if arg in ('-h', '--help', '-?', '--usage'): _usage()
            elif arg in ('-V', '--version'): _version()
            elif arg in ('--debug',): DEBUG = 1
            else:
                _error('Unknown parameter %r' % arg)
        elif arg[-4:] == '.prj': projects.append(arg)
        else: modules.append(arg)

    if len(projects) > 1:
        _error('Too many projects')
    if len(projects) == 1:
        if len(modules) > 0:
            _error('You must specify either a project or a list of modules')
        if not os.path.exists(projects[0]):
            _error('Cannot open project file %s' % projects[0])
        gui = EpydocGUI()
        gui.open(projects[0])
        gui.mainloop()
    else:
        gui = EpydocGUI()
        for module in modules: gui.add_module(module, check=1)
        gui.mainloop()

if __name__ == '__main__': gui()

