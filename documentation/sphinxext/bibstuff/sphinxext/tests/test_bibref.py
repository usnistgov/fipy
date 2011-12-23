""" Testing sphinxext module
"""

import os
from os.path import abspath, dirname, join as pjoin
import sys
from subprocess import check_call

_my_path = dirname(__file__)
doc_path = abspath(pjoin(_my_path, 'doc'))
# path to extension
sys.path.insert(0, abspath(pjoin(_my_path, '..')))
# path to custom style
sys.path.insert(0, abspath(_my_path))

import bibref as bs

from nose.tools import assert_true, assert_equal, assert_raises


def test_custom_styles():
    cs = bs.custom_styles({})
    assert_equal(cs, {})
    assert_raises(ImportError, bs.custom_styles, {'something':
                                                  'implausible.something'})
    cs = bs.custom_styles({'test-style': 'style1'})
    assert_true(hasattr(cs['test-style'], 'CitationManager'))
    # Test nested package.module
    cs = bs.custom_styles({'test-style': 'tests.style1'})
    assert_true(hasattr(cs['test-style'], 'CitationManager'))


def test_builds():
    # just does dumb builds to check they work without error
    pwd = os.getcwd()
    try:
        os.chdir(doc_path)
        check_call('make clean', shell=True)
        check_call('make html', shell=True)
        check_call('make clean', shell=True)
        # might not work on windows or systems without latex support
        check_call('make latex', shell=True)
    finally:
        os.chdir(pwd)
