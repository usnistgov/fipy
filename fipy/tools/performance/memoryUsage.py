"""

This python script is ripped from
http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/286222/index_txt

"""

__all__ = []

import os

_proc_status = '/proc/%d/status' % os.getpid()

_scale = {'kB': 1024.0, 'mB': 1024.0*1024.0,
          'KB': 1024.0, 'MB': 1024.0*1024.0}

def _VmB(VmKey, pid = None):
    '''Private.
    '''
    global _proc_status, _scale
    if pid is not None:
        _proc_status = '/proc/%d/status' % pid

     # get pseudo file  /proc/<pid>/status
    try:
        t = open(_proc_status)
        v = t.read()
        t.close()
    except:
        return 0.0  # non-Linux?
     # get VmKey line e.g. 'VmRSS:  9999  kB\n ...'
    i = v.index(VmKey)
    v = v[i:].split(None, 3)  # whitespace
    if len(v) < 3:
        return 0.0  # invalid format?
     # convert Vm value to bytes
    return float(v[1]) * _scale[v[2]]


def _memory(since=0.0, pid = None):
    '''Return memory usage in bytes.
    '''
    return _VmB('VmSize:', pid) - since


def _resident(since=0.0):
    '''Return resident memory usage in bytes.
    '''
    return _VmB('VmRSS:') - since


def _stacksize(since=0.0):
    '''Return stack size in bytes.
    '''
    return _VmB('VmStk:') - since

def _peak(since=0.0):
    '''Return stack size in bytes.
    '''
    return _VmB('VmPeak:') - since
