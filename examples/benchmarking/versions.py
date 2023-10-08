from __future__ import print_function
from __future__ import unicode_literals
import os
import re
import shutil
from subprocess import Popen, PIPE
import sys
import tempfile

from fipy.tools.parser import parse

from examples.benchmarking.utils import monitor

def main():
    url = parse('--svnURL', action='store',
                  type='string', default=None)

    revision = parse('--svnRevision', action='store',
                     type='string', default="HEAD:0")

    args = sys.argv[1:]

    import pysvn

    m = re.match(r"(\d+|HEAD):(\d+|HEAD)", revision)

    def str2rev(s):
        if s == "HEAD":
            return pysvn.Revision(pysvn.opt_revision_kind.head)
        else:
            return pysvn.Revision(pysvn.opt_revision_kind.number, int(s))

    if m is None:
        revisionStart = revisionEnd = str2rev(revision)
    else:
        revisionStart = str2rev(m.group(1))
        revisionEnd = str2rev(m.group(2))

    client = pysvn.Client()

    # svn manipulations on the working copy in-place are dangerous

    if url is None:
        info = client.info('.')
        url = info.url

    dir = tempfile.mkdtemp()

    env = os.environ.copy()
    env['PYTHONPATH'] = dir

    benchmarker = os.path.join(os.path.dirname(__file__),
                               "benchmarker.py")

    pyth = sys.executable or "python"

    try:
        client.checkout(url, dir)

        scanf_e = "[-+]?(\d+(\.\d*)?|\d*\.\d+)([eE][-+]?\d+)?"

        reCPU = re.compile("cpu time: (%s) s / step / cell" % scanf_e)
        reRSZ = re.compile("max resident memory: (%s) B / cell" % scanf_e)
        reVSZ = re.compile("max virtual memory: (%s) B / cell" % scanf_e)

        for entry in client.log(url,
                                revision_start=revisionStart,
                                revision_end=revisionEnd):
            try:
                client.update(dir, revision=entry.revision)

                p = Popen([pyth, benchmarker] + args
                          + ["--numberOfSteps=0"],
                          stdout=PIPE,
                          stderr=PIPE)

                cpu0, rsz0, vsz0 = monitor(p)

                p = Popen([pyth, benchmarker] + args
                          + ["--cpuBaseLine=%f" % cpu0],
                          stdout=PIPE,
                          stderr=PIPE,
                          env=env)

                cpu, rsz, vsz = monitor(p)

                print(entry.revision.number, cpu, rsz, vsz)

            except Exception as e:
                print(entry.revision.number, e)
    except Exception as e:
        print(e)

    shutil.rmtree(dir)

if __name__ == "__main__":
    main()
