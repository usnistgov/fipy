====================
Buildbot maintenance
====================

---------------
Website address
---------------

http://build.cmi.kent.edu:8010/

---------------------------
Modifying the configuration
---------------------------

The buildbot configuration file is hosted in a Git repository. To clone the
repository::

    $ git clone ssh://buildbot@build.cmi.kent.edu/fipy-buildbot.git

This should create a directory called `fipy-buildbot`::

    $ cd fipy-buildbot/
    $ ls
    buildbot.tac  master.cfg  public_html/

Edit `master.cfg` as needed::

    $ $EDITOR master.cfg

If possible, check the configuration locally::

    $ buildbot checkconfig

Push the updated config back to the repository::

    $ git add master.cfg
    $ git commit -m 'fixed flux capacitor'
    $ git push origin master

Buildbot will automatically reconfigure when the push is received.

If the Buildbot website is down after pushing, that probably means that
your change introduced an error into the config file. To debug, check the
log files at http://build.cmi.kent.edu/fipy-buildbot.log.

-------------------
Adding a buildslave
-------------------

Choose a name and a password for the new buildslave. Modify `master.cfg` as per
the instructions above (on any machine) to recognize the new buildslave by
adding a new element to the `c['slaves']` list like so::

    newslave = BuildSlave('newSlaveName', 'newSlavePasswd')

    c['slaves'] = [
                   ..., # existing slaves
                   newslave,
                  ]

Now, define one or more builders on this slave. In most cases, you will
want to define `quick` and `full` builders with default tests, which can be
done with just::

    c['builders'].extend(make_builders(name="SlaveType",
                                       slaves=[newslave]))

where `"SlaveType"` should be a generic description of the slave, e.g., it's 
OS, rather than the slave's name. 

If you need to customize the tests (perhaps because the slave is missing 
some prerequisites for the standard tests), define customized versions of 
:class:`QuickFactory` and :class:`FullFactory` and then call::

    c['builders'].extend(make_builders(name="SlaveType",
                                       slaves=[newslave],
                                       quick_factory_class=MyQuickFactory,
                                       full_factory_class=MyFullFactory))

A single, arbitrary builder can be defined by::

    c['builders'].append(make_builder(name="SlaveType",
                                      slaves=[newslave],
                                      category="BuildCategory",
                                      factory=buildFactory))

`"BuildCategory"` should be be a brief description of what type of build
this is ("quick", "full", "docs", etc.). `buildFactory` should be a
:class:`~buildbot.process.factory.BuildFactory` instance that defines the
:class:`~buildbot.process.buildsteps.BuildStep`\s to take.

Push `master.cfg` and make sure buildbot reconfigures successfully::
    
    $ git add master.cfg
    $ git commit -m "added buildslave 'newSlaveName'"
    $ git push origin master

Ensure that the buildslave can communicate on port 9989. This may entail
forwarding ports on the buildslave's router. The buildslave will receive
commands over this port from the master to begin a build and will relay results
to the master over the same port.

Ensure that `buildslave` is installed on the slave::

    slave$ buildslave --version

If not, get it::

    slave$ easy_install buildbot-slave

Then, set up the buildslave with all of FiPy's dependencies.

Ask `buildslave` to create the slave::

    slave$ buildslave create-slave [newDir] build.cmi.kent.edu:9989 \
                                   [slaveName] [slavePasswd]
                                   
Modify admin and host information::

    slave$ uname -a > newDir/info/host
    slave$ $EDITOR newDir/info/admin
    slave$ $EDITOR newDir/info/host

Start the slave::

    slave$ cd newDir
    slave$ buildslave start

Add buildslave command to the crontab::

    slave$ crontab -e
    # m  h dom mon dow command
    */10 * *   *   *   cd .../newDir ; .../buildslave start > /dev/null 2>&1 &

Py3k
====

Setting up a Python 3 build slave is a little tricky because buildbot
doesn't install under Python 3 (at least [http://twistedmatrix.com twistd]
doesn't).

What worked for me was to set up a `virtualenv` using Python 2.x to 
install `buildbot-slave` (as well as everything needed for a Python 2.x 
FiPy tester). Then I set up another `virtualenv` using `python3` and 
installed :term:`FiPy`\'s Py3k prerequisites there (`numpy` and `scipy` for now). 
Then I added `${Python2xVirtualEnv}/bin:${Python3xVirtualEnv}/bin` to the 
$PATH used by `buildslave`.

Mac OS X
========

I used [http://matforge.org/fipy/wiki/InstallFiPy/MacOSX/HomeBrew Homebrew]
and [http://matforge.org/fipy/wiki/InstallFiPy/PipInstallsPython pip] to
set up my test `virtualenv`\s and then I basically followed
http://trac.buildbot.net/wiki/UsingLaunchd, with the modifications (see
Py3k_)::

    <dict>
      <key>PATH</key>
      <string>/Applications/Gmsh.app/Contents/MacOS:/Users/fipy/.virtualenvs/buildbot/bin:/Users/fipy/.virtualenvs/buildb\
    ot3k/bin:/Users/fipy/.homebrew/share/python:/Users/fipy/.homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr\
    /X11/bin</string>
      <key>DYLD_LIBRARY_PATH</key>
      <string>/Users/fipy/.virtualenvs/buildbot/lib</string>
      </dict>

and::

    <key>UserName</key>
    <string>fipy</string>
    <key>WorkingDirectory</key>
    <string>/Users/fipy/.virtualenvs/buildbot/build/fipy-buildbot/fipy-build</string>
    <key>Nice</key>
    <integer>20</integer>

---------------
Adding a branch
---------------

Modify `master.cfg` to include the name of the branch in the `branches` list at
the top of the file like so::

    branches = [
                'trunk',
                ...,
                'branches/newBranch',
               ]

Push `master.cfg` and ensure that buildbot reconfigured correctly by checking
the website.

You can also test a branch without adding it (or test before committing) 
by using `buildbot try`.

.. attention::
   `buildbot try` is broken in Buildbot v0.8.3p1.  
   `build.cmi.kent.edu` needs to be upgraded.

