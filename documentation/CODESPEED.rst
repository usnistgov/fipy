=====================
Codespeed maintenance
=====================

---------------
Website address
---------------

https://github.com/tobami/codespeed

---------------------------
Modifying the configuration
---------------------------

Using Codespeed 0.8.1

To download it, go to https://github.com/tobami/codespeed and use the
downloads button on the right side to get the latest stable
version. Below that on the page is the README document which does a
pretty good job at letting you know how to get started. Once
the package is downloaded, you'll need to cd into the speedcenter/
directory and do three commands::

    $ python manage.py syncdb

In this command they'll prompt you to make a superuser.  Take this
opportunity because I don't think you'll get another one (at least, I
couldn't find one after typing no)::

    $ python manage.py migrate

After this step, you'll have created your database where the data will
be stored data.db.

Once this is setup you'll be ready to run the Codespeed server.  For
testing, it's just the localhost 8000::

    $ python manage.py runserver 8000

You can't view any of the pages until you go to /admin and create an
environment, and then you still can't really access the meat of the
Codespeed pages until there's some actual data to display.

You can manually create a piece of data in the admin page by creating
a project, executable, result, revision, benchmark, and branch (must
be named default).  Once all of that is created, you can view the data
on the Changes page and see the data as a plot on the Timeline view.

The benchmarks we've chosen for our benchmarking are:
    1. Initialization time, which times from the beginning through the
       creation of mesh and the other variables up to the first loop.
    2. First timestep; the first iteration can take significantly longer than the rest,
       so this needs to be measured separately from the rest
    3. Average of the rest of the timesteps; be advised, this benchmark is misleading in that
       the way that the benchmarks are measured is such that it may actually be
       'hitting the stopwatch' multiple times per loop iteration,
       so while the average may be, for example, 5 seconds,
       it may turn out that each iteration actually takes 20 seconds.
    4. Total runtime of example.

Codespeed uploads this data using urllib and urllib2 modules and a
dictionary with the data.  The format of the dictionary is key to
properly uploading information.  I've found that one of the most
sensitive pieces is the date.  In general, the format must be the same
as given in datetime.today().

Continuing on uploading data, the makers of the program tried to make
Codespeed user friendly, but there are some things that need to be
done that won't happen automatically.  Here is an example dictionary
for uploading data::

     data = {
     	  'commitid': '14',
    	  'branch': 'default',#Always use default for trunk/master/tip
    	  'project': 'Test',
	  'executable': 'myexe O3 64bits',
          'benchmark': 'float',
	  'environment': "Test",
          'result_value': 4000,
	  }

	  # Optional fields
	  data.update({
    	  	'revision_date': None, # Optional. Default is taken either
                           # from VCS integration or from current date
    		'result_date': current_date, # Optional, default is current date
    		'std_dev': 1.11111, # Optional. Default is blank
    		'max': 4001.6, # Optional. Default is blank
    		'min': 3995.1, # Optional. Default is blank
	  })

As I've mentioned, the date has to be in the form of current_date.
Luckily, commitid, branch, project, executable, and benchmark will all
be created for you if they do not already exist.  However, environment
must already exist or there will be an error message. Set environment to
whatever computer is running the test. All of these can
be created in the admin page of Codespeed.

Codespeed is quite customizable but all of that is in Java, so beyond
adding the FiPy logo and changing the titles of pages, we may not want
to do much customization at all.  Refer to the README for instructions
for customization.
