.. _installation:

Installation
============

Portcullis is primarily a C++ application with some python scripts.  We use the 
GNU build system ``Autotools`` to assist with package management and to make the 
software portable across UNIX type operating systems.  Installation of portcullis
therefore follows a similar incantation to other autotools based projects::

  ```./configure; make; sudo make install;```

However, there are a few caveats.  If you cloned the software directly from the 
git repository you must first run "./autogen.sh" to create the configure and make 
files for your project.  If you downloaded a source code distribution tarball those
scripts are already present so you can skip this step.

External Dependencies
---------------------

Portcullis depends on some external software, specifically boost, pthreads 
and zlib.  Please make sure these programs are correctly configured and installed 
on your system prior to building portcullis.  Consult the each program's installation
guide separately for instructions on how to do this.  Should you install these dependencies
into non-standard locations you can direct portcullis to them by using the following
options when running the configure script.

  - ```--with-boost``` - for specifiying a custom boost installation directory
  - ```--with-zlib``` - for specifying a custom zlib installation directory

Note there is not option for specifying a custom pthreads location.  We assume 
this is correctly installed and configured for your system already.  In most cases
it will be.

Boost is statically linked and doesn't need to be available at runtime.  zlib is 
still dynamically linked as V0.8.2 so will need to be on your LD_LIBRARY_PATH,
or in one of the automatically searched lib directories in order for portcullis 
to dynamically link them at runtime.  We plan to make changes in the future for
zlib to be statically, linked removing the final runtime dependency.


Internal Dependencies
---------------------

Portcullis contains HTSlib, samtools and dlib in the source tree.  The user does
not need to do anything specical to handle these dependencies as they are automatically
built and managed inside portcullis.  However, it is important to note that portcullis
will create the samtools executable in it's bin directory after installation, which
may conflict with your own samtools executable if it was already installed on your
PATH.  If you do not want portcullis to potentially override or conflict with an 
existing samtools installation you might want to consider installing portcullis 
to a custom location.  You can do this with the ```--prefix``` option when 
running the configure script.  We might revisit this in the future to remove
this potential issue.


Compilation and Installation
----------------------------

First change into the portcullis root directory and run ```./configure```, providing
any options you feel are appropriate.  By default the installation directory is "/usr/local", 
so the portcullis executable would be found at "/usr/local/bin" by default.  If you
want to change this use the ``--prefix`` option as previously described.  For a full
list of configuration options type ```./configure --help```.

Next compile the software.  This can be done by typing ```make```.  The compiles
all internal dependencies and portcullis itself.

To check the code compiled correct and is operating as expected you can optionally
type  ```make check``` to runs some tests.  This includes unit tests for HTSlib 
and samtools which are embedded in the portcullis source tree.  To run only portcullis 
unit tests go into the ``tests`` subdirectory and run ```make check``` there.

Finally to install the compiled code to the specified (or default) installation
directory type ```make install```.
