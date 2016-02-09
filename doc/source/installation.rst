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

Portcullis depends on some external software:
 * boost
 * samtools
 * dlib
 * pthreads
 * zlib
 * sphinx (optional)

Please make sure these programs are correctly configured and installed 
on your system prior to building portcullis.  Consult the each program's installation
guide separately for instructions on how to do this.  Should you install these dependencies
into non-standard locations you can direct portcullis to them by using the following
options when running the configure script.

  - ```--with-boost``` - for specifying a custom boost installation directory
  - ```--with-dlib``` - for specifying a custom dlib installation directory
  - ```--with-zlib``` - for specifying a custom zlib installation directory

Note there is not option for specifying a custom pthreads or samtools location.  
We assume pthreads is correctly installed and configured for your system already.  In most cases
it will be.  For samtools, we just require the executable to be on the path.

If the user has sphinx installed then documentation will also be built along with
the software.  If sphinx is not detected then the documentation building stage is
skipped and documentation won't be available locally, although it can still be 
found at: https://kat.readthedocs.org/en/latest/

Boost is statically linked and doesn't need to be available at runtime.  zlib and pthreads are 
dynamically linked so will need to be on your LD_LIBRARY_PATH,
or in one of the automatically searched lib directories in order for portcullis 
to dynamically link them at runtime.  No other non-system libraries need linking at runtime.


Internal Dependencies
---------------------

Portcullis contains HTSlib in the source tree.  The user does
not need to do anything special to handle htslib as it is automatically
built and managed inside portcullis.


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
type  ```make check``` to runs some tests.  This includes unit tests for HTSlib, 
which are embedded in the portcullis source tree.  To run only portcullis 
unit tests go into the ``tests`` subdirectory and run ```make check``` there.

Finally to install the compiled code to the specified (or default) installation
directory type ```make install```.
