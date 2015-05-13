
#Portcullis



##Installation:

  - If you cloned the git repository you must first run "./autogen.sh" to create the configure and make files for your project.  Do not worry if this fails due to missing dependencies at this stage.  If you downloaded a source code distribution tarball then you can skip this step.
  - Portcullis depends on some external software: bamtools, boost, pthreads and zlib.  Please make sure these programs are correctly configured and installed on your system prior to building portcullis.
  - Now, for a typical installation on a machine where you have root access type ```./configure; make; sudo make install;```

The configure script can take several options as arguments.  One commonly modified option is ```--prefix```, which will install portcullis to a custom directory.  By default this is "/usr/local", so the portcullis executable would be found at "/usr/local/bin" by default.  In addition, some options specific to managing portcullis dependencies located in non-standard locations are:

  - ```--with-bamtools``` - for specifiying a custom bamtools installation directory
  - ```--with-boost``` - for specifiying a custom boost installation directory
  - ```--with-zlib``` - for specifying a custom zlib installation directory

Note there is not option for specifying a custom pthreads location.  We assume this is correctly installed and configured for your system already.

Type ```./configure --help``` for full details.

The Makefile for portcullis can take several goals.  Full details of common make goals can be found in the INSTALL file.  Typically, the following options can optionally used by KAT:

  - ```make check``` - runs all unit tests.  This includes unit tests for htslib and samtools which are embedded in the portcullis source tree.  To run only portcullis unit tests go into the ``tests`` subdirectory and run ``make check`` there.
  - ```make dist``` - packages the installation into a tarballed distributable.
  - ```make distcheck``` - runs some sanity tests to ensure the tarballed distributable is likely to work.


##Operating Instructions:

After portcullis has been installed, the `portcullis` executable should be available.

Typing `portcullis` or `portcullis --help` at the command line will present you with the portcullis help message.

There are 4 modes available:

    - prep    - Prepares input data so that it is suitable for junction analysis
    - junc    - Calculates junction metrics for the prepared data
    - filter  - Separates alignments based on whether they are likely to represent genuine splice junctions or not
    - bamfilt - Filters a BAM to remove any reads associated with invalid junctions
    - full    - Runs prep, junc, filter and bamfilt as a complete pipeline


##Licensing:

GNU GPL V3.  See COPYING file for more details.


##Authors:

Daniel Mapleson
Luca Venturini
David Swarbreck

See AUTHORS file for more details.


##Acknowledgements:

Affiliation: The Genome Analysis Centre (TGAC)
Funding: The Biotechnology and Biological Sciences Research Council (BBSRC)
