
#Portculis



##Installation:

  - If you cloned the git repository you must first run "./autogen.sh" to create the configure and make files for your project.  Do not worry if this fails due to missing dependencies at this stage.  If you downloaded a source code distribution tarball then you can skip this step.
  - Portculis depends on samtools.  Instead of having to install this separately we provide the source code bundled and customised for portculis here.  To compile samtools change into the samtools directory and simply type ```make```
  - Now, for a typical installation on a machine where you have root access type ```./configure; make; sudo make install;```

The configure script can take several options as arguments.  One commonly modified option is ```--prefix```, which will install portculis to a custom directory.  By default this is "/usr/local", so the portculis executable would be found at "/usr/local/bin" by default.  In addition, some options specific to managing KAT dependencies located in non-standard locations are:

  - ```--with-doxygen``` - for specifying a custom doxygen directory (doxygen is only required for generating code documention.

Type ```./configure --help``` for full details.

The Makefile for portculis can take several goals.  Full details of common make goals can be found in the INSTALL file.  Typically, the following options can optionally used by KAT:

  - ```make check``` - runs unit tests.  Requires boost unit test framework to be installed and available.
  - ```make dist``` - packages the installation into a tarballed distributable.
  - ```make distcheck``` - runs some sanity tests to ensure the tarballed distributable is likely to work.


##Operating Instructions:

After portculis has been installed, the following tools should be available:


##Licensing:

GNU GPL V3.  See COPYING file for more details.


##Authors:

Daniel Mapleson
David Swarbreck

See AUTHORS file for more details.


##Acknowledgements:

Affiliation: The Genome Analysis Centre (TGAC)
Funding: The Biotechnology and Biological Sciences Research Council (BBSRC)
