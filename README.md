
#Portcullis

Portcullis predicts junctions from aligned RNA-seq data.  We expect
the user to have already generated a BAM file using a splice aware aligner of their
choice.  For example, Tophat, Gsnap, STAR or HISAT will work fine.  Portcullis is designed to
be as portable as possible so where possible does not rely on esoteric SAM tags and other
artifacts that are not consistently present in all SAM/BAMs.  Portcullis
will then analyse the BAM file to look for alignments containing gaps (REFSKIP 'N' cigar ops)
and create a detailed analysis of all distinct gaps found in the BAM file, these
are considered as potential junctions.  Portcullis provides various means (both 
manual and/or automatic) of filtering these potential junctions in order to remove 
false positives.  Portcullis can also filter the original BAM file removing alignments 
associated with `bad` junctions.  Both the filtered junctions and BAM files are cleaner
and more usable resources which can more effectively be used to assist in downstream 
analyses such as gene prediction and genome annotation. 

##Installation:

  - If you cloned the git repository you must first run ```./autogen.sh``` to create the configure and make files for your project, this step requires automake and autoconf installed.  If you downloaded a source code distribution tarball then you can skip this step.
  - Portcullis depends on some external software: samtools, boost and zlib.  Please make sure these programs are correctly configured and installed on your system prior to building portcullis.
  - Optionally install python3 and sphinx to generate documentation and use supplementary scripts.
  - Now, for a typical installation on a machine where you have root access type ```./configure; make; sudo make install;```

The configure script can take several options as arguments.  One commonly modified option is ```--prefix```, which will install portcullis to a custom directory.  By default this is "/usr/local", so the portcullis executable would be found at "/usr/local/bin" by default.  In addition, some options specific to managing portcullis dependencies located in non-standard locations are:

  - ```--with-boost``` - for specifying a custom boost installation directory
  - ```--with-zlib``` - for specifying a custom zlib installation directory

Type ```./configure --help``` for full details.

The Makefile for portcullis can take several goals.  Full details of common make goals can be found in the INSTALL file.  Typically, the following options can optionally be used:

  - ```make check``` - runs all unit tests.  This includes unit tests for htslib and samtools which are embedded in the portcullis source tree.  To run only portcullis unit tests go into the ``tests`` subdirectory and run ``make check`` there.
  - ```make pdf``` - generates a PDF copy of the manual (requires sphinx to be installed)
  

##Operating Instructions:

After portcullis has been installed, the ```portcullis``` executable should be available.

Typing ```portcullis``` or ```portcullis --help``` at the command line will present you with the portcullis help message.

There are 4 modes available:

    - prep    - Prepares input data so that it is suitable for junction analysis
    - junc    - Calculates junction metrics for the prepared data
    - filter  - Separates alignments based on whether they are likely to represent genuine splice junctions or not
    - bamfilt - Filters a BAM to remove any reads associated with invalid junctions
    - full    - Runs prep, junc, filter and bamfilt as a complete pipeline

Typing ```portcullis <mode> --help``` will bring up help and usage information specific to that mode.

##Licensing:

GNU GPL V3.  See COPYING file for more details.


##Authors:

 * Daniel Mapleson
 * Luca Venturini
 * David Swarbreck

See AUTHORS file for more details.


##Acknowledgements:

Affiliation: The Genome Analysis Centre (TGAC)
Funding: The Biotechnology and Biological Sciences Research Council (BBSRC)
