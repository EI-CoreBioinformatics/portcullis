![alt text](doc/source/images/portcullis_logo.png "Portcullis")

Portcullis
==========

Portcullis stands for PORTable CULLing of Invalid Splice junctions from pre-aligned
RNA-seq data.  It is known that RNAseq mapping tools generate many invalid junction
predictions, particularly in deep datasets with high coverage over splice sites.
In order to address this, instead for creating a new RNAseq mapper, with a focus
on SJ accuracy we created a tool that takes in a BAM
file generated by an RNAseq mapper of the user's own choice (e.g. Tophat2, Gsnap,
STAR2 or HISAT2) as input (i.e. it's portable).  It then, analyses and quantifies
all splice junctions in the file before, filtering (culling) those which are unlikely
to be genuine.  Portcullis output's junctions in a variety of formats making it
suitable for downstream analysis (such as differential splicing analysis and gene
modelling) without additional work.  Portcullis can also filter the original BAM file removing alignments
associated with `bad` junctions.  Both the filtered junctions and BAM files are cleaner
and more usable resources which can more effectively be used to assist in downstream
analyses such as gene prediction and genome annotation.

Installation
------------

**From package manager**

The simplist way to install portcullis is via package manager.  We support both brew and bioconda.  We try to keep these reciepes as up to date as possible but sometimes there can be a lag between updating the source code and ensuring the package managers have the latest reciepes.

For brew type: ```brew install brewsci/bio/portcullis```

For bioconda type: ```conda install portcullis```



**From source**

Installing from source will ensure you have the latest version of the software.  If you wish to install from source please first confirm that first you have these dependencies are installed and configured:

 - **GCC** V4.8+
 - **autoconf** V2.53+
 - **automake** V1.11+
 - **make**
 - **libtool** V2.4.2+
 - **zlib**
 - **pthreads**
 - **samtools** V1.2+
 - **Python3** V3.5+ (including python3 development libraries and the *pandas*, *numpy*, *scipy*, *matplotlib*, and *sklearn* packages)
 - **Sphinx-doc** V1.3+ (Optional: only required for building the documentation.)

With regards to python3 and sphinx we recommend installing anaconda3 as this contains all packages and programs required by portcullis.  Also, if you have installed python to a custom location please verify that the *bin* directors on the *PATH* environment variable.

Then proceed with the following steps:

 - Clone the git repository (For ssh: ```git clone git@github.com:maplesond/portcullis.git```; or for https: ```git clone https://github.com/maplesond/portcullis.git```), into a directory on your machine.
 - "cd" into root directory of the installation  
 - Build boost by tying ```./build_boost.sh```.
 - Create configuration script by typing: ```./autogen.sh```.
 - Generate makefiles and confirm dependencies: ```./configure```
 - Compile software: ```make```
 - Run tests (optional) ```make check```
 - Install: ```sudo make install```

The configure script can take several options as arguments.  One commonly modified
option is ```--prefix```, which will install portcullis to a custom directory.  By
default this is "/usr/local", so the portcullis executable would be found at "/usr/local/bin"
by default.  Type ```./configure --help``` for full details and available options.

NOTE: if KAT is failing at the ```./autogen.sh``` step you will likely need to install autotools.  The following command should do this on MacOS: ```brew install autoconf automake libtool```.  On a debian system this can be done with: ```sudo apt-get install autoconf automake libtool```.



Operating Instructions
----------------------

After portcullis has been installed, the ```portcullis``` executable should be available.

Typing ```portcullis``` or ```portcullis --help``` at the command line will present you with the portcullis help message.

These modes are available:

 - **prep**    - Prepares input data so that it is suitable for junction analysis
 - **junc**    - Calculates junction metrics for the prepared data
 - **filter**  - Separates alignments based on whether they are likely to represent genuine splice junctions or not
 - **bamfilt** - Filters a BAM to remove any reads associated with invalid junctions
 - **full**    - Runs prep, junc, filter and optionally bamfilt as a complete pipeline

Typing ```portcullis <mode> --help``` will bring up help and usage information specific to that mode.

In addition to portcullis, we provide a tool-suite for manipulating junction files

An online version of the manual can be found here: [https://portcullis.readthedocs.org/en/latest/](https://portcullis.readthedocs.org/en/latest/).

Licensing
---------

GNU GPL V3.  See COPYING file for more details.


Authors
-------

 * Daniel Mapleson
 * Luca Venturini
 * David Swarbreck

See AUTHORS file for more details.


Acknowledgements
----------------

Affiliation: The Genome Analysis Centre (TGAC)
Funding: The Biotechnology and Biological Sciences Research Council (BBSRC)
