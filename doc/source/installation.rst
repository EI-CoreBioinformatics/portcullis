.. _installation:

Installation
============

Before installing portcullis please first confirm these dependencies are installed and configured:

 - **GCC** V4.8+
 - **autoconf** V2.53+
 - **automake** V1.11+
 - **make**
 - **libtool** V2.4.2+
 - **zlib**
 - **pthreads**
 - **samtools** V1.2+
 - **Python3** V3.5+ (including python3 development libraries and the *numpy*, *scipy*, *matplotlib*, and *sklearn* packages)
 - **Sphinx-doc** V1.3+ (Optional: only required for building the documentation.)

With regards to python3 and sphinx we recommend installing anaconda3 as this contains all packages and programs required by portcullis.
If you have installed python to a custom location please verify that the *bin* and *lib* directories are on python_full_ver
*PATH* and *LD_LIBRARY_PATH* environment variables respectively.

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


Internal Dependencies
---------------------

Portcullis contains *HTSlib* and *Ranger* (a random forest implementation)  in the source tree.  The user does
not need to do anything special to handle *htslib* and *ranger* as these are automatically
built and managed inside portcullis.

Portcullis also comes with a python package for analysing, comparing and converting
junction files, called junctools.  Should you not wish to build / install this
you can add the ``--disable-junctools`` option to the ``configure`` script.  For more
information about junctools see `junctools <junctools.html>`_ for more information.
