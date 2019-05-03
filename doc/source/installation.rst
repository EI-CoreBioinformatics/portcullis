.. _installation:

Installation
============

From container
~~~~~~~~~~~~~~

We support both docker and singularity containers and portcullis should be accessible directly from docker hub or singularity hub.

For singularity: ``singularity pull --name portcullis.img shub://maplesond/portcullis:master``.  You can then run portcullis or junctools from within the singularity image.  e.g.: ``singularity exec portcullis.img portcullis --help``.

To run from docker container, keep in mind you need to mount in any working directories to the container with the `-v` option.  Ideally, mount these into the /data directory which is the container's working directory.  For example:

``docker run --it --rm -v /abspath/to/data/on/host:/data maplesond/portcullis:stable portcullis full genome.fa rnaalignments.bam``


From brew
~~~~~~~~~

If you have brew installed on your system you should be able to install a recent version of Portcullis by simply typing:

``brew install brewsci/bio/portcullis``


From bioconda
~~~~~~~~~~~~~

If you use bioconda you can install Portcullis using :

``conda install portcullis --channel=bioconda``


From source
~~~~~~~~~~~


Before installing Portcullis from source please first confirm these dependencies are installed and configured:

 - **GCC** V4.8+
 - **autoconf** V2.53+
 - **automake** V1.11+
 - **make**
 - **libtool** V2.4.2+
 - **zlib**
 - **pthreads**
 - **boost** V1.52+
 - **samtools** V1.2+
 - **Python3** V3.5+ (including python3 development libraries and the *pandas*, *numpy*, *tabulate* packages)
 - **Sphinx-doc** V1.3+ (Optional: only required for building the documentation.)

NOTE ON INSTALLING PYTHON: Many system python installations do not come with the C API immediately available, which prevents Portcullis from embedding python code.  We typically would recommend installing anaconda3 as this would include the latest version of python, all required python packages as well as the C API.  If you are running a debian system and the C libraries are not available by default and you wish to use the system python installation the you can install them using: ``sudo apt-get install python-dev``. Also, if you have installed python to a custom location please verify that the *bin* directors on the *PATH* environment variable, and the lib (or lib64) directory is on the *LD_LIBRARY_PATH* or *LD_RUN_PATH* as appropriate.

Then proceed with the following steps:

 - Clone the git repository (For ssh: ```git clone git@github.com:maplesond/portcullis.git```; or for https: ```git clone https://github.com/maplesond/portcullis.git```), into a directory on your machine.
 - "cd" into root directory of the installation
 - Create configuration script by typing: ```./autogen.sh```.
 - Generate makefiles and confirm dependencies: ```./configure```
 - Compile software: ```make```
 - Run tests (optional) ```make check```
 - Install: ```sudo make install```

The configure script can take several options as arguments.  One commonly modified
option is ```--prefix```, which will install portcullis to a custom directory.  By
default this is "/usr/local", so the portcullis executable would be found at "/usr/local/bin"
by default.  Type ```./configure --help``` for full details and available options.

NOTE: if Portcullis is failing at the ```./autogen.sh``` step you will likely need to install autotools.  The following command should do this on MacOS: ```brew install autoconf automake libtool```.  On a debian system this can be done with: ```sudo apt-get install autoconf automake libtool```.


Internal Dependencies
---------------------

Portcullis contains *HTSlib* and *Ranger* (a random forest implementation)  in the source tree.  The user does
not need to do anything special to handle *htslib* and *ranger* as these are automatically
built and managed inside portcullis.

Portcullis also comes with a python package for analysing, comparing and converting junction files, called junctools.  This stands alone from portcullis so is not strictly required.  Should you not wish to install this you can add the ``--disable-py-install`` option to the ``configure`` script.  You can manually install this by going into the ``./scripts/junctools`` directory and typing ``python3 setup.py install``.  For more information about junctools see `junctools <junctools.html>`_ for more information.  Please note however that the portcullis python package is required for the filtering stage of Portcullis to run successfully.
