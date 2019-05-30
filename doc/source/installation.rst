.. _installation:

Installation
============

We support multiple methods for installing and running portcullis.  Hopefully your favourite container or package manager is supported below.  If not let us know and we'll try to work to get it integrated there.


Docker
~~~~~~

.. code-block:: bash
  
  # Keep in mind you need to mount in any working directories to the container with the `-v` option.
  # Ideally, mount these into the /data directory which is the container's working directory.
  docker run --it --rm -v /abspath/to/data/on/host:/data maplesond/portcullis:stable portcullis --help


Singularity
~~~~~~~~~~~

.. code-block:: bash

  # First download the container:
  singularity pull --name portcullis.img shub://maplesond/portcullis:master
  
  # Then to execute commands in the container:
  singularity exec portcullis.img portcullis --help


Conda
~~~~~

``conda install portcullis --channel=bioconda``


From brew
~~~~~~~~~

``brew install brewsci/bio/portcullis``


From source
~~~~~~~~~~~


Before installing Portcullis from source please first confirm these dependencies are installed and configured:

 - **GCC** V4.8+
 - **autoconf** V2.53+
 - **automake** V1.11+
 - **make**
 - **libtool** V2.4.2+
 - **zlib-dev**
 - **pthreads**
 - **boost-dev** V1.52+
 - **samtools** V1.2+
 - **Python3** V3.5+ (Make sure the following packages are installed: *pandas*, *matplotlib*, *setuptools*, *sphinx*, *tabulate*)

.. code-block:: bash

  # Clone the repo
  git clone git@github.com:maplesond/portcullis.git
  
  # Move into repo directory
  cd portcullis
  
  # Generate configure script
  ./autogen.sh
  
  # Confirm dependencies and generate makefiles
  # Adding --prefix <dir> will tell make install to put everything in a
  # particular directory.  Default is /usr/local.
  ./configure
  
  # Compile (increasing -j will make it go faster!
  make -j 2
  
  # Run some unit tests (you can increase -j here too)
  make -j 2 check
  
  # Install to prefix dir
  make install


*Common problems*

 - Many system python installations do not come with the C API immediately available, which prevents Portcullis from embedding python code.  We typically would recommend installing anaconda3 as this would include the latest version of python, all required python packages as well as the C API.  If you are running a debian system and the C libraries are not available by default and you wish to use the system python installation the you can install them using: ``sudo apt-get install python-dev``.  Also, if you have installed python to a custom location please verify that the *bin* directors on the *PATH* environment variable, and the lib (or lib64) directory is on the *LD_LIBRARY_PATH* or *LD_RUN_PATH* as appropriate.
 - If Portcullis is failing at the ```./autogen.sh``` step you will likely need to install autotools.  The following command should do this on MacOS: ```brew install autoconf automake libtool```.  On a debian system this can be done with: ```sudo apt-get install autoconf automake libtool```.



Internal Dependencies
---------------------

Portcullis contains *HTSlib* and *Ranger* (a random forest implementation)  in the source tree.  The user does
not need to do anything special to handle *htslib* and *ranger* as these are automatically
built and managed inside portcullis.

Portcullis also comes with a python package for analysing, comparing and converting junction files, called junctools.  This stands alone from portcullis so is not strictly required.  Should you not wish to install this you can add the ``--disable-py-install`` option to the ``configure`` script.  You can manually install this by going into the ``./scripts/junctools`` directory and typing ``python3 setup.py install``.  For more information about junctools see `junctools <junctools.html>`_ for more information.  Please note however that the portcullis python package is required for the filtering stage of Portcullis to run successfully.
