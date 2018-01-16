BootStrap: debootstrap
OSVersion: trusty
MirrorURL: http://us.archive.ubuntu.com/ubuntu/


%labels
	Maintainer daniel.maplesond@earlham.ac.uk
	Version 2.0


%runscript
	echo "This is what happens when you run the container..."

%post
	echo "Updating repos"
	apt-get update -qq

	echo "Installing C libraries"
	apt-get -y --force-yes install libstdc++6 libc6-dev build-essential dh-autoreconf

	echo "Installing misc apps"
	apt-get -y --force-yes install vim git

	echo "Installing internet access apps"
	apt-get -y --force-yes install wget php5-curl libcurl4-gnutls-dev libssl-dev

	echo "Installing compression packages"
	apt-get -y --force-yes install libzip2 zlib1g-dev libbz2-dev bzip2 liblzma-dev

	echo "Install python and packages"
	apt-get -y --force-yes install python3 python3-dev python3-setuptools sphinx-doc
	update-alternatives --install /usr/bin/python python /usr/bin/python3 1

	echo "Install boost properly"	
	wget -q http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.gz/download
	mv download boost.tar.gz
 	tar -xf boost.tar.gz
 	cd boost_1_59_0
 	./bootstrap.sh --with-libraries=chrono,timer,program_options,filesystem,system
 	./b2 -d0 install; 
 	cd ..

	echo "Installing Samtools"
	wget -q https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2
	tar -xf samtools-1.6.tar.bz2
	cd samtools-1.6
	./configure --without-curses
	make -j2
	make install
	cd ..

	echo "Build portcullis and junctools"
	rm -rf portcullis
	git clone https://github.com/maplesond/portcullis.git
	cd portcullis
	./autogen.sh
	./configure --disable-silent-rules
	make V=1 -j2
	make V=1 -j2 check
	cat tests/test-suite.log
	make install
	

%apphelp portcullis
	portcullis --help

%apprun portcullis
	portcullis "$@"


%apphelp junctools
	junctools --help

%apprun junctools
	junctools "$@"
