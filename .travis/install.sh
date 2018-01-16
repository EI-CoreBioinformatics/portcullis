#!/bin/bash


- mkdir -p exdeps
- cd exdeps

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then

	# Install boost -- looks like boost 1.55.0_2 is already installed with brew
	#brew install boost

	# Download anaconda
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;

	# Install samtools
	brew tap homebrew/science
	brew install homebrew/science/samtools

else

	# Download conda
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;

	# Samtools installation (just use an old version for now... ideally use V1.3 to be consistent with internal htslib)
	sudo apt-get install samtools
fi

bash miniconda.sh -b -p $TRAVIS_BUILD_DIR/exdeps/miniconda; export PATH="$TRAVIS_BUILD_DIR/exdeps/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda install libgcc
conda info -a
conda create -q -n test-environment python=3.5 anaconda
source activate test-environment;

cd ..
