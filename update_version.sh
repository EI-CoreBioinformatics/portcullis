#!/bin/bash

VERSION=$1

MAJOR="$(cut -d '.' -f 1 <<< "$VERSION")"
MINOR="$(cut -d '.' -f 2 <<< "$VERSION")"

# Configure.ac
sed -i "s/AC_INIT(\[portcullis\],\[.*\]/AC_INIT(\[portcullis\],\[$VERSION\]/" configure.ac

# Docs
sed -i "s/version = .*/version = '$MAJOR.$MINOR'/" doc/source/conf.py
sed -i "s/release = .*/release = '$VERSION'/" doc/source/conf.py

# Singularity
if [ -d Singularity ]; then sed -i "s/Version.*/Version $VERSION/" Singularity; fi

# Docker
if [ -d Dockerfile ]; then sed -i "s/VERSION:.*/VERSION:-$VERSION\}/" Dockerfile; fi

