#!/bin/sh -e

# Ensure the boost submodule is present
git submodule sync
git submodule update --recursive --init deps/boost

# Build boost
cd deps/boost
./bootstrap.sh --prefix=build --with-libraries=chrono,exception,program_options,timer,filesystem,system,stacktrace
./b2 link=static install

cd ../..
