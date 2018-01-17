#!/bin/sh -e

# Ensure the boost submodule is present
git submodule sync
git submodule update --recursive --init deps/boost

# Build boost
cd deps/boost
./bootstrap.sh --prefix=../boost_build --with-libraries=chrono,exception,program_options,timer,filesystem,system,stacktrace
./b2 --prefix=../boost_build headers
./b2 --prefix=../boost_build install

cd ../..
