#!/bin/sh -e

# Using BCP
#bcp algorithm build chrono exception filesystem phoenix program_options property_tree spirit stacktrace system timer deps/boost


# Build boost
cd deps/boost
./bootstrap.sh --prefix=build --with-libraries=chrono,exception,program_options,timer,filesystem,system,stacktrace
./b2 headers
./b2 install
cd ../..
