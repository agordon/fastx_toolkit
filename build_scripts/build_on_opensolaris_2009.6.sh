#!/bin/sh

export CC=gcc-4.3.2
export CXX=g++-4.3.2

#Tell pkg-config to look for libraries in /usr/local/lib, too.
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH

# Configure, with static binaries and target directroy
./configure --enable-static --enable-all-static --prefix="$PWD/build" || exit 1

# build the programs
make || exit 1

# Install them to the target directory (doesn't require root)
make install || exit 1

# Create a package of the static binaries
VERSION=$(grep '#define VERSION' config.h | sed 's/.*"\(.*\)".*/\1/')
TARBALL=fastx_toolkit_${VERSION}_binaries_$(uname -s)_$(uname -r)_$(uname -m).tar.bz2
cd build || exit 1
tar -cjvf "../$TARBALL" ./bin/* || exit 1
cd ..

echo "Static Binaries are installed in:"
echo "  $PWD/build/bin"
echo
echo "Static binaries tarball:"
echo "  $TARBALL"


