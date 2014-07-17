#!/bin/sh
mkdir -p config
autoreconf --force --install -I config  -I m4
echo "Now run ./configure --prefix=$HOME && make install" && \
echo "Add --disable-shared to the configure line if building on Mac OS X"

