#!/bin/bash
# for use in bigdata machine only
rm -r build
mkdir build
cd build
cmake ../source -DCMAKE_INSTALL_PREFIX=${LOCAL} -DCMAKE_INSTALL_LIBDIR=${LOCAL}/lib
make install
