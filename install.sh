#!/bin/bash
# for use in bigdata machine only
cmake ../source -DCMAKE_INSTALL_PREFIX=${LOCAL} -DCMAKE_INSTALL_LIBDIR=${LOCAL}/lib
make install
