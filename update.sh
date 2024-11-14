#!/bin/bash
# for use in bigdata machine only

cd $(dirname "$0")
git add .
git commit -m "${2}"
git push origin $1
cd build
cmake ../source -DCMAKE_INSTALL_PREFIX=${LOCAL} -DCMAKE_INSTALL_LIBDIR=${LOCAL}/lib
make install