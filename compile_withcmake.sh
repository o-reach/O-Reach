#!/bin/bash

INC=0
if [ "$1" == "-i" ]
then
	INC=1
fi

NCORES=4
unamestr=`uname`
if [[ "$unamestr" == "Linux" ]]; then
        NCORES=`grep -c ^processor /proc/cpuinfo`
fi

if [[ "$unamestr" == "Darwin" ]]; then
        NCORES=`sysctl -n hw.ncpu`
fi

if [ ${INC} -eq 0 ]
then
	rm -rf deploy
	rm -rf build
	mkdir build
else
	test -d build || mkdir build
fi
cd build
cmake ../
make -j $NCORES VERBOSE=1
cd ..

if [ ${INC} -eq 0 ]
then
	mkdir deploy
else
  test -d deploy || mkdir deploy
fi
cp ./build/reachability deploy/

