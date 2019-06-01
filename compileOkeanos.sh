#!/usr/bin/env bash

cd cmake-build-debug
cp ../CMakeListsOkeanos.txt ../CMakeLists.txt
cmake ..
make

rm -rf tests
cp -r ../tests .

cp -f ../batch.batch .
cp -f ../runOkeanos.sh .
cp -f ../createBuildDirs.sh .