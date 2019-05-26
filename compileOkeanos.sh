#!/usr/bin/env bash

cp CMakeListsOkeanos.txt CMakeLists.txt
cd cmake-build-debug
cmake ..
make

rm -rf tests
cp -r ../tests .

cp -f ../batch.batch .