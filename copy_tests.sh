#!/usr/bin/env bash
pwd
rm -rf tests
cp -r ../tests .
cp -f batch.batch cmake-build-debug/batch.batch
