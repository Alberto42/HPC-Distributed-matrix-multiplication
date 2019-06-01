#!/usr/bin/env bash

../compileOkeanos.sh
rm output/*
rm error/*
rm logger/*
rm -rf matrixmul+*
sbatch batch.batch --account=gc72-18 --reservation=rzadca