#!/usr/bin/env bash

../compileOkeanos.sh
rm output/*
rm error/*
rm logger/*
sbatch batch.batch --account=gc72-18 --reservation=rzadca