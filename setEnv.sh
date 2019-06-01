#!/usr/bin/env bash

cat < ../welcomeMessage
module swap PrgEnv-cray PrgEnv-gnu
module swap gcc/4.9.3 gcc/7.3.0

#module swap PrgEnv-gnu PrgEnv-cray
#module swap gcc/7.3.0 gcc/4.9.3