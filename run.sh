#!/bin/bash
# hahahaah
make clean
make
cd ..
cd ./spin1_1D_OL_NC_ct
make clean
make
./myappname -0.01 > pso.file 2>&1 &
cd ..
cd ./spin1_1D_OL_NC
./myappname -0.01 > pso.file 2>&1 &
