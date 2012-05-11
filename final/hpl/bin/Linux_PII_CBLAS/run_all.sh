#!/bin/sh
echo "Starting 128"
cp HPL_128.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
echo "Starting 144"
cp HPL_144.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
echo "Starting 160"
cp HPL_160.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
echo "Starting 176"
cp HPL_176.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
echo "Starting 192"
cp HPL_192.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
echo "Starting 208"
cp HPL_208.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
echo "Starting 224"
cp HPL_224.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
echo "Starting 240"
cp HPL_240.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
echo "Starting 256"
cp HPL_256.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
