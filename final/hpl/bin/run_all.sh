#!/bin/sh
cp HPL_128.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
cp HPL_144.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
cp HPL_160.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
cp HPL_176.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
cp HPL_192.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
cp HPL_208.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
cp HPL_224.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
cp HPL_240.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl
cp HPL_256.dat HPL.dat
mpiexec -np 8 --hostfile hostfile xhpl