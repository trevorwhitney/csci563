#!/bin/sh

mpiexec -np 8 --hostfile hostfile cg.A.2 > results/cg.A.2.txt
mpiexec -np 8 --hostfile hostfile cg.S.2 > results/cg.S.2.txt
mpiexec -np 8 --hostfile hostfile ep.B.2 > results/ep.B.2.txt
mpiexec -np 8 --hostfile hostfile ep.W.2 > results/ep.W.2.txt
mpiexec -np 8 --hostfile hostfile is.C.2 > results/is.C.2.txt
mpiexec -np 8 --hostfile hostfile lu.A.2 > results/lu.A.2.txt
mpiexec -np 8 --hostfile hostfile lu.S.2 > results/lu.S.2.txt
mpiexec -np 8 --hostfile hostfile mg.B.2 > results/mg.B.2.txt
mpiexec -np 8 --hostfile hostfile mg.W.2 > results/mg.W.2.txt
mpiexec -np 8 --hostfile hostfile cg.B.2 > results/cg.B.2.txt
mpiexec -np 8 --hostfile hostfile cg.W.2 > results/cg.W.2.txt
mpiexec -np 8 --hostfile hostfile ep.C.2 > results/ep.C.2.txt
mpiexec -np 8 --hostfile hostfile is.A.2 > results/is.A.2.txt
mpiexec -np 8 --hostfile hostfile is.S.2 > results/is.S.2.txt
mpiexec -np 8 --hostfile hostfile lu.B.2 > results/lu.B.2.txt
mpiexec -np 8 --hostfile hostfile lu.W.2 > results/lu.W.2.txt
mpiexec -np 8 --hostfile hostfile mg.C.2 > results/mg.C.2.txt
mpiexec -np 8 --hostfile hostfile cg.C.2 > results/cg.C.2.txt
mpiexec -np 8 --hostfile hostfile ep.A.2 > results/ep.A.2.txt
mpiexec -np 8 --hostfile hostfile ep.S.2 > results/ep.S.2.txt
mpiexec -np 8 --hostfile hostfile is.B.2 > results/is.B.2.txt
mpiexec -np 8 --hostfile hostfile is.W.2 > results/is.W.2.txt
mpiexec -np 8 --hostfile hostfile lu.C.2 > results/lu.C.2.txt
mpiexec -np 8 --hostfile hostfile mg.A.2 > results/mg.A.2.txt
mpiexec -np 8 --hostfile hostfile mg.S.2 > results/mg.S.2.txt
