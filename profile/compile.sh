#!/bin/bash
module purge
module load autoload
module load intelmpi/oneapi-2022--binary
source /cineca/prod/opt/compilers/intel/oneapi-2022/binary/setvars.sh

mpiifort jacobi-mpi-sendrecv.f90 -o jacobi-mpi-bare.exe

#export LD_PRELOAD=/g100/prod/opt/compilers/intel/oneapi-2022/binary/itac/2021.5.0/slib/libVT.so
mpiifort -trace -O2 jacobi-mpi-sendrecv.f90 -o jacobi-mpi-instr.exe
