# OpenRBC - A Fast Simulator of Red Blood Cells at Protein-Resolution

This repository contains the binary for performance evaluation of the OpenPOWER Challenge Hackathon project 'OpenRBC'.

## To run the OpenRBC demo on SuperVessel

```
cd <working_copy>/openrbc
OMP_NUM_THREADS=32 GOMP_CPU_AFFINITY=0-31 ./openrbc
```

The result data file can be visualized with [VMD](http://www.ks.uiuc.edu/Research/vmd/) using the following TCL command.

```
# Initial structure
topo readlammpsdata lipid.data bond
topo readlammpsdata protein.data bond
# Trajectory
mol new trajectory.lammpstrj autobonds 0
```

[Video on YouTube](https://youtu.be/ahhvixWfRpM)

## To run the legacy solver

```
cd <working_copy>/legacy
mpirun -np 8 ./md_rbc
```

The result files in Cfg/ can also be viewed with VMD.

## To compile from source

```
cd <working_copy>/src
make
```
Optional parameters for make
* `ARCH=[power8|x86]` compilation target architecture
* `CXX=[g++|xlC_r|icpc]` compiler
* `DEBUG=[0|1]` whether to turn on 'optimization for debug'
* `AALLOC=[0|1]` fall-back option for older linux kernels that do not provide the `std::aligned_alloc` API.
