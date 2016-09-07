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
mpirun -np 8 ./Art_md
```

The result files in Cfg/ can also be viewed with VMD.
