# Quick Start Guide for: OpenRBC - A Fast Simulator of Red Blood Cells at Protein-Resolution

## Compilation
x86 processor, recent Linux distribution, `g++` newer than 5.3.0:
```
cd <working_copy>/src
make
```
To override the compiler used
```
make CXX=path_to_compiler
```
For IBM Power 8 architecture
```
make ARCH=power8
```
For older Linux systems whose glibc do not provide an `std::aligned_alloc` function.
```
make AALLOC=1
```
A summary of arguments to `make`:

* `ARCH=[power8|x86]` compilation target architecture
* `CXX=[g++|xlC_r|icpc]` compiler
* `DEBUG=[0|1]` whether to turn on 'optimization for debug'
* `AALLOC=[0|1]` fall-back option for older linux kernels that do not provide the `std::aligned_alloc` function.

## Small example

```
cd <working_copy>
mkdir example-small
cd example-small
ln -s ../src/openrbc
OMP_NUM_THREADS=... ./openrbc -E 100 -t 10
```
To see a list of command line options:
```
./openrbc --help
```
OpenMP thread binding is vital to achieve optimal performance acorss multiple NUMA domains. This can be done by setting the environment variables:

| Compiler/Runtime     | Environment variable |
|----------------------|----------------------|
| Generic, OpenMP 4.5+ | OMP_PROC_BIND        |
| G++/GOMP             | GOMP_CPU_AFFINITY    |
| Intel/IOMP           | KMP_AFFINITY         |

In most cases placing consecutive threads on adjacent core/hardware threads would be optimal, however it is highly recommended that you test and tune thread binding before performing production runs. To figure out your socket/core topology refer to the `lstopo`, `hwloc-ls`, `hwloc-info` commnands from the `hwloc` package.

## Full-size example

```
cd <working_copy>/example-large
ln -s ../src/openrbc
OMP_NUM_THREADS=... ./openrbc -i trimesh -m rbc -E 100 -t ...
```
where
`-i` instruct the program to initialize using mesh files specified by the `-m` argument.

## Visualization

The initial structure file can be visualized with [VMD](http://www.ks.uiuc.edu/Research/vmd/) using the following TCL command.
```
topo readlammpsdata cell.data bond
```
The trajectory can be convert to a LAMMPS trajectory format by
```
cd <working_copy>/src
make orbc-util
cd <working_copy>/<case_directory>
../orbc-util convert cell.orbc cell.lammpstrj
```
which can then be appended to the initial structure in VMD by
```
mol addfile cell.lammpstrj autobonds 0
```
The LAMMPS trajectory file, however, is text-based and works poorly for whole-cell simulations containing millions of particles. In this case a binary format can be generated and appended in vmd:
```bash
../orbc-util convert cell.orbc cell.%06d.namdbin
```
```TCL
% adjust iteration bound and increment accordingly
for {set i 0} {$i < 1000000} {set i [expr $i + 100]} {
	mol addfile cell.[format "%06d" $i].namdbin type namdbin
}
```

[Video on YouTube](https://youtu.be/ahhvixWfRpM)
