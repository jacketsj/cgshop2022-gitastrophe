# CG:SHOP 2022 - gitastrophe

This repository contains the source code of team gitastrophe in [CG:SHOP 2022](https://cgshop.ibr.cs.tu-bs.de/competition/cg-shop-2022/).

For more information, see our paper in the SoCG 2022 proceedings.

# Usage instructions

## Compiling/Building

We rely on `cmake` to build our code. If you don't have `cmake`, you can obtain it on an Ubuntu-based system with:
```sh
sudp apt-get install cmake
```

There should be no additional hidden dependencies (except a modern version of `g++`) for the primary executables,
and you should be able to compile with:
```
./build.sh
```
from the main directory.

### Compiling GPU Code

While it is not discussed in our paper, we did also explore running our primary conflict-optimizer algorithm on a GPU.
No performance improvements over the CPU-only code were obtained,
but only a handful of different possible methods for parallelization were explored.
This part of the code depends on [hipSYCL](https://github.com/illuhad/hipSYCL).
As of writing, the best way to install hipSYCL is through their repositories,
instructions for this method and others available [here](https://github.com/illuhad/hipSYCL/blob/develop/install/scripts/README.md).
Note that in order to run on GPU, the `/opt/hipSYCL/etc/hipSYCL/syclcc.json` file must be edited to specify a `default-platform`, and possibly some other parameters depending on the platform used and system configuration.

## Obtaining Instances

In order to run our application,
the `instances` and `dimacs-instances` folders should first be populated.

The `instances` folder should be populated with the [data available for download from the CG:SHOP 2022 web page](https://cgshop.ibr.cs.tu-bs.de/static/competitions/cgshop_2022/instances.zip).

The `dimacs-instances` folder should be populated with the [archived data for the second DIMACS challenge colouring instances](http://archive.dimacs.rutgers.edu/pub/challenge/graph/benchmarks/color/).

## Running

To run on the contest geometric instances, do one of the following:
```sh
mkdir temp                                                 # Create a folder to store the solutions (can be changed in encodings.h)
ls instances/*.json | build/simple                         # Initialize every instance with a simple greedy colouring
ls instances/*.json | build/simple -c                      # Run the conflict optimizer on each instance
ls instances/*.json | build/simple -t 5 -c                 # Run the conflict optimizer on each instance, with 5 threads
ls instances/sqrpecn3218.instance.json | build/simple -o   # Run a simple local search algorithm on the instance sqrpecn3218
```
There are a large number of other command-line options available.
A full list of them can be found in the source code `src/simple.cpp`.
Regardless of the algorithm run, improvements over the current best solution (if any) will be stored in the `temp` directory (make sure to create this directory before running!).

To run on the dimacs vertex-colouring instances, do
```sh
ls dimacs-instances/*.col | build/dimacs -t <NUM_THREADS>
```
