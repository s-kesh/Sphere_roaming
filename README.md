# Roaming Simulation on the Surface of a Sphere

This repository provides a Monte Carlo simulation framework for modeling the movement of particles constrained to the surface of a sphere.

## Dependencies

To build and run the simulation, the following dependencies are required:

- **C Compiler**  
- **Meson Build System** (can be converted to a simple Makefile if needed)  
- **HDF5 Library** (for storing individual trajectories)  
- **OpenMP** (optional, for parallel execution of multiple trajectories)  
- **gnuplot** (optional, for visualizing trajectories)  

## Building the Project

To build the simulation, run:

```sh
meson setup build
meson compile -C build
```

## Usage

The compilation generates two executables: ```HeRoam``` and ```h5read```.

- **HeRoam**: The main simulation binary. It takes a configuration file as input to execute the Monte Carlo simulation. The configuration file is well-documented with inline comments.
- **h5read**: A utility to read and visualize HDF5 trajectory data. This tool uses gnuplot to plot particle positions at the specified time interval.