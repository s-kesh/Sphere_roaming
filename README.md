# Roaming Simulation on the Surface of a Sphere

This repository provides a Monte Carlo simulation framework for modeling the movement of particles constrained to the surface of a Helium droplet. It assume Helium droplet to be just a solid sphere.

## Simulation Steps

1. Calculate the radius \( r \) for the given number of He atoms.  
2. Determine the average number of excitations in the droplet based on the given laser intensity, pulse width, and cross-section.  
3. Generate the number of excitations using a Poissonian distribution.  
4. Distribute their positions uniformly on the surface of the sphere.  
5. Assign initial velocities based on the Maxwell-Boltzmann distribution.  
6. Let the excitations roll on the surface under the influence of force.  
7. Stop the simulation when all pairs of particles have met.  

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

- **HeRoam**: The main simulation binary. It takes a configuration file as input to execute the Monte Carlo simulation.

```sh
./HeRoam <config.conf>
```

- **h5read**: A utility to read and visualize HDF5 trajectory data. This tool uses gnuplot to plot particle positions at the specified time interval.

```sh
./h5read <filename> <time_interval>
```

## Configuration File

Most of the configuration file is well-documented with inline comments.  
However, ensure that the force file is placed at the correct path before running the simulation.  

**Note:** The force file is not included in this repository, as its size exceeds 100MB.