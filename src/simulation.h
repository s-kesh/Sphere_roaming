#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "vector.h"
#include <stdio.h>

/*
Structure to hold the configuration read from the configuration file.
*/
struct config {
  int print;
  int no; // no of particles for each droplet
  int number;
  int saveflag;
  int start;
  int step;
  int stop;
  int helium_number;
  long seed;
  double mass;
  double velocity;
  double dt;
  double max_time;
  double icd_dist;
  double force_grid;
  double force_grid_start;
  long force_grid_length;
  char force_file[50];
};

/*
Structure to hold the particle properties.
*/
typedef struct {
  unsigned int index;
  unsigned int success;
  double time;
  double force_mag;
  double velocity_sq;
  Quat orientation;
  Vector3D position;
  Vector3D velocity;
  Vector3D ang_velocity;
  Vector3D force;
} Particle;

/*
Structure to hold the particles.
*/
typedef struct {
  unsigned int no; // Number of particles to hold
  unsigned int index; // Keep track of MC trajectory
  Particle *particle; // Particle pointer to first particle.
} Particles;

/*
Function to read the configuration file and store the values in the configuration structure.
@param filename: Name of the configuration file.
@param conf: Configuration structure to store the values.
@return: 1 if successful, 0 otherwise.
*/
int read_config(const char *filename, struct config *conf);

/*
Function to generate random velocities according to the Maxwell-Boltzmann distribution.
@param conf: Configuration structure containing the velocity.
@return: Random velocity.
*/
double generate_velocity(const struct config *conf);

/*
Function to calculate the force between the particles.
@param particles: Pointer to the particles structure.
@param conf: Configuration structure containing the force values.
@param Force_list: Array containing the force values for different distances.
*/
void calculate_force(Particles *particles,
                     const struct config *conf,
                     const double *Force_list);

/*
Function to initialize the particles.
@param radius: Radius of the droplet.
@param conf: Configuration structure containing the configuration read from the file.
@param Force_list: Array containing the force values for different distances.
@param particles: Pointer to the particles structure.
@return: 0 if successful, 1 otherwise.
*/
int initialize_particles(const double radius,
                             const struct config *conf,
                             const double *Force_list,
                             Particles *particles);

/*
Function to free the memory allocated for the particles.
@param particles: Pointer to the particles structure.
*/                             
void free_particles(Particles *particles);

/*
Function to simulate the particles using the velocity verlet algorithm.
@param conf: Configuration structure containing the configuration read from the file.
@param radius: Radius of the droplet.
@param force_list: Array containing the force values for different distances.
@param particles: Pointer to the particles structure.
*/
void simulate_particles(const struct config *conf,
                       const double radius,
                       const double *force_list,
                       Particles *particles);

/*
Function to find the acceleration of the particle.
@param mass: Mass of the particle.
@param radius: Radius of the droplet.
@param pos: Pointer to the position vector of the particle.
@param force: Pointer to the force vector on the particle.
@param acc: Pointer to the acceleration vector of the particle.
*/
void find_accelration(const double mass, const double radius,
                      const Vector3D *pos, const Vector3D *force,
                      Vector3D *acc);

/*
Function to update the angular velocity of the particle.
@param timestep: Timestep of the simulation.
@param oldangvel: Pointer to the old angular velocity vector.
@param acc: Pointer to the acceleration vector of the particle.
@param newangvel: Pointer to the new angular velocity vector.
*/
void update_angvel(const double timestep, const Vector3D *oldangvel,
                   const Vector3D *acc, Vector3D *newangvel);

/*
Function to update the orientation of the particle.
@param radius: Radius of the droplet.
@param timestep: Timestep of the simulation.
@param orient: Pointer to the orientation quaternion of the particle.
@param angvel: Pointer to the angular velocity vector of the particle.
@param pos: Pointer to the position vector of the particle.
*/
void update_orientation(const double radius, const double timestep,
                        Quat *orient, const Vector3D *angvel,
                        Vector3D *pos);

/*
Function to save the particle properties in a file.
@param fpointer: Pointer to the file to save the particle properties.
@param particles: Pointer to the particles structure.
*/
void save_particles(FILE *fpointer, const Particles *particles);
#endif
