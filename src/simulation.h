#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "vector.h"
#include <stdio.h>

/*
Structure to hold the configuration read from the configuration file.
*/
struct config {
  int print;
  float hv;
  float intensity;
  float cross_section;
  float pulse_width;
  int number;
  int saveflag;
  int start;
  int step;
  int stop;
  int helium_number;
  long seed;
  float mass;
  float velocity;
  float dt;
  float max_time;
  float icd_dist;
  float force_grid;
  float force_grid_start;
  long force_grid_length;
  char force_file[50];
};

/*
Structure to hold the particle properties.
*/
typedef struct {
  unsigned int index;
  unsigned int success;
  float time;
  float force_mag;
  float velocity_sq;
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
Function to initialize the particles.
@param he_number: Number of helium atoms
@param conf: Configuration structure containing the configuration read from the file.
@param Force_list: Array containing the force values for different distances.
@param particles: Pointer to the particles structure.
@return: 0 if successful, 1 otherwise.
*/
int initialize_particles(const unsigned long long he_number,
                             const struct config *conf,
                             const float *Force_list,
                             Particles *particles);

/*
Function to free the memory allocated for the particles.
@param particles: Pointer to the particles structure.
*/
void free_particles(Particles *particles);

/*
Function to simulate the particles using the velocity verlet algorithm.
@param conf: Configuration structure containing the configuration read from the file.
@param he_number: Number of helium atoms
@param force_list: Array containing the force values for different distances.
@param particles: Pointer to the particles structure.
*/
void simulate_particles(const struct config *conf,
                        const unsigned long long he_number,
                       const float *force_list,
                       Particles *particles);

/*
Function to save the particle properties in a file.
@param fpointer: Pointer to the file to save the particle properties.
@param particles: Pointer to the particles structure.
*/
void save_particles(FILE *fpointer, const Particles *particles);
#endif
