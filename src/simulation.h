#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "vector.h"
#include <stdio.h>

struct config {
  int print;
  int no; // no of particles for each droplet
  int number;
  int saveflag;
  int xyzflag;
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

typedef struct {
  unsigned int no;
  unsigned int index;
  unsigned int *success;
  double *time;

  double *velocity_sq;
  double *force_mag;

  Quat *orientation;
  Vector3D *position;
  Vector3D *velocity;
  Vector3D *ang_velocity;
  Vector3D *force;
} Particle ;

int read_config(const char *filename, struct config *conf);
double generate_velocity(const struct config *conf);
void calculate_force(Particle *particle,
                     const struct config *conf,
                     const double *Force_list);
int initialize_particle_pair(const double radius,
                             const struct config *conf,
                             const double *Force_list,
                             Particle *particle);
void free_particles(Particle *particle);
void simulate_particle(const struct config *conf,
                       const double radius,
                       const double *force_list,
                       Particle *particle);
void find_accelration(const double mass, const double radius,
                      const Vector3D *pos, const Vector3D *force,
                      Vector3D *acc);
void update_angvel(const double timestep, const Vector3D *oldangvel,
                   const Vector3D *acc, Vector3D *newangvel);
void update_orientation(const double radius, const double timestep,
                        Quat *orient, const Vector3D *angvel,
                        Vector3D *pos);

void save_particle_pair(FILE *fpointer, const Particle *particle);
void create_xyz_file(const int tindex, const Particle *particle);
#endif
