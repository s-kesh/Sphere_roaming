#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "vector.h"
#include <stdio.h>

struct config {
  int print;
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
  double velratio;
  double dt;
  double max_time;
  double icd_dist;
  double force_grid;
  double force_grid_start;
  long force_grid_length;
  char force_file[50];
};

struct Particle_Pair {
  unsigned int success;
  unsigned int index;
  double time;
  double distance_sq;

  double p1_velocity_sq;
  double p2_velocity_sq;
  double force;

  Vector3D p1_pos;
  Vector3D p2_pos;
  Quat p1_qua;
  Quat p2_qua;
  Vector3D p1_angvel;
  Vector3D p2_angvel;

  Vector3D p1_vel;
  Vector3D p2_vel;

  Vector3D p1_force;
  Vector3D p2_force;
};

int read_config(const char *filename, struct config *conf);
double distancesq_between_pair(const struct Particle_Pair *particle);
double generate_velocity(const struct config *conf);
int initialize_particle_pair(const double radius,
                             const struct config *conf,
                             const double *Force_list,
                             struct Particle_Pair *particle);
void simulate_particle(const struct config *conf,
                       const double radius,
                       const double *force_list,
                       struct Particle_Pair *particle);
void find_accelration(const double mass, const double radius,
                      const Vector3D *pos, const Vector3D *force,
                      Vector3D *acc);
void update_angvel(const double timestep, const Vector3D *oldangvel,
                   const Vector3D *acc, Vector3D *newangvel);
void update_orientation(const double radius, const double timestep,
                        Quat *orient, const Vector3D *angvel,
                        Vector3D *pos);

void save_particle_pair(FILE *fpointer, const struct Particle_Pair *particle);
void create_xyz_file(const int tindex, const struct Particle_Pair *particle);
#endif
