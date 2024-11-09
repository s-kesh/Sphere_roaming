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
  double mass;
  double velocity;
  double dt;
  double max_time;
  double icd_dist;
  char force_file[50];
};

struct Particle_Pair {
    unsigned int success;
    unsigned int index;
    double time;
    double distance;

    double p1_velocity;
    double p2_velocity;

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
double distance_between_pair(const struct Particle_Pair *particle);
int initialize_particle_pair(const double radius,
                             const struct config *conf,
                             struct Particle_Pair *particle);
void simulate_particle(const struct config *conf,
                       const double radius, const int list_size,
                       const double *rlist, const double *force_list,
                       struct Particle_Pair *particle);
double find_force(const double distance, const int list_size,
                  const double *rlist, const double *force_list);

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
