#include "simulation.h"
#include "vector.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>


int read_config(const char *filename, struct config *conf) {
  FILE *ff = fopen(filename, "r");
  if (!ff) {
    perror("Error opening file");
    return 0;
  }


  char line[256];
  while (fgets(line, sizeof(line), ff)) {
  // Ignore comments and empty lines
    if (line[0] == '#' || line[0] == '\n') {
    continue;
    }

    // Parse key and value
    char key[50];
    char value[50];
    if (sscanf(line, "%49s = %49s", key, value) == 2) {
      if (strcmp(key, "print") == 0)
        conf->print = atoi(value);
      else if (strcmp(key, "number") == 0)
        conf->number = atoi(value);
      else if (strcmp(key, "saveflag") == 0)
        conf->saveflag = atoi(value);
      else if (strcmp(key, "xyzflag") == 0)
        conf->xyzflag = atoi(value);
      else if (strcmp(key, "start") == 0)
        conf->start = atoi(value);
      else if (strcmp(key, "step") == 0)
        conf->step = atoi(value);
      else if (strcmp(key, "stop") == 0)
        conf->stop = atoi(value);
      else if (strcmp(key, "helium_number") == 0)
        conf->helium_number = atoi(value);
      else if (strcmp(key, "mass") == 0)
        conf->mass = atof(value);
      else if (strcmp(key, "velocity") == 0)
        conf->velocity = atof(value);
      else if (strcmp(key, "dt") == 0)
        conf->dt = atof(value);
      else if (strcmp(key, "max_time") == 0)
        conf->max_time = atof(value);
      else if (strcmp(key, "icd_dist") == 0)
        conf->icd_dist = atof(value);
      else if (strcmp(key, "force_file") == 0)
        strcpy(conf->force_file, value);
    }
  }
  return 1;
  fclose(ff);
}

double distance_between_pair(const struct Particle_Pair *particle) {
    // Position of particle 1
    double x1   = particle->p1_pos.x;
    double y1   = particle->p1_pos.y;
    double z1   = particle->p1_pos.z;

    // Position of particle 2
    double x2   = particle->p2_pos.x;
    double y2   = particle->p2_pos.y;
    double z2   = particle->p2_pos.z;

    double xx = (x1 - x2);
    double yy = (y1 - y2);
    double zz = (z1 - z2);

    return sqrt(xx*xx + yy*yy + zz*zz);
}

int initialize_particle_pair(const double radius,
                             const struct config *conf,
                             struct Particle_Pair *particle) {
  double dist = 0.0;
  struct Particle_Pair *par = NULL;
  Quat q1 = {0, 0, 0, 0};
  Vector3D pos1 = {0, 0, 0};
  Quat q2 = {0, 0, 0, 0};
  Vector3D pos2 = {0, 0, 0};
  Vector3D a = {0, 0, 0};
  Vector3D res = {0, 0, 0};
  double ang_vel = (1E-5 * conf->velocity) / radius; // Unit fs¯¹

  // Setup random number generator
  srand(100401430);
  // srand(time(NULL));

  for (int i = 0; i < conf->number; ++i) {
      par = particle + i;
      // Generate random Quan
      quat_random(&q1);
      quat_to_pos(&q1, &pos1);
      scalar_multiply(radius, &pos1, &pos1);

      do {
        quat_random(&q2);
        quat_to_pos(&q2, &pos2);
        scalar_multiply(radius, &pos2, &pos2);
        dist = distance(&pos1, &pos2);
      } while (dist < conf->icd_dist);


      // Set pair position
      par->index = i;
      par->success = 0;
      par->time = 0;
      par->distance = dist;

      par->p1_qua = q1;
      par->p2_qua = q2;
      par->p1_pos = pos1;
      par->p2_pos = pos2;

      // Generate random angular velocity
      par->p1_velocity = 1E-5 * conf->velocity; // Convert m/s to Å/fs
      par->p2_velocity = 1E-5 * conf->velocity; // Convert m/s to Å/fs

      // Need to find the perpenducular vector to position
      // Create a vector 'a' which is not parallel to position vector
      if (par->p1_pos.x != 0 || par->p1_pos.y != 0) {
        a.x = -1 * par->p1_pos.y;
        a.y = par->p1_pos.x;
        a.z = 0;
      } else {
        a.x = 0;
        a.y = -1 * par->p1_pos.z;
        a.z = par->p1_pos.y;
      }
      // Take a cross product
      cross_product(&a, &(par->p1_pos), &res);
      double mag = norm(&res);
      // Create angular Velocity vector
      scalar_multiply(ang_vel/mag, &res, &(par->p1_angvel));

      // Need to find the perpenducular vector to position
      // Create a vector 'a' which is not parallel to position vector
      if (par->p2_pos.x != 0 || par->p2_pos.y != 0) {
        a.x = -1 * par->p2_pos.y;
        a.y = par->p2_pos.x;
        a.z = 0;
      } else {
        a.x = 0;
        a.y = -1 * par->p2_pos.z;
        a.z = par->p2_pos.y;
      }
      // Take a cross product
      cross_product(&a, &(par->p2_pos), &res);
      mag = norm(&res);
      // Create angular Velocity vector
      scalar_multiply(ang_vel/mag, &res, &(par->p2_angvel));

      // Calculate velocity
      cross_product(&(par->p1_pos), &(par->p1_angvel), &(par->p1_vel));
      cross_product(&(par->p2_pos), &(par->p2_angvel), &(par->p2_vel));

      // Set force to zero
      par->p1_force.x = 0;
      par->p1_force.y = 0;
      par->p1_force.z = 0;
      par->p2_force.x = 0;
      par->p2_force.y = 0;
      par->p2_force.z = 0;
  }

  return 0;
}

double find_force(const double distance, const int list_size,
                  const double *rlist, const double *force_list) {
  int indx = 0;
  if (distance < rlist[0]) indx = 0;
  else if (distance > rlist[list_size - 1]) indx = list_size - 2;
  else indx = (int)((distance - rlist[0]) / 0.001);
  return -1 * force_list[indx] * 9.6485332E-03;
}

void find_accelration(const double mass,
                      const double radius,
                      const Vector3D *pos,
                      const Vector3D *force,
                      Vector3D *acc) {

  Vector3D torque = {0, 0, 0};
  cross_product(pos, force, &torque);
  scalar_multiply(1.0 / mass / radius / radius, &torque, acc);
}

void update_angvel(const double timestep,
                   const Vector3D *oldangvel,
                   const Vector3D *acc,
                   Vector3D *newangvel) {
  Vector3D temp;
  scalar_multiply(timestep, acc, &temp);
  add_vectors(oldangvel, &temp, newangvel);
}

void update_orientation(const double radius,
                        const double timestep,
                        Quat *orient,
                        const Vector3D *angvel,
                        Vector3D *pos) {
  Quat Omega = {0, angvel->x, angvel->y, angvel->z};
  Quat temp = {0, 0, 0, 0};
  quat_product(&Omega, orient, &temp);
  quat_scalar_multiply(&temp, timestep/2.0, &temp);
  quat_add(orient, &temp, orient);
  quat_scalar_multiply(orient,
                       1.0 / quat_norm(orient),
                       orient);
  quat_to_pos(orient, pos);
  scalar_multiply(radius, pos, pos);
}

void save_particle_pair(FILE *fpointer, const struct Particle_Pair *particle) {
      fprintf(fpointer, "%d\t%d\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\n",
              particle->index, particle->success,
              particle->time, particle->distance,
              particle->p1_velocity, particle->p2_velocity,
              particle->p1_qua.a,
              particle->p1_qua.x,
              particle->p1_qua.y,
              particle->p1_qua.z,
              particle->p2_qua.a,
              particle->p2_qua.x,
              particle->p2_qua.y,
              particle->p2_qua.z,
              particle->p1_pos.x,
              particle->p1_pos.y,
              particle->p1_pos.z,
              particle->p2_pos.x,
              particle->p2_pos.y,
              particle->p2_pos.z,
              particle->p1_vel.x,
              particle->p1_vel.y,
              particle->p1_vel.z,
              particle->p2_vel.x,
              particle->p2_vel.y,
              particle->p2_vel.z,
              particle->p1_force.x,
              particle->p1_force.y,
              particle->p1_force.z,
              particle->p2_force.x,
              particle->p2_force.y,
              particle->p2_force.z,
              particle->p1_angvel.x,
              particle->p1_angvel.y,
              particle->p1_angvel.z,
              particle->p2_angvel.x,
              particle->p2_angvel.y,
              particle->p2_angvel.z);

}

void create_xyz_file(const int tindex, const struct Particle_Pair *particle) {
  int index = particle->index;

  char filename[80];
  sprintf(filename, "./results/trajectories_%06d/traj_%09d.xyz", index, tindex);
  FILE *ffile = fopen(filename, "w");
  fprintf(ffile, "2\n");
  fprintf(ffile, "Trajectory of particles\n");
  fprintf(ffile, "C %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
          particle->p1_pos.x,
          particle->p1_pos.y,
          particle->p1_pos.z,
          particle->p1_vel.x,
          particle->p1_vel.y,
          particle->p1_vel.z,
          particle->p1_force.x,
          particle->p1_force.y,
          particle->p1_force.z);

  fprintf(ffile, "C %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n",
          particle->p2_pos.x,
          particle->p2_pos.y,
          particle->p2_pos.z,
          particle->p2_vel.x,
          particle->p2_vel.y,
          particle->p2_vel.z,
          particle->p2_force.x,
          particle->p2_force.y,
          particle->p2_force.z);
  fclose(ffile);
}

void simulate_particle(const struct config *conf,
                       const double radius, const int list_size,
                       const double *rlist, const double *force_list,
                       struct Particle_Pair *particle) {

  double dt = conf->dt;
  double forcemag = 0.0;
  Vector3D ang_vel1 = {0., 0., 0.};
  Vector3D acc1 = {0., 0., 0.};
  Vector3D ang_vel2 = {0., 0., 0.};
  Vector3D acc2 = {0., 0., 0.};

  char ffname[80];
  FILE *ffile = NULL;
  int tindex = 0;


  if (conf->xyzflag) {
    char dirname[80];
    sprintf(dirname, "./results/trajectories_%06d", particle->index);

    // Create result directory to keep results
    struct stat st = {0};
    if (stat(dirname, &st) == -1) {
      mkdir(dirname, 0700);
    }
  }

  if (conf->saveflag) {
    sprintf(ffname, "./results/trajectory_%06d.txt", particle->index);
    ffile = fopen(ffname, "w");
    fprintf(ffile,
            "Index\tSuccess\tTime\tDistance\t"
            "Vel_1\tVel_2\t"
            "QA_1\tQX_1\tQY_1\tQZ_1\t"
            "QA_2\tQX_2\tQY_2\tQZ_2\t"
            "X_1\tY_1\tZ_1\t"
            "X_2\tY_2\tZ_2\t"
            "VX_1\tVY_1\tVZ_1\t"
            "VX_2\tVY_2\tVZ_2\t"
            "FX_1\tFY_1\tFZ_1\t"
            "FX_2\tFY_2\tFZ_2\t"
            "WX_1\tWY_1\tWZ_1\t"
            "WX_2\tWY_2\tWZ_2\n");
  }

  particle->time = 0;
  while (particle->time < conf->max_time && !particle->success) {
    // Save trajectory in text file
    if (conf->saveflag)
      save_particle_pair(ffile, particle);
    if (conf->xyzflag)
      create_xyz_file(tindex, particle);

    // Find Force
    forcemag = find_force(particle->distance, list_size, rlist, force_list);
    // If force is less than 1E-10 amu.Å.fs¯² (10 neV/Å, 1E-17 N)
    // Consider it zero
    if (fabs(forcemag) - 1E-10 > 0) {
      dt = 1.0;
      // Step 1: Find the torque on particles
      particle->p1_force.x = forcemag * (particle->p1_pos.x - particle->p2_pos.x) / particle->distance;
      particle->p1_force.y = forcemag * (particle->p1_pos.y - particle->p2_pos.y) / particle->distance;
      particle->p1_force.z = forcemag * (particle->p1_pos.z - particle->p2_pos.z) / particle->distance;
      scalar_multiply(-1, &(particle->p1_force), &(particle->p2_force));

      find_accelration(conf->mass, radius, &(particle->p1_pos), &(particle->p1_force), &acc1);
      find_accelration(conf->mass, radius, &(particle->p2_pos), &(particle->p2_force), &acc2);

      // Step 2: Calculate angular velocity at timestamp Δt/2
      update_angvel(dt/2.0, &(particle->p1_angvel), &acc1, &ang_vel1);
      update_angvel(dt/2.0, &(particle->p2_angvel), &acc2, &ang_vel2);
    } else {
      ang_vel1.x = particle->p1_angvel.x;
      ang_vel1.y = particle->p1_angvel.y;
      ang_vel1.z = particle->p1_angvel.z;
      ang_vel2.x = particle->p2_angvel.x;
      ang_vel2.y = particle->p2_angvel.y;
      ang_vel2.z = particle->p2_angvel.z;
    }

    // Step 3: Update Orientation
    update_orientation(radius, dt, &(particle->p1_qua),
                       &ang_vel1, &(particle->p1_pos));
    update_orientation(radius, dt, &(particle->p2_qua),
                       &ang_vel2, &(particle->p2_pos));

    // Check if they met
    particle->distance = distance_between_pair(particle);
    if (particle->distance < conf->icd_dist) {
      particle->success = 1;
      break;
    }
 
    // Check force again
    forcemag = find_force(particle->distance, list_size, rlist, force_list);

    if (fabs(forcemag) - 1E-10 > 0) {
      // Step 4: Find the torque at new position
      particle->p1_force.x = forcemag * (particle->p1_pos.x - particle->p2_pos.x) / particle->distance;
      particle->p1_force.y = forcemag * (particle->p1_pos.y - particle->p2_pos.y) / particle->distance;
      particle->p1_force.z = forcemag * (particle->p1_pos.z - particle->p2_pos.z) / particle->distance;
      scalar_multiply(-1, &(particle->p1_force), &(particle->p2_force));

      find_accelration(conf->mass, radius, &(particle->p1_pos), &(particle->p1_force), &acc1);
      find_accelration(conf->mass, radius, &(particle->p2_pos), &(particle->p2_force), &acc2);

      // Step 5: Calculate angular velocity at Δt
      update_angvel(dt/2.0, &ang_vel1, &acc1, &(particle->p1_angvel));
      update_angvel(dt/2.0, &ang_vel2, &acc2, &(particle->p2_angvel));
    } else {
      particle->p1_angvel.x = ang_vel1.x;
      particle->p1_angvel.y = ang_vel1.y;
      particle->p1_angvel.z = ang_vel1.z;
      particle->p2_angvel.x = ang_vel2.x;
      particle->p2_angvel.y = ang_vel2.y;
      particle->p2_angvel.z = ang_vel2.z;
    }

    // Calculate new velocity magnitute
    cross_product(&(particle->p1_pos), &(particle->p1_angvel), &(particle->p1_vel));
    particle->p1_velocity = norm(&(particle->p1_vel));
    cross_product(&(particle->p2_pos), &(particle->p2_angvel), &(particle->p2_vel));
    particle->p2_velocity = norm(&(particle->p2_vel));

    particle->time += dt;
    tindex += 1;
  }

  if (conf->saveflag) {
    save_particle_pair(ffile, particle);
    fclose(ffile);
  }

  if (conf->xyzflag)
    create_xyz_file(tindex, particle);
}
