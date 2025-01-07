#include "simulation.h"
#include "vector.h"
#include "h5file.h"

#include <H5Ipublic.h>
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
      else if (strcmp(key, "no") == 0)
        conf->no = atoi(value);
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
      else if (strcmp(key, "seed") == 0)
        conf->seed = atoi(value);
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
      else if (strcmp(key, "force_grid") == 0)
        conf->force_grid = atof(value);
      else if (strcmp(key, "force_grid_start") == 0)
        conf->force_grid_start = atof(value);
      else if (strcmp(key, "force_grid_length") == 0)
        conf->force_grid_length = atoi(value);
      else if (strcmp(key, "force_file") == 0)
        strcpy(conf->force_file, value);
    }
  }
  fclose(ff);
  return 1;
}

void calculate_force(Particles *par,
                     const struct config *conf,
                     const double *Force_list) {

  const double inv_grid_size = 1.0 / conf->force_grid;
  const double force_start = conf->force_grid_start;
  const double icd_dist2 = conf->icd_dist*conf->icd_dist;
  const double max_dist_sq = conf->force_grid_start + (conf->force_grid_length - 2) * conf->force_grid;

  unsigned int num = par->no;
  double *dist_sq = (double *)malloc((num * (num-1) / 2) * sizeof(double));
  for (unsigned int j = 0; j < par->no; ++j) {
    for (unsigned int k = j + 1; k < par->no; ++k) {
      unsigned int index = j*num + k - (j + 1)*(j + 2)/2;
      double d = distance_sq(par->position + j, par->position + k);
      *(dist_sq + index) = d;
      if(d < icd_dist2) {
        par->success[j] = 1;
        par->success[k] = 1;
      }
      par->force[j].x = 0;
      par->force[j].y = 0;
      par->force[j].z = 0;

      par->force[k].x = 0;
      par->force[k].y = 0;
      par->force[k].z = 0;
    }
  }

  for (unsigned int j = 0; j < par->no; ++j) {
    if (par->success[j]) continue;
    Vector3D f = {0, 0, 0};
    for (unsigned int k = j + 1; k < par->no; ++k) {
      if (par->success[k]) continue;
      // Pair j,k
      unsigned int index = j*num + k - (j + 1)*(j + 2)/2;
      double d = *(dist_sq + index);
      if (d < max_dist_sq) {
        double scale = Force_list[(int)((d - force_start) * inv_grid_size)] / sqrt(d);
        f.x = scale * (par->position[j].x - par->position[k].x);
        f.y = scale * (par->position[j].y - par->position[k].y);
        f.z = scale * (par->position[j].z - par->position[k].z);

        par->force[j].x += f.x;
        par->force[j].y += f.y;
        par->force[j].z += f.z;


        par->force[k].x -= f.x;
        par->force[k].y -= f.y;
        par->force[k].z -= f.z;
      }
    }
    par->force_mag[j] = norm(par->force + j);
  }
  free(dist_sq);
}

int initialize_particle_pair(const int number,
                             const double radius,
                             const struct config *conf,
                             const double *Force_list,
                             Particles *particle) {
  Vector3D res = {0, 0, 0};

  const double ang_vel = (1E-5 * conf->velocity) / radius; // Unit fs¯¹
  const double velocity_sq = conf->velocity * conf->velocity * 1E-10; // Convert m/s to Å/fs

  // Setup random number generator
  if (conf->seed)
    srand(conf->seed);
  else
    srand(time(NULL));

  Particles *par = NULL;
  for (int i = 0; i < conf->number; ++i) {
      par = particle + i;

      // Initialise particles
      par->no = number;
      par->index = i;
      par->success = (unsigned int *)malloc(number * sizeof(unsigned int));
      par->time = (double *)malloc(number * sizeof(double));
      par->velocity_sq = (double *)malloc(number * sizeof(double));
      par->force_mag = (double *)malloc(number * sizeof(double));
      par->orientation = (Quat *)malloc(number * sizeof(Quat));
      par->position = (Vector3D *)malloc(number * sizeof(Vector3D));
      par->velocity = (Vector3D *)malloc(number * sizeof(Vector3D));
      par->ang_velocity = (Vector3D *)malloc(number * sizeof(Vector3D));
      par->force = (Vector3D *)malloc(number * sizeof(Vector3D));

      // Generate single particles
      for (int j = 0; j < number; ++j) {
        par->success[j] = 0;

        // Generate position
        quat_random(par->orientation + j);
        quat_to_pos(par->orientation + j, par->position + j);
        scalar_multiply(radius, par->position + j, par->position + j);

        // Generate velocity
        par->velocity_sq[j] = velocity_sq;
        create_perpend_vector(par->position + j, &res);
        scalar_multiply(ang_vel/norm(&res), &res, par->ang_velocity + j);
        cross_product(par->position + j, par->ang_velocity + j, par->velocity + j);

        par->force[j].x = 0;
        par->force[j].y = 0;
        par->force[j].z = 0;
      }

      // Calculate force
      calculate_force(par, conf, Force_list);
  }
  return 0;
}

void free_particles(Particles *particle) {
  free(particle->success);
  free(particle->time);
  free(particle->velocity_sq);
  free(particle->force_mag);
  free(particle->orientation);
  free(particle->position);
  free(particle->velocity);
  free(particle->ang_velocity);
  free(particle->force);
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

/*
// SPIRAL
void update_orientation(const double radius,
                        const double timestep,
                        Quat *orient,
                        const Vector3D *angvel,
                        Vector3D *pos) {
  Quat tempquat = {0, 0, 0, 0};
  Quat res = {0, 0, 0, 0};
  Vector3D temp = {0, 0, 0};
  double mag = norm(angvel);
  double theta = timestep/2.0 * mag;

  scalar_multiply(sin(theta)/mag, angvel, &temp);
  tempquat.a = cos(theta);
  tempquat.x = temp.x;
  tempquat.y = temp.y;
  tempquat.z = temp.z;

  quat_product(orient, &tempquat, &res);
  orient->a = res.a;
  orient->x = res.x;
  orient->y = res.y;
  orient->z = res.z;

  quat_to_pos(orient, pos);
  scalar_multiply(radius, pos, pos);
}
*/

// Velocity Verlet
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

void save_particle_pair(FILE *fpointer, const Particles *particle) {
  for (unsigned int i = 0; i < particle->no; ++i) {
    fprintf(fpointer,
            "%d\t%d\t%d\t%.10e\t"
            "%.10e\t%.10e\t%.10e\t"
            "%.10e\t%.10e\t%.10e\t"
            "%.10e\t%.10e\t%.10e\n",
            particle->index, i,
            particle->success[i],
            particle->time[i],
            particle->position[i].x,
            particle->position[i].y,
            particle->position[i].z,
            particle->velocity[i].x,
            particle->velocity[i].y,
            particle->velocity[i].z,
            particle->force[i].x,
            particle->force[i].y,
            particle->force[i].z);
  }
}

void create_xyz_file(const int tindex, const Particles *particle) {
  int index = particle->index;

  char filename[80];
  sprintf(filename, "./results/trajectories_%06d/traj_%09d.xyz", index, tindex);
  FILE *ffile = fopen(filename, "w");
  fprintf(ffile, "2\n");
  fprintf(ffile, "Trajectory of particles\n");
  for (unsigned int i = 0; i < particle->no; ++i) {
    fprintf(ffile, "C %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
            particle->position[i].x,
            particle->position[i].y,
            particle->position[i].z,
            particle->velocity[i].x,
            particle->velocity[i].y,
            particle->velocity[i].z,
            particle->force[i].x,
            particle->force[i].y,
            particle->force[i].z);
  }
  fclose(ffile);
}

void simulate_particle(const struct config *conf,
                       const double radius,
                       const double *force_list,
                       Particles *particle) {

  double dt = conf->dt;
  hid_t file_id;
  const int parindex = particle->index;
  const int saveflag = conf->saveflag;
  const int xyzflag = conf->xyzflag;
  const double mass = conf->mass;

  Vector3D ang_vel = {0., 0., 0.};
  Vector3D acc = {0., 0., 0.};

  char ffname[80];
  // FILE *ffile = NULL;
  int tindex = 0;

  printf("Simulating particle %d\n", parindex);


  if (xyzflag) {
    char dirname[80];
    sprintf(dirname, "./results/trajectories_%06d", parindex);

    // Create result directory to keep results
    struct stat st = {0};
    if (stat(dirname, &st) == -1) {
      mkdir(dirname, 0700);
    }
  }

  if (saveflag) {
    sprintf(ffname, "./results/trajectory_%06d.h5", parindex);
    file_id = initialize_file(ffname);
  }

  int sucess = 0;
  double time = 0;

  while (sucess || time < conf->max_time) {
    // if (tindex % 100) printf("%d\t", tindex);
    if (saveflag)
      save_timestep_to_hdf5(file_id, tindex, particle);
    if (xyzflag)
      create_xyz_file(tindex, particle);

    for (unsigned int i = 0; i < particle->no; ++i) {
      if (!particle->success[i]) {
        find_accelration(mass, radius, particle->position + i, particle->force + i,
                        &acc);
        // ω(t + Δt/2) =  ω(t) + Δt/2 ɑ(t)
        update_angvel(dt/2.0, particle->ang_velocity + i, &acc, &ang_vel);
        // q(t+Δt) = q(t) + Δt/2 q(t)ω(t + Δt/2)
        update_orientation(radius, dt, particle->orientation + i, &ang_vel, particle->position + i);
      }
    }

    // Find the force at half-step
    calculate_force(particle, conf, force_list);
    for (unsigned int i = 0; i < particle->no; ++i) {
      if (!particle->success[i]) {
        find_accelration(mass, radius, particle->position + i, particle->force + i,
                        &acc);
        update_angvel(dt/2.0, &ang_vel, &acc, particle->ang_velocity + i);
        cross_product(particle->position + i, particle->ang_velocity + i,
                      particle->velocity + i);
        particle->velocity_sq[i] = norm_sq(particle->velocity + i);
        particle->time[i] += dt;
      }
    }
    tindex += 1;

    time = particle->time[0];
    sucess = particle->success[0];
    for (unsigned int i = 1; i < particle->no; ++i) {
      sucess *= particle->success[i];
      if (time < particle->time[i])
        time = particle->time[i];
    }
  }

  if (saveflag)
    close_file(file_id);
}
