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

    // Parse key-value pairs
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

void calculate_force(Particles *pars,
                     const struct config *conf,
                     const float *Force_list) {

  const float inv_grid_size = 1.0 / conf->force_grid;
  const float force_start = conf->force_grid_start;
  const float icd_dist2 = conf->icd_dist*conf->icd_dist;
  // const float max_dist_sq = conf->force_grid_start + (conf->force_grid_length - 2) * conf->force_grid;

  unsigned int num = pars->no;
  float *dist_sq = (float *)malloc((num * (num - 1) / 2) * sizeof(float));
 
  // Initialize forces on each particle to zero
  for (unsigned int j = 0; j < pars->no; ++j) {
    Particle *particle_j = pars->particle + j;
    particle_j->force.x = 0.0;
    particle_j->force.y = 0.0;
    particle_j->force.z = 0.0;
    particle_j->force_mag = 0.0;
  }

  for (unsigned int j = 0; j < pars->no; ++j) {
    // Take a pointer to particle j
    Particle *particle_j = pars->particle + j;
    if (particle_j->success) continue;

    // If pair of j,k is close to each other
    // Set success flag to 1
    for (unsigned int k = j + 1; k < pars->no; ++k) {
      Particle *particle_k = pars->particle + k;
      if (particle_k->success) continue;
      unsigned int index = j*num + k - (j + 1)*(j + 2)/2;
      dist_sq[index] = distance_sq(&(particle_j->position), &(particle_k->position));
      if (dist_sq[index] < icd_dist2) {
        particle_j->success = 1;
        particle_k->success = 1;
      }
    }
  }

  // Calculate forces
  for (unsigned int j = 0; j < pars->no - 1; ++j) {
    Particle *particle_j = pars->particle + j;
    if (particle_j->success) continue;
    Vector3D f;
    f.x = 0;
    f.y = 0;
    f.z = 0;

    for (unsigned int k = j + 1; k < pars->no; ++k) {
      Particle *particle_k = pars->particle + k;
      if (particle_k->success) continue;
      unsigned int index = j*num + k - (j + 1)*(j + 2)/2;
      float d = dist_sq[index];
      int d_index = (int)((d - force_start) * inv_grid_size);
      if (d_index < conf->force_grid_length) {
        float scale = Force_list[d_index] / sqrt(d);
        for (unsigned index = 0; index < 3; ++index) {
          f.v[index] = scale * (particle_j->position.v[index] - particle_k->position.v[index]);
          particle_j->force.v[index] += f.v[index];
          particle_k->force.v[index] -= f.v[index];
        }
      }
    }
    particle_j->force_mag = norm(&(particle_j->force));
  }

  free(dist_sq);
} 

float generate_velocity(const struct config *conf) {
  float pi = 4 * atan(1);

  // Generate 4 uniformly distributed numbers
  float u1 = (float)rand() / RAND_MAX;
  float u2 = (float)rand() / RAND_MAX;
  float u3 = (float)rand() / RAND_MAX;
  float u4 = (float)rand() / RAND_MAX;

  // Box-Muller algorithim
  // Generate 3 normally distributed numbers
  float z1 = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2);
  float z2 = sqrt(-2.0 * log(u1)) * sin(2.0 * pi * u2);
  float z3 = sqrt(-2.0 * log(u3)) * cos(2.0 * pi * u4);
  // float z4 = sqrt(-2.0 * log(u3)) * sin(2.0 * pi * u4);

  return sqrt(z1*z1 + z2*z2 + z3*z3) * conf->velocity;
}

int initialize_particles(const float radius,
                             const struct config *conf,
                             const float *Force_list,
                             Particles *particles) {

  // Setup random number generator
  if (conf->seed)
    srand(conf->seed);
  else
    srand(time(NULL));

  Particles *pars = NULL;
  for (int i = 0; i < conf->number; ++i) {
      pars = particles + i;

      pars->no = conf->no;
      pars->index = i;
      // Initializing indivisual particles
      pars->particle = (Particle *)malloc(conf->no * sizeof(Particle));
      for (int j = 0; j < conf->no; ++j) {
        Particle *par = pars->particle + j;
        par->index = j;
        par->success = 0;

        // Generate Position
        quat_random(&(par->orientation));

        quat_to_pos(&(par->orientation), &(par->position));
        scalar_multiply(radius, &(par->position), &(par->position));

        // Generate velocity
        float vel = generate_velocity(conf) * 1E-5; // Convert m/s to Å/fs
        float ang_vel = vel / radius; // Unit fs¯¹
        par->velocity_sq = vel*vel;
        Vector3D res;
        create_perpend_vector(&(par->position), &res);
        scalar_multiply(ang_vel/norm(&res), &res, &(par->ang_velocity));
        cross_product(&(par->position), &(par->ang_velocity),
                      &(par->velocity));
      }

      // Calculate forces
      calculate_force(pars, conf, Force_list);
  }
  return 0;
}

void free_particles(Particles *particles) {
  free(particles->particle);
}

void find_accelration(const float mass,
                      const float radius,
                      const Vector3D *pos,
                      const Vector3D *force,
                      Vector3D *acc) {

  Vector3D torque;
  cross_product(pos, force, &torque);
  scalar_multiply(1.0 / mass / radius / radius, &torque, acc);
}

void update_angvel(const float timestep,
                   const Vector3D *oldangvel,
                   const Vector3D *acc,
                   Vector3D *newangvel) {
  Vector3D temp;
  scalar_multiply(timestep, acc, &temp);
  add_vectors(oldangvel, &temp, newangvel);
}

/*
// SPIRAL
void update_orientation(const float radius,
                        const float timestep,
                        Quat *orient,
                        const Vector3D *angvel,
                        Vector3D *pos) {
  Quat tempquat = {0, 0, 0, 0};
  Quat res = {0, 0, 0, 0};
  Vector3D temp = {0, 0, 0};
  float mag = norm(angvel);
  float theta = timestep/2.0 * mag;

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
void update_orientation(const float radius,
                        const float timestep,
                        Quat *orient,
                        const Vector3D *angvel,
                        Vector3D *pos) {
  Quat Omega;
  Omega.a = 0;
  Omega.x = angvel->x;
  Omega.y = angvel->y;
  Omega.z = angvel->z;
  Quat temp;
  quat_product(&Omega, orient, &temp);
  quat_scalar_multiply(&temp, timestep/2.0, &temp);
  quat_add(orient, &temp, orient);
  quat_scalar_multiply(orient,
                       1.0 / quat_norm(orient),
                       orient);
  quat_to_pos(orient, pos);
  scalar_multiply(radius, pos, pos);
}

void save_particles(FILE *fpointer, const Particles *particles) {
  for (unsigned int i = 0; i < particles->no; ++i) {
    Particle *par = particles->particle + i;
    fprintf(fpointer,
            "%d\t%d\t%d\t%.10e\t"
            "%.10e\t%.10e\t%.10e\t"
            "%.10e\t%.10e\t%.10e\t"
            "%.10e\t%.10e\t%.10e\n",
            particles->index, i,
            par->success,
            par->time,
            par->position.x,
            par->position.y,
            par->position.z,
            par->velocity.x,
            par->velocity.y,
            par->velocity.z,
            par->force.x,
            par->force.y,
            par->force.z);
  }
}

void simulate_particles(const struct config *conf,
                       const float radius,
                       const float *force_list,
                       Particles *particles) {

  float dt = conf->dt;
  hid_t file_id;
  const int parindex = particles->index;
  const int saveflag = conf->saveflag;
  const float mass = conf->mass;

  Vector3D acc;
  Vector3D *ang_velocities = (Vector3D *)malloc(particles->no * sizeof(Vector3D));
  for (unsigned int i = 0; i < particles->no; ++i) {
    ang_velocities[i].x = 0;
    ang_velocities[i].y = 0;
    ang_velocities[i].z = 0;
  }

  char ffname[100];
  int tindex = 0;

  if (saveflag) {
    sprintf(ffname, "./results/trajectory_%06d.h5", parindex);
    file_id = initialize_file(ffname);
  }

  int sucess = 0;
  float time = 0;

  while (!sucess && time < conf->max_time) {
    /*
    if (tindex % 1000 == 0) printf("%d\t%f\t%f\t%d\t%d\n",
      tindex, time, conf->max_time, sucess, time < conf->max_time);
    */
    
    if (saveflag)
      save_timestep_to_hdf5(file_id, tindex, particles);

    for (unsigned int i = 0; i < particles->no; ++i) {
      Particle *particle = particles->particle + i;

      // First half-step
      if (!particle->success) {
        find_accelration(mass, radius, &(particle->position), &(particle->force),
                        &acc);
        // ω(t + Δt/2) =  ω(t) + Δt/2 ɑ(t)
        update_angvel(dt/2.0, &(particle->ang_velocity), &acc, ang_velocities + i);
        // q(t+Δt) = q(t) + Δt/2 q(t)ω(t + Δt/2)
        update_orientation(radius, dt, &(particle->orientation), ang_velocities + i, &(particle->position));
      }
    }

    // Find the force at half-step
    calculate_force(particles, conf, force_list);

    // Second half-step
    for (unsigned int i = 0; i < particles->no; ++i) {
      Particle *particle = particles->particle + i;
      if (!particle->success) {
        find_accelration(mass, radius, &(particle->position), &(particle->force),
                        &acc);
        // ω(t + Δt) = w(t + Δt/2) + Δt/2 ɑ(t)
        update_angvel(dt/2.0, ang_velocities + i, &acc, &(particle->ang_velocity));
        // Set velocity; v = ω x r
        cross_product(&(particle->position), &(particle->ang_velocity),
                      &(particle->velocity));
        particle->velocity_sq= norm_sq(&(particle->velocity));
        particle->time += dt;
      }
    }

    tindex += 1;
    sucess = particles->particle->success;
    time = particles->particle->time;
    for (unsigned int i = 1; i < particles->no; ++i) {
      Particle *particle = particles->particle + i;
      sucess *=  particle->success;
      if (time < particle->time)
        time = particle->time;
    }
  }

  if (saveflag)
    close_file(file_id);

  free(ang_velocities);
}