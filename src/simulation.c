#include "simulation.h"
#include "vector.h"
#include <math.h>
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
      else if (strcmp(key, "seed") == 0)
        conf->seed = atoi(value);
      else if (strcmp(key, "mass") == 0)
        conf->mass = atof(value);
      else if (strcmp(key, "velocity") == 0)
        conf->velocity = atof(value);
      else if (strcmp(key, "velratio") == 0)
        conf->velratio = atof(value);
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

double distancesq_between_pair(const struct Particle_Pair *particle) {
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

    return xx*xx + yy*yy + zz*zz;
}

double generate_velocity(const struct config *conf) {
  double pi = 4 * atan(1);

  // Generate 4 uniformly distributed numbers
  double u1 = (double)rand() / RAND_MAX;
  double u2 = (double)rand() / RAND_MAX;
  double u3 = (double)rand() / RAND_MAX;
  double u4 = (double)rand() / RAND_MAX;

  // Box-Muller algorithim
  double z1 = sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2);
  double z2 = sqrt(-2.0 * log(u1)) * sin(2.0 * pi * u2);
  double z3 = sqrt(-2.0 * log(u3)) * cos(2.0 * pi * u4);
  // double z4 = sqrt(-2.0 * log(u3)) * sin(2.0 * pi * u4);

  return sqrt(z1*z1 + z2*z2 + z3*z3) * conf->velocity;
}

int initialize_particle_pair(const double radius,
                             const struct config *conf,
                             const double *force_list,
                             struct Particle_Pair *particle) {
  double dist = 0.0;
  double inv_grid_size = 1.0 / conf->force_grid;
  double force_start = conf->force_grid_start;

  struct Particle_Pair *par = NULL;
  Quat q1 = {0, 0, 0, 0};
  Vector3D pos1 = {0, 0, 0};
  Quat q2 = {0, 0, 0, 0};
  Vector3D pos2 = {0, 0, 0};
  Vector3D res = {0, 0, 0};

  // const double velocity_sq = conf->velocity * conf->velocity * 1E-10; // Convert m/s to Å/fs
  const double icd_dist2 = conf->icd_dist*conf->icd_dist;
  const double max_dist_sq = conf->force_grid_start + (conf->force_grid_length - 2) * conf->force_grid;

  // Setup random number generator
  if (conf->seed)
    srand(conf->seed);
  else
    srand(time(NULL));

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
        dist = distance_sq(&pos1, &pos2);
      } while (dist < icd_dist2);


      // Set pair position
      par->index = i;
      par->success = 0;
      par->time = 0;
      par->distance_sq = dist;

      par->p1_qua = q1;
      par->p2_qua = q2;
      par->p1_pos = pos1;
      par->p2_pos = pos2;

      // Generate two random speed
      double v1 = generate_velocity(conf) * 1E-5; // Convert m/s to Å/fs
      double v2 = generate_velocity(conf) * 1E-5; // Convert m/s to Å/fs
      par->p1_velocity_sq = v1*v1;
      par->p2_velocity_sq = v2*v2;
      double ang_vel1 = v1 / radius; // Unit fs¯¹
      double ang_vel2 = v2 / radius; // Unit fs¯¹

      // Generate random angular velocity
      // double ratio = conf->velratio;
      // par->p1_velocity_sq = (2.0 / (ratio*ratio + 1)) * velocity_sq;
      // par->p2_velocity_sq = (2.0 * ratio * ratio / (ratio*ratio + 1)) * velocity_sq;
      // double ang_vel1 = sqrt(par->p1_velocity_sq) / radius; // Unit fs¯¹
      // double ang_vel2 = sqrt(par->p2_velocity_sq) / radius; // Unit fs¯¹

      // Need to find the perpenducular vector to position
      create_perpend_vector(&(par->p1_pos), &res);
      scalar_multiply(ang_vel1/norm(&res), &res, &(par->p1_angvel));
      create_perpend_vector(&(par->p2_pos), &res);
      scalar_multiply(ang_vel2/norm(&res), &res, &(par->p2_angvel));

      // Calculate velocity
      cross_product(&(par->p1_pos), &(par->p1_angvel), &(par->p1_vel));
      cross_product(&(par->p2_pos), &(par->p2_angvel), &(par->p2_vel));

      // Find force initially
      if (particle->distance_sq > max_dist_sq) {
        particle->force = 0;
        particle->p1_force.x = 0.0;
        particle->p1_force.y = 0.0;
        particle->p1_force.z = 0.0;
        particle->p2_force.x = 0.0;
        particle->p2_force.y = 0.0;
        particle->p2_force.z = 0.0;
      }
      else {
        particle->force = force_list[(int)((particle->distance_sq - force_start) * inv_grid_size)];
        particle->p1_force.x = particle->force * (particle->p1_pos.x - particle->p2_pos.x) / sqrt(particle->distance_sq);
        particle->p1_force.y = particle->force * (particle->p1_pos.y - particle->p2_pos.y) / sqrt(particle->distance_sq);
        particle->p1_force.z = particle->force * (particle->p1_pos.z - particle->p2_pos.z) / sqrt(particle->distance_sq);
        scalar_multiply(-1, &(particle->p1_force), &(particle->p2_force));
      }
  }
  return 0;
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

void save_particle_pair(FILE *fpointer, const struct Particle_Pair *particle) {
  fprintf(fpointer,
          "%d\t%d\t%.10lf\t%.10lf\t"
          "%.10e\t%.10e\t%.10e\t"
          "%.10e\t%.10e\t%.10e\t%.10e\t"
          "%.10e\t%.10e\t%.10e\t%.10e\t"
          "%.10e\t%.10e\t%.10e\t"
          "%.10e\t%.10e\t%.10e\t"
          "%.10e\t%.10e\t%.10e\t"
          "%.10e\t%.10e\t%.10e\t"
          "%.10e\t%.10e\t%.10e\t"
          "%.10e\t%.10e\t%.10e\t"
          "%.10e\t%.10e\t%.10e\t"
          "%.10e\t%.10e\t%.10e\n",
          particle->index, particle->success, particle->time,
          particle->distance_sq,
          particle->p1_velocity_sq,
          particle->p2_velocity_sq,
          particle->force,
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
  fprintf(ffile, "C %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
          particle->p1_pos.x,
          particle->p1_pos.y,
          particle->p1_pos.z,
          particle->p1_vel.x,
          particle->p1_vel.y,
          particle->p1_vel.z,
          particle->p1_force.x,
          particle->p1_force.y,
          particle->p1_force.z);
  fprintf(ffile, "C %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e %.10e\n",
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
                       const double radius,
                       const double *force_list,
                       struct Particle_Pair *particle) {

  double dt = conf->dt;
  const int parindex = particle->index;
  const int saveflag = conf->saveflag;
  const int xyzflag = conf->xyzflag;
  const double icd_dist2 = conf->icd_dist * conf->icd_dist;
  const double max_time = conf->max_time;
  const double max_dist_sq = conf->force_grid_start + (conf->force_grid_length - 2) * conf->force_grid;
  const double inv_grid_size = 1.0/conf->force_grid;
  const double grid_start = conf->force_grid_start;
  const double mass = conf->mass;

  Vector3D ang_vel1 = {0., 0., 0.};
  Vector3D acc1 = {0., 0., 0.};
  Vector3D ang_vel2 = {0., 0., 0.};
  Vector3D acc2 = {0., 0., 0.};

  char ffname[80];
  FILE *ffile = NULL;
  int tindex = 0;


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
    sprintf(ffname, "./results/trajectory_%06d.%02d", parindex, (int)(conf->velratio*10));
    ffile = fopen(ffname, "w");
    fprintf(ffile,
            "Index\tSuccess\tTime\tDistance_Sq\t"
            "Vel_Sq_1\tVel_Sq_2\tForce\t"
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

  // Calculate initial torque and acceleration
  if (fabs(particle->force) > 1E-10) {
    dt = 1.0;
    find_accelration(mass, radius, &(particle->p1_pos), &(particle->p1_force), &acc1);
    find_accelration(mass, radius, &(particle->p2_pos), &(particle->p2_force), &acc2);
  }

  while (particle->time < max_time && !particle->success) {
    // Save trajectory in text file
    if (saveflag)
      save_particle_pair(ffile, particle);
    if (xyzflag)
      create_xyz_file(tindex, particle);

    if (fabs(particle->force) > 1E-10) {
      dt = 1.0;
      // Step 1: Calculate angular velocity at timestamp Δt/2
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

    // Step 2: Update Orientation
    update_orientation(radius, dt, &(particle->p1_qua),
                       &ang_vel1, &(particle->p1_pos));
    update_orientation(radius, dt, &(particle->p2_qua),
                       &ang_vel2, &(particle->p2_pos));

    // Check if they met
    particle->distance_sq = distancesq_between_pair(particle);
    if (particle->distance_sq < icd_dist2)
      particle->success = 1;

    // Check force again
    if (particle->distance_sq > max_dist_sq)
      particle->force = 0;
    else
      particle->force = force_list[(int)((particle->distance_sq - grid_start) * inv_grid_size)];
    if (fabs(particle->force) > 1E-10) {
      // Step 4: Find the torque at new position
      particle->p1_force.x = particle->force * (particle->p1_pos.x - particle->p2_pos.x) / sqrt(particle->distance_sq);
      particle->p1_force.y = particle->force * (particle->p1_pos.y - particle->p2_pos.y) / sqrt(particle->distance_sq);
      particle->p1_force.z = particle->force * (particle->p1_pos.z - particle->p2_pos.z) / sqrt(particle->distance_sq);
      scalar_multiply(-1, &(particle->p1_force), &(particle->p2_force));

      find_accelration(mass, radius, &(particle->p1_pos), &(particle->p1_force), &acc1);
      find_accelration(mass, radius, &(particle->p2_pos), &(particle->p2_force), &acc2);

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

    // Calculate new velocity
    cross_product(&(particle->p1_pos), &(particle->p1_angvel), &(particle->p1_vel));
    particle->p1_velocity_sq = norm_sq(&(particle->p1_vel));
    cross_product(&(particle->p2_pos), &(particle->p2_angvel), &(particle->p2_vel));
    particle->p2_velocity_sq = norm_sq(&(particle->p2_vel));

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
