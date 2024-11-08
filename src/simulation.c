#include "simulation.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define ICD_DIST 6.0 // Å
#define M_HE 4.0 // amu
#define MAX_TIME 100000000 // fs

double distance_between_pair(struct Particle_Pair *particle) {
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

double distance(struct Vector3D *v1, struct Vector3D *v2) {

    double xx = (v1->x - v2->x);
    double yy = (v1->y - v2->y);
    double zz = (v1->z - v2->z);

    return sqrt(xx*xx + yy*yy + zz*zz);
}

double distance_values(double x1, double y1, double z1,
                       double x2, double y2, double z2) {

    double xx = (x1 - x2);
    double yy = (y1 - y2);
    double zz = (z1 - z2);

    return sqrt(xx*xx + yy*yy + zz*zz);
}

void cross_product(struct Vector3D *v1,
                   struct Vector3D *v2,
                   struct Vector3D *res) {
  res->x = v1->y * v2->z - v1->z * v2->y;
  res->y = v1->z * v2->x - v1->x * v2->z;
  res->z = v1->x * v2->y - v1->y * v2->x;
}

double norm(struct Vector3D *v1) {
  return sqrt(v1->x*v1->x + v1->y*v1->y + v1->z*v1->z);
}

void normalize(struct Vector3D *vec, struct Vector3D *res) {
    double nor = norm(vec);
    res->x = vec->x / nor;
    res->y = vec->y / nor;
    res->z = vec->z / nor;
}

void add_vectors(struct Vector3D *v1,
                 struct Vector3D *v2,
                 struct Vector3D *res) {
  res->x = v1->x + v2->x;
  res->y = v1->y + v2->y;
  res->z = v1->z + v2->z;
}

void scalar_multiply(double scalar,
                     struct Vector3D *vec,
                     struct Vector3D *res) {
  res->x = vec->x * scalar;
  res->y = vec->y * scalar;
  res->z = vec->z * scalar;
}

double dot_product(struct Vector3D *v1,
                   struct Vector3D *v2) {
  return v1->x*v2->x + v1->y*v2->y + v1->z*v2->z;
}

void rotate_vector(struct Vector3D *v,
                   struct Vector3D *axis,
                   double angle,
                   struct Vector3D* result) {
    double cosTheta = cos(angle);
    double sinTheta = sqrt(1 - cosTheta*cosTheta);
    double dot = dot_product(axis, v); // Compute dot product of axis and vector

    // Rodrigues' rotation formula
    result->x = cosTheta * v->x + sinTheta * (axis->y * v->z - axis->z * v->y) + (1 - cosTheta) * dot * axis->x;
    result->y = cosTheta * v->y + sinTheta * (axis->z * v->x - axis->x * v->z) + (1 - cosTheta) * dot * axis->y;
    result->z = cosTheta * v->z + sinTheta * (axis->x * v->y - axis->y * v->x) + (1 - cosTheta) * dot * axis->z;
}

int initialize_particle_pair(double radius,
                             double velocity,
                             struct Particle_Pair *particle,
                             int no_of_particles) {
  const double pi = 4 * atan(1);
  double distance = 0.0;
  double theta = 0.0;
  double phi = 0.0;
  double x1, y1, z1;
  double x2, y2, z2;
  struct Particle_Pair *par = NULL;

  // Setup random number generator
  srand(100401430);
  // srand(time(NULL));

  for (int i = 0; i < no_of_particles; ++i) {
      par = particle + i;
      // Generate random position
      theta = pi * ((double)rand() / RAND_MAX);
      phi = 2 * pi * ((double)rand() / RAND_MAX);
      x1 = radius * cos(phi) * sin(theta);
      y1 = radius * sin(phi) * sin(theta);
      z1 = radius * cos(theta);

      do {
        theta = pi * ((double)rand() / RAND_MAX);
        phi = 2 * pi * ((double)rand() / RAND_MAX);
        x2 = radius * cos(phi) * sin(theta);
        y2 = radius * sin(phi) * sin(theta);
        z2 = radius * cos(theta);
        distance = distance_values(x1, y1, z1, x2, y2, z2);
      } while (distance < ICD_DIST);


      // Set pair position
      par->index = i;
      par->success = 0;
      par->time = 0;
      par->distance = distance;

      par->p1_pos.x = x1;
      par->p1_pos.y = y1;
      par->p1_pos.z = z1;
      par->p2_pos.x = x2;
      par->p2_pos.y = y2;
      par->p2_pos.z = z2;

      // Generate random angular velocity
      velocity = 1E-5 * velocity; // Convert m/s to Å/fs
      par->p1_velocity = velocity;
      par->p2_velocity = velocity;
      double ang_vel = velocity / radius; // Unit: fs¯¹

      // Need to find the perpenducular vector to position
      // Create a vector 'a' which is not parallel to position vector
      struct Vector3D a = {0, 0, 0};
      struct Vector3D res = {0, 0, 0};
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
  }

  return 0;
}

double find_force(double distance, int list_size,
                  double *rlist, double *force_list) {
  int indx = 0;
  if (distance < rlist[0]) indx = 0;
  else if (distance > rlist[list_size - 1]) indx = list_size - 2;
  else indx = (int)((distance - rlist[0]) / 0.001);
  return -1 * force_list[indx] * 9.6485332E-03;
}

void simulate_particle(double radius,
                       int list_size,
                       double *rlist, double *force_list,
                       struct Particle_Pair *particle) {

  int saveflag = 0;
  double DT = 1;

  double forcemag = 0.0;
  double norm_angvel1 = 0.0;
  double angle1 = 0.0;
  double norm_angvel2 = 0.0;
  double angle2 = 0.0;
  struct Vector3D force1 = {0., 0., 0.};
  struct Vector3D force2 = {0., 0., 0.};
  struct Vector3D torque1 = {0., 0., 0.};
  struct Vector3D torque2 = {0., 0., 0.};
  struct Vector3D ang_vel1 = {0., 0., 0.};
  struct Vector3D acc1 = {0., 0., 0.};
  struct Vector3D ang_vel2 = {0., 0., 0.};
  struct Vector3D acc2 = {0., 0., 0.};
  struct Vector3D rotvector1 = {0., 0., 0.};
  struct Vector3D rotvector2 = {0., 0., 0.};
  struct Vector3D res = {0., 0., 0.};

  char ffname[80];
  FILE *ffile = NULL;

  if (saveflag) {
    sprintf(ffname, "./results/trajectory_%06d.txt", particle->index);
    ffile = fopen(ffname, "w");
    fprintf(ffile, "Index\tSuccess\tTime\tDistance\t"
            "Vel_1\tVel_2\t"
             "X_1\tY_1\tZ_1\t"
             "X_2\tY_2\tZ_2\t"
             "WX_1\tWY_1\tWZ_1\t"
             "WX_2\tWY_2\tWZ_2\n");
  }

  particle->time = 0;
  while (particle->time < MAX_TIME && !particle->success) {
    // Save trajectory in text file
    if (saveflag) {
      fprintf(ffile, "%d\t%d\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\n",
              particle->index, particle->success,
              particle->time, particle->distance,
              particle->p1_velocity, particle->p2_velocity,
              particle->p1_pos.x,
              particle->p1_pos.y,
              particle->p1_pos.z,
              particle->p2_pos.x,
              particle->p2_pos.y,
              particle->p2_pos.z,
              particle->p1_angvel.x,
              particle->p1_angvel.y,
              particle->p1_angvel.z,
              particle->p2_angvel.x,
              particle->p2_angvel.y,
              particle->p2_angvel.z);
    }

    forcemag = find_force(particle->distance, list_size, rlist, force_list);
    if (fabs(forcemag - 1E-10) > 0) {
      // Step 1: Find the torque on particles
      force1.x = forcemag * (particle->p1_pos.x - particle->p2_pos.x) / particle->distance;
      force1.y = forcemag * (particle->p1_pos.y - particle->p2_pos.y) / particle->distance;
      force1.z = forcemag * (particle->p1_pos.z - particle->p2_pos.z) / particle->distance;
      scalar_multiply(-1, &force1, &force2);
      cross_product (&(particle->p1_pos), &force1, &(torque1));
      cross_product (&(particle->p2_pos), &force2, &(torque2));

      // Step 2: Calculate angular velocity at timestamp Δt/2
      scalar_multiply(DT / 2.0 / M_HE / radius / radius, &torque1, &acc1);
      add_vectors(&(particle->p1_angvel), &acc1, &ang_vel1);
      scalar_multiply(DT / 2.0 / M_HE / radius / radius, &torque2, &acc2);
      add_vectors(&(particle->p2_angvel), &acc2, &ang_vel2);
    }

    // Find rotation angle and rotation axis
    norm_angvel1 = norm(&ang_vel1);
    angle1 = norm_angvel1 * DT;
    norm_angvel2 = norm(&ang_vel2);
    angle2 = norm_angvel2 * DT;
    scalar_multiply(1.0/norm_angvel1, &ang_vel1, &rotvector1);
    scalar_multiply(1.0/norm_angvel2, &ang_vel2, &rotvector2);

    // Step 3: Update the position (Rotate it)
    rotate_vector(&(particle->p1_pos), &rotvector1, angle1, &res);
    particle->p1_pos.x = res.x;
    particle->p1_pos.y = res.y;
    particle->p1_pos.z = res.z;
    rotate_vector(&(particle->p2_pos), &rotvector2, angle2, &res);
    particle->p2_pos.x = res.x;
    particle->p2_pos.y = res.y;
    particle->p2_pos.z = res.z;

    forcemag = find_force(particle->distance, list_size, rlist, force_list);
    if (fabs(forcemag - 1E-10) > 0) {
      // Step 4: Find the torque at new position
      force1.x = forcemag * (particle->p1_pos.x - particle->p2_pos.x) / particle->distance;
      force1.y = forcemag * (particle->p1_pos.y - particle->p2_pos.y) / particle->distance;
      force1.z = forcemag * (particle->p1_pos.z - particle->p2_pos.z) / particle->distance;
      scalar_multiply(-1, &force1, &force2);
      cross_product (&(particle->p1_pos), &force1, &(torque1));
      cross_product (&(particle->p2_pos), &force2, &(torque2));

      // Step 5: Calculate angular velocity at Δt
      scalar_multiply(DT / 2.0 / M_HE / radius / radius, &torque1, &acc1);
      add_vectors(&ang_vel1, &acc1, &(particle->p1_angvel));
      scalar_multiply(DT / 2.0 / M_HE / radius / radius, &torque2, &acc2);
      add_vectors(&ang_vel2, &acc2, &(particle->p2_angvel));
    }

    // Calculate new velocity magnitute
    cross_product(&(particle->p1_pos), &(particle->p1_angvel), &res);
    particle->p1_velocity = norm(&res);
    cross_product(&(particle->p2_pos), &(particle->p2_angvel), &res);
    particle->p2_velocity = norm(&res);

    // Find distance between particles
    particle->distance = distance_between_pair(particle);
    particle->time += DT;
    // Check if they met
    if (particle->distance < ICD_DIST) particle->success = 1;
  }

  if (saveflag) {
    // Write particle properties in a file
    fprintf(ffile, "%d\t%d\t%.10lf\t%.10lf\t"
            "%.10lf\t%.10lf\t"
            "%.10lf\t%.10lf\t%.10lf\t"
            "%.10lf\t%.10lf\t%.10lf\t"
            "%.10lf\t%.10lf\t%.10lf\t"
            "%.10lf\t%.10lf\t%.10lf\n",
            particle->index, particle->success,
            particle->time, particle->distance,
            particle->p1_velocity, particle->p2_velocity,
            particle->p1_pos.x,
            particle->p1_pos.y,
            particle->p1_pos.z,
            particle->p2_pos.x,
            particle->p2_pos.y,
            particle->p2_pos.z,
            particle->p1_angvel.x,
            particle->p1_angvel.y,
            particle->p1_angvel.z,
            particle->p2_angvel.x,
            particle->p2_angvel.y,
            particle->p2_angvel.z);

    fclose(ffile);
  }

}
