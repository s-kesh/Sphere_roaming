#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include "simulation.h"

#include <sys/types.h>
#include <sys/stat.h>

void print_configuration(const struct config *conf) {
    printf("Configuration File Loaded\n");
    printf("no = %d\n", conf->no);
    printf("number = %d\n", conf->number);
    printf("saveflag = %d\n", conf->saveflag);
    printf("xyzflag = %d\n", conf->xyzflag);
    printf("start = %d\n", conf->start);
    printf("step = %d\n", conf->step);
    printf("stop = %d\n", conf->stop);
    printf("helium_number = %d\n", conf->helium_number);
    printf("seed = %ld\n", conf->seed);
    printf("mass = %lf\n", conf->mass);
    printf("velocity = %lf\n", conf->velocity);
    printf("velratio = %lf\n", conf->velratio);
    printf("dt = %lf\n", conf->dt);
    printf("max_time = %lf\n", conf->max_time);
    printf("icd_dist = %lf\n", conf->icd_dist);
    printf("force_grid = %lf\n", conf->force_grid);
    printf("force_grid_start = %lf\n", conf->force_grid_start);
    printf("force_grid_length = %ld\n", conf->force_grid_length);
    printf("force_file = %s\n", conf->force_file);
}


void simulate_for_number(const int helium_number, const struct config *conf,
                         const double *Flist) {
  // Calculate radius of droplet
  double radius = 2.22 * pow((double)helium_number, 1.0/3.0);

  // Array to hold particles
  Particles *particles = (Particles *)malloc(conf->number * sizeof(Particles));

  // Initialize random particleSimulation
  initialize_particle_pair(conf->no, radius, conf, Flist, particles);

  // Write particle properties in a file
  char ffname[80];
  FILE *ffile;
  sprintf(ffname, "./results/particlefile_%d.%02d", helium_number, (int)(conf->velratio*10));
  ffile = fopen(ffname, "w");
  fprintf(ffile,
          "Index\tNumber\tSuccess\tTime\t"
          "X\tY\tZ\t"
          "VX\tVY\tVZ\t"
          "FX\tFY\tFZ\n");
  for (int i = 0; i < conf->number; ++i) {
    save_particle_pair(ffile, particles + i);
  }
  fclose(ffile);

  // Simulate using velocity verlet

  #pragma omp parallel for
  for (int i = 0; i < conf->number; ++i)
    simulate_particle(conf, radius, Flist, particles + i);

  // Write particle properties in a file
  sprintf(ffname, "./results/particlefile_after_%d.%02d", helium_number, (int)(conf->velratio*10));
  ffile = fopen(ffname, "w");
  fprintf(ffile,
          "Index\tNumber\tSuccess\tTime\t"
          "X\tY\tZ\t"
          "VX\tVY\tVZ\t"
          "FX\tFY\tFZ\n");
  for (int i = 0; i < conf->number; ++i)
    save_particle_pair(ffile, particles + i);

  for (int i = 0; i < conf->number; ++i)
    free_particles(particles + i);

  fclose(ffile);
  free(particles);

}

int main(int argc, char *argv[]) {

  if (argc < 2) {
    fprintf(stderr, "Usage: %s <config_file>\n", argv[0]);
    return 1;
  }

  struct config con;
  int success = read_config(argv[1], &con);

  // Print configuration
  if (success) {
    if (con.print)
      print_configuration(&con);
  } else {
    perror("Error reading configuration file");
    return 1;
  }

  FILE *fitfile = fopen(con.force_file, "r");
  if (fitfile == NULL) {
      perror("Error opening file");
      return 1;
  }

  // Read the potential file and store the values into arrays
  double* Flist = (double *)malloc(con.force_grid_length * sizeof(double)); // Singlet_Single array
  double rF;
  double yst; // Singlet_Triplet array
  double ytt; // Triplet_Triplet array

  int count = 0;
  while (fscanf(fitfile, "%lf\t%lf\t%lf\t%lf\n", &rF, Flist + count, &yst, &ytt) == 4 && count < con.force_grid_length) {
    count++;
  }
  fclose(fitfile);

  // Create result directory to keep results
  struct stat st = {0};
  if (stat("./results", &st) == -1) {
    mkdir("./results", 0700);
  }

  // Save forcelist
  char fffname[80] = "./results/forcelist.txt";
  FILE *fffile = fopen(fffname, "w");
  fprintf(fffile, "Distance_Squared\tForce\n");
  for (int i = 0; i < count; ++i) {
    fprintf(fffile, "%lf\t%.10e\n", con.force_grid_start + con.force_grid*i, Flist[i]);
  }

  int helium_number = con.helium_number;
  if (con.start)
    for (helium_number = con.start; helium_number < con.stop;
         helium_number += con.step)
      simulate_for_number(helium_number, &con, Flist);
  else
    simulate_for_number(helium_number, &con, Flist);

  free(Flist);
  return 0;
}
