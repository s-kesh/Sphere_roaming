/*
Simulation entry point file.
This file reads the configuration file and the force file and then simulates the particles using the velocity verlet algorithm.
The results are saved in the results directory.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include "simulation.h"

#include <sys/types.h>
#include <sys/stat.h>

/*
Prints the configuration read from the configuration file.
@param conf: Configuration structure containing the configuration read from the file.
*/
void print_configuration(const struct config *conf) {
    printf("Configuration File Loaded\n");
    printf("hv = %f\n", conf->hv);
    printf("intensity = %f\n", conf->intensity);
    printf("cross_section = %f\n", conf->cross_section);
    printf("pulse_width = %f\n", conf->pulse_width);
    printf("number = %d\n", conf->number);
    printf("saveflag = %d\n", conf->saveflag);
    printf("start = %d\n", conf->start);
    printf("step = %d\n", conf->step);
    printf("stop = %d\n", conf->stop);
    printf("helium_number = %d\n", conf->helium_number);
    printf("seed = %ld\n", conf->seed);
    printf("mass = %f\n", conf->mass);
    printf("velocity = %f\n", conf->velocity);
    printf("dt = %f\n", conf->dt);
    printf("max_time = %f\n", conf->max_time);
    printf("icd_dist = %f\n", conf->icd_dist);
    printf("force_grid = %f\n", conf->force_grid);
    printf("force_grid_start = %f\n", conf->force_grid_start);
    printf("force_grid_length = %ld\n", conf->force_grid_length);
    printf("force_file = %s\n", conf->force_file);
}

/*
Simulates the particles for a given helium number.
@param helium_number: Number of helium atoms in the droplet.
@param conf: Configuration structure containing the configuration read from the file.
@param Flist: Array containing the force values for different distances.
*/
void simulate_for_number(const int helium_number, const struct config *conf,
                         const float *Flist) {

  // Array to hold particles
  // Number of particles = number_of_simulations
  Particles *particles = (Particles *)malloc(conf->number * sizeof(Particles));

  // Initialize random particleSimulation
  initialize_particles(helium_number, conf, Flist, particles);

  // Write particle properties in a file
  char ffname[80];
  FILE *ffile;
  sprintf(ffname, "./results/particlefile_%d", helium_number);
  ffile = fopen(ffname, "w");
  fprintf(ffile,
          "Index\tNumber\tSuccess\tTime\t"
          "X\tY\tZ\t"
          "VX\tVY\tVZ\t"
          "FX\tFY\tFZ\n");
  for (int i = 0; i < conf->number; ++i)
    save_particles(ffile, particles + i);
  fclose(ffile);

  // Simulate using velocity verlet

  #pragma omp parallel for
  for (int i = 0; i < conf->number; ++i)
    simulate_particles(conf, helium_number, Flist, particles + i);

  // Write particle properties in a file
  sprintf(ffname, "./results/particlefile_after_%d", helium_number);
  ffile = fopen(ffname, "w");
  fprintf(ffile,
          "Index\tNumber\tSuccess\tTime\t"
          "X\tY\tZ\t"
          "VX\tVY\tVZ\t"
          "FX\tFY\tFZ\n");
  for (int i = 0; i < conf->number; ++i)
    save_particles(ffile, particles + i);

  for (int i = 0; i < conf->number; ++i)
    free_particles(particles + i);

  fclose(ffile);
  free(particles);

}

/*
Main function to read the configuration file and the force file and simulate the particles.
*/
int main(int argc, char *argv[]) {

  if (argc < 2) {
    fprintf(stderr, "Usage: %s <config_file>\n", argv[0]);
    return 1;
  }

  // Read configuration file
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

  // Open force file
  FILE *fitfile = fopen(con.force_file, "r");
  if (fitfile == NULL) {
      perror("Error opening file");
      return 1;
  }

  // Read the force file and store the values into arrays
  float* Flist = (float *)malloc(con.force_grid_length * sizeof(float)); // Singlet_Single array
  float rF;
  float yst; // Singlet_Triplet array
  float ytt; // Triplet_Triplet array

  int count = 0;
  while (fscanf(fitfile, "%f\t%f\t%f\t%f\n", &rF, Flist + count, &yst, &ytt) == 4 && count < con.force_grid_length) {
    count++;
  }
  fclose(fitfile);

  // Create result directory to keep results
  struct stat st = {0};
  if (stat("./results", &st) == -1) {
    mkdir("./results", 0700);
  }

  // Save forcelist to check if force file is read correctly
  char fffname[80] = "./results/forcelist.txt";
  FILE *fffile = fopen(fffname, "w");
  fprintf(fffile, "Distance_Squared\tForce\n");
  for (int i = 0; i < count; ++i) {
    fprintf(fffile, "%f\t%.10e\n", con.force_grid_start + con.force_grid*i, Flist[i]);
  }

  // Simulate for the given helium number
  // If start, stop and step are given, simulate for all helium numbers in the range
  // Else simulate for the given helium number
  int helium_number = con.helium_number;
  if (con.start)
    for (helium_number = con.start; helium_number < con.stop;
         helium_number += con.step)
      simulate_for_number(helium_number, &con, Flist);
  else
    simulate_for_number(helium_number, &con, Flist);

  // Free memory
  free(Flist);
  return 0;
}
