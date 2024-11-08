#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include "simulation.h"

#include <sys/types.h>
#include <sys/stat.h>

void printUsage(const char *progName) {
  fprintf(stderr, "Usage: %s -s <start_number> -a <step_number> -e <end_number> -v <velocity> -n <number>\n", progName);
  exit(EXIT_FAILURE);
}

int main(int argc, char *argv[]) {
  int opt;
  int helium_number = -1;
  double radius = -1.0;
  double velocity = -1.0;
  int number = -1;

  // int start = -1;
  // int step = -1;
  // int end = -1;

  while ((opt = getopt(argc, argv, "v:n:")) != -1) {
    switch (opt) {
     // case 's':
     //   start = atof(optarg);
     //   break;
     // case 'a':
     //   step = atof(optarg);
     //   break;
     // case 'e':
     //   end = atof(optarg);
     //   break;
      case 'v':
        velocity = atof(optarg);
        break;
      case 'n':
        number = atoi(optarg);
        break;
      default:
        printUsage(argv[0]);
    }
  }

  // Check if both arguments were provided
  // if (start == -1 || step == -1 || end == -1 || velocity == -1.0 || number == -1) {
  if (velocity == -1.0 || number == -1) {
    printUsage(argv[0]);
  }

  char fitfilename[80] = "./data/fit_force.txt";
  FILE *fitfile = fopen(fitfilename, "r");

  if (fitfile == NULL) {
      perror("Error opening file");
      return 1;
  }

  // Read the potential file and store the values into arrays
  double* rF = (double *)malloc(1000000 * sizeof(double)); // radius array
  double* Flist = (double *)malloc(1000000 * sizeof(double)); // Singlet_Single array
  double yst; // Singlet_Triplet array
  double ytt; // Triplet_Triplet array

  int count = 0;
  while (fscanf(fitfile, "%lf\t%lf\t%lf\t%lf\n", rF + count, Flist + count, &yst, &ytt) == 4) {
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
  fprintf(fffile, "Radius\tForce\n");
  for (int i = 0; i < count; ++i) {
    fprintf(fffile, "%lf\t%lf\n", rF[i], Flist[i]);
  }

  // Array to hold twin particles
  int noP = number; // No of particle pairs
  struct Particle_Pair *particles = (struct Particle_Pair *)malloc(noP * sizeof(struct Particle_Pair));
  char ffname[80];
  FILE *ffile = NULL;

  // for (helium_number = start; helium_number < end; helium_number += step) {
    helium_number = 10000;
    radius = 2.22 * pow((double)helium_number, 1.0/3.0);

    // Initialize random particleSimulation
    initialize_particle_pair(radius, velocity, particles, noP);

    // Write particle properties in a file
    sprintf(ffname, "./results/particlefile_%d.txt", helium_number);
    ffile = fopen(ffname, "w");
    fprintf (ffile,"Index\tSuccess\tTime\tDistance\t"
             "X_1\tY_1\tZ_1\t"
             "X_2\tY_2\tZ_2\t"
             "WX_1\tWY_1\tWZ_1\t"
             "WX_2\tWY_2\tWZ_2\n");
    for (int i = 0; i < noP; ++i) {
      fprintf(ffile, "%d\t%d\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\n",
              particles[i].index, particles[i].success,
              particles[i].time, particles[i].distance,
              particles[i].p1_pos.x,
              particles[i].p1_pos.y,
              particles[i].p1_pos.z,
              particles[i].p2_pos.x,
              particles[i].p2_pos.y,
              particles[i].p2_pos.z,
              particles[i].p1_angvel.x,
              particles[i].p1_angvel.y,
              particles[i].p1_angvel.z,
              particles[i].p2_angvel.x,
              particles[i].p2_angvel.y,
              particles[i].p2_angvel.z);
    }
    fclose(ffile);

    for (int i = 0; i < noP; ++i)
      simulate_particle(radius, count-1, rF, Flist, particles + i);

    // Write particle properties in a file
    sprintf(ffname, "./results/particlefile_after_%d.txt", helium_number);
    ffile = fopen(ffname, "w");
    fprintf (ffile,"Index\tSuccess\tTime\tDistance\t"
             "X_1\tY_1\tZ_1\t"
             "X_2\tY_2\tZ_2\t"
             "WX_1\tWY_1\tWZ_1\t"
             "WX_2\tWY_2\tWZ_2\n");
    for (int i = 0; i < noP; ++i) {
      fprintf(ffile, "%d\t%d\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\t"
              "%.10lf\t%.10lf\t%.10lf\n",
              particles[i].index, particles[i].success,
              particles[i].time, particles[i].distance,
              particles[i].p1_pos.x,
              particles[i].p1_pos.y,
              particles[i].p1_pos.z,
              particles[i].p2_pos.x,
              particles[i].p2_pos.y,
              particles[i].p2_pos.z,
              particles[i].p1_angvel.x,
              particles[i].p1_angvel.y,
              particles[i].p1_angvel.z,
              particles[i].p2_angvel.x,
              particles[i].p2_angvel.y,
              particles[i].p2_angvel.z);
    }
    fclose(ffile);
    free(particles);
  // }

  return 0;
}
