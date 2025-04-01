#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include "h5file.h"

#define PI 3.14159265358979323846
#define N_LAT 20
#define N_LON 20


void setup_gnuplot(FILE *gnuplot) {
    fprintf(gnuplot, "set terminal pngcairo enhanced size 800,800\n");
    fprintf(gnuplot, "unset border\n");
    fprintf(gnuplot, "unset key\n");
    fprintf(gnuplot, "set xzeroaxis\n");
    fprintf(gnuplot, "set yzeroaxis\n");
    fprintf(gnuplot, "set zzeroaxis\n");
    fprintf(gnuplot, "set xyplane at 0\n");
    fprintf(gnuplot, "unset tics\n");
    fprintf(gnuplot, "set parametric\n");
    fprintf(gnuplot, "set xrange [-1.2:1.2]\n");
    fprintf(gnuplot, "set yrange [-1.2:1.2]\n");
    fprintf(gnuplot, "set zrange [-1.2:1.2]\n");
    fprintf(gnuplot, "set view equal xyz\n");
    fprintf(gnuplot, "set size ratio -1\n");
    fprintf(gnuplot, "set autoscale fix\n");
}

void generate_sphere_data(char **lat_data, size_t *lat_size, char **lon_data, size_t *lon_size) {
    // Latitude lines
    FILE *lat_stream = open_memstream(lat_data, lat_size);
    for (int k = 1; k < N_LAT; k++) {
        const double theta = (PI * k) / N_LAT;
        const double z = cos(theta);
        const double r = sin(theta);
        for (int j = 0; j <= N_LON; j++) {
            const double phi = (2 * PI * j) / N_LON;
            fprintf(lat_stream, "%f %f %f\n", r * cos(phi), r * sin(phi), z);
        }
        fprintf(lat_stream, "\n");
    }
    fprintf(lat_stream, "e\n");
    fclose(lat_stream);

    // Longitude lines
    FILE *lon_stream = open_memstream(lon_data, lon_size);
    for (int j = 0; j < N_LON; j++) {
        const double phi = (2 * PI * j) / N_LON;
        for (int k = 0; k <= N_LAT; k++) {
            const double theta = (PI * k) / N_LAT;
            fprintf(lon_stream, "%f %f %f\n",
                   sin(theta) * cos(phi),
                   sin(theta) * sin(phi),
                   cos(theta));
        }
        fprintf(lon_stream, "\n");
    }
    fprintf(lon_stream, "e\n");
    fclose(lon_stream);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <filename> <time_step>\n", argv[0]);
        return 1;
    }

    const char *filename = argv[1];
    const int plot_time_step = atoi(argv[2]);

    // Open HDF5 file
    hid_t file_id = open_file(filename);
    if (file_id < 0) {
        fprintf(stderr, "Error opening HDF5 file '%s'.\n", filename);
        return 1;
    }

    // Read file metadata
    clock_t start = clock();
    hsize_t num_time;
    H5Gget_num_objs(file_id, &num_time);
    printf("Time steps in file: %llu\n", num_time);

    // Read particle data
    Particles *particles = malloc(num_time * sizeof(Particles));
    for (hsize_t i = 0; i < num_time; ++i) {
        particles[i].index = i;
        particles[i].no = 0;
        particles[i].particle = NULL;
        if (read_timestep_from_hdf5(file_id, i, &particles[i]) < 0) {
            fprintf(stderr, "Error reading timestep %llu\n", i);
            for (hsize_t j = 0; j < i; j++) free(particles[j].particle);
            free(particles);
            H5Fclose(file_id);
            return 1;
        }
    }
    H5Fclose(file_id);
    printf("Reading time: %.2fs\n", (double)(clock() - start) / CLOCKS_PER_SEC);

    // Calculate normalization radius
    unsigned int num_particles = particles[0].no;
    const Vector3D first_pos = particles[0].particle[0].position;
    const double radius = sqrt(first_pos.x*first_pos.x +
                             first_pos.y*first_pos.y +
                             first_pos.z*first_pos.z);

    // Write normalized trajectories to file
    start = clock();
    FILE *traj_file = fopen("trajectories.bin", "wb");
    if (!traj_file) {
        fprintf(stderr, "Failed to create trajectories.dat\n");
        for (hsize_t i = 0; i < num_time; i++) free(particles[i].particle);
        free(particles);
        return 1;
    }

    #pragma omp parallel for
    for (unsigned int p = 0; p < num_particles; p++) {
        char traj_filename[256];
        snprintf(traj_filename, sizeof(traj_filename), "trajectory_%u.bin", p);

        FILE *traj_file = fopen(traj_filename, "wb");
        if (!traj_file) {
            fprintf(stderr, "Error creating file for particle %u\n", p);
            continue;
        }

        // Write normalized positions directly
        for (hsize_t t = 0; t < num_time; t++) {
            const Vector3D pos = particles[t].particle[p].position;
            const float normalized[3] = {
                pos.x / radius,
                pos.y / radius,
                pos.z / radius
            };
            fwrite(normalized, sizeof(float), 3, traj_file);
        }
        fclose(traj_file);
    }
    printf("Trajectory write time: %.2fs\n", (double)(clock() - start) / CLOCKS_PER_SEC);

    // Cleanup particle data
    for (hsize_t i = 0; i < num_time; i++) free(particles[i].particle);
    free(particles);

    // Generate sphere geometry
    start = clock();
    char *lat_data, *lon_data;
    size_t lat_size, lon_size;
    generate_sphere_data(&lat_data, &lat_size, &lon_data, &lon_size);

    // Plot configuration
    const char *colors[] = {"#FF0000", "#FF007F", "#FF00FF", "#7F00FF", "#0000FF"};
    const int num_colors = sizeof(colors)/sizeof(colors[0]);

    // Parallel plotting
    #pragma omp parallel for schedule(dynamic)
    for (hsize_t t = 0; t < num_time; t += plot_time_step) {
        FILE *gnuplot = popen("gnuplot -p", "w");
        if (!gnuplot) {
            #pragma omp critical
            fprintf(stderr, "Failed to open gnuplot for timestep %llu\n", t);
            continue;
        }

        setup_gnuplot(gnuplot);
        fprintf(gnuplot, "set output 'images/step_%010llu.png'\n", t);
        fprintf(gnuplot, "set title 'Time %llu fs'\n", t);

        // Build plot command
        char cmd[16384];
        char *ptr = cmd;
        ptr += sprintf(ptr, "splot '-' w l lc 'gray', '-' w l lc 'gray'");


        for (unsigned int p = 0; p < num_particles; p++) {
            const char *color = colors[p % num_colors];
            char filename[256];
            snprintf(filename, sizeof(filename), "trajectory_%u.bin", p);

            // Trajectory line
            if (t > 0) {
                ptr += sprintf(ptr,
                    ", '%s' binary format='%%float%%float%%float'"
                    " every ::0::%llu using 1:2:3 w l lc rgb '%s'",
                    filename, t-1, color);
            }

            // Current position
            ptr += sprintf(ptr,
                ", '%s' binary format='%%float%%float%%float'"
                " every ::%llu::%llu using 1:2:3 w p pt 7 lc rgb '%s'",
                filename, t, t, color);
        }

        fprintf(gnuplot, "%s\n", cmd);
        fprintf(gnuplot, "%s%s", lat_data, lon_data);
        fprintf(gnuplot, "unset output\n");
        fflush(gnuplot);
        pclose(gnuplot);
    }

    free(lat_data);
    free(lon_data);
    printf("Total plotting time: %.2fs\n", (double)(clock() - start) / CLOCKS_PER_SEC);

    return 0;
}
