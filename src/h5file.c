#include "h5file.h"
#include <H5Dpublic.h>
#include <H5Fpublic.h>
#include <H5Gpublic.h>
#include <H5Ipublic.h>
#include <H5Ppublic.h>
#include <H5Spublic.h>
#include <H5Tpublic.h>
#include <H5public.h>
#include <H5version.h>
#include <stdio.h>
#include <stdlib.h>

// Global HDF5 datatypes (defined once)
static hid_t vector3d_type = -1;
static hid_t quat_type = -1;
static hid_t particle_type = -1;


// Function to initialize HDF5 datatypes
static void initialize_hdf5_types() {
    if (vector3d_type < 0) {
        hsize_t dim_vector3d = 3;
        vector3d_type = H5Tarray_create(H5T_NATIVE_FLOAT, 1, &dim_vector3d);
    }
    if (quat_type < 0) {
        hsize_t dim_quat = 4;
        quat_type = H5Tarray_create(H5T_NATIVE_FLOAT, 1, &dim_quat);
    }
    if (particle_type < 0) {
      particle_type = H5Tcreate(H5T_COMPOUND, sizeof(Particle));
 
      H5Tinsert(particle_type, "index", HOFFSET(Particle, index),
                H5T_NATIVE_UINT);
      H5Tinsert(particle_type, "success", HOFFSET(Particle, success),
                H5T_NATIVE_UINT);
      H5Tinsert(particle_type, "time", HOFFSET(Particle, time),
                H5T_NATIVE_FLOAT);
      H5Tinsert(particle_type, "force_mag", HOFFSET(Particle, force_mag),
                H5T_NATIVE_FLOAT);
      H5Tinsert(particle_type, "velocity_sq", HOFFSET(Particle, velocity_sq),
                H5T_NATIVE_FLOAT);
      H5Tinsert(particle_type, "orientation", HOFFSET(Particle, orientation.q), quat_type);
      H5Tinsert(particle_type, "position", HOFFSET(Particle, position.v), vector3d_type);
      H5Tinsert(particle_type, "velocity", HOFFSET(Particle, velocity.v), vector3d_type);
      H5Tinsert(particle_type, "ang_velocity", HOFFSET(Particle, ang_velocity.v), vector3d_type);
      H5Tinsert(particle_type, "force", HOFFSET(Particle, force.v), vector3d_type);
    }
}

// Function to free HDF5 datatypes
static void cleanup_hdf5_types() {
    if (vector3d_type >= 0) {
      H5Tclose(quat_type);
      vector3d_type = -1;
    }
    if (quat_type >= 0) {
      H5Tclose(quat_type);
      quat_type = -1;
    }
    if (particle_type >= 0) {
      H5Tclose(particle_type);
      particle_type = -1;
    }
}

hid_t initialize_file(const char *filename) {
  hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  // Initialize datatypes once
  initialize_hdf5_types();

  return file_id;
}

// Save timestep data
herr_t save_timestep_to_hdf5(const hid_t file_id, const unsigned int step,
                           const Particles *data) {
  // printf("Step %d\n", (int)step);

  // Create a group
  char group_name[32];
  snprintf(group_name, sizeof(group_name), "/Step_%u", step);
  hid_t group_id =
      H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // Save particles
  hsize_t dims[1];
  dims[0] = data->no;
  hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
  hid_t dataset_id =
      H5Dcreate(group_id, "particles", particle_type, dataspace_id, H5P_DEFAULT,
                H5P_DEFAULT, H5P_DEFAULT);
  herr_t status = H5Dwrite(dataset_id, particle_type, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, (void *)(data->particle));
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Gclose(group_id);
  return(status);
}

// Close the HDF5 file
void close_file(const hid_t file_id) {
    cleanup_hdf5_types(); // Free HDF5 datatypes
    H5Fclose(file_id);
}
