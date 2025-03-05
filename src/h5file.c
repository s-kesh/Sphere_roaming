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

herr_t write_1d_array(hid_t groupid, const char *name, hid_t type, void *data, hsize_t size) {
    hsize_t dims[1] = {size};
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    hid_t dataset_id = H5Dcreate(groupid, name, type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    herr_t status = H5Dwrite(dataset_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    return status;
}

hid_t initialize_file(const char *filename) {
    hid_t file_id;
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if (file_id < 0) {
        fprintf(stderr, "Error: Could not create HDF5 file %s\n", filename);
        return 0;
    }
    return file_id;
}


void close_file(const hid_t file_id) {
    H5Fclose(file_id);
}

void save_timestep_to_hdf5(const hid_t file_id, const unsigned int step,
                          const Particle *particles) {
    // Create a group
    char group_name[32];
    snprintf(group_name, sizeof(group_name), "/Step_%u", step);
    hid_t group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Save particle
    unsigned int no = particles->no;
    write_1d_array(group_id, "no",
                   H5T_NATIVE_UINT, (void *)&(particles->no), 1);
    write_1d_array(group_id, "index",
                   H5T_NATIVE_UINT, (void *)&(particles->index), 1);
    write_1d_array(group_id, "sucess",
                   H5T_NATIVE_UINT, (void *)(particles->success), no);
    write_1d_array(group_id, "time",
                   H5T_NATIVE_DOUBLE, (void *)(particles->time), no);
    write_1d_array(group_id, "velocity_sq",
                   H5T_NATIVE_DOUBLE, (void *)(particles->velocity_sq), no);
    write_1d_array(group_id, "force_mag",
                   H5T_NATIVE_DOUBLE, (void *)(particles->force_mag), no);
    write_1d_array(group_id, "orientation",
                   H5T_NATIVE_DOUBLE, (void *)(particles->orientation), 4*no);
    write_1d_array(group_id, "position",
                   H5T_NATIVE_DOUBLE, (void *)(particles->position), 3*no);
    write_1d_array(group_id, "velocity",
                   H5T_NATIVE_DOUBLE, (void *)(particles->velocity), 3*no);
    write_1d_array(group_id, "ang_velocity",
                   H5T_NATIVE_DOUBLE, (void *)(particles->ang_velocity), 3*no);
    write_1d_array(group_id, "force",
                   H5T_NATIVE_DOUBLE, (void *)(particles->force), 3*no);

    // Close the group and file
    H5Gclose(group_id);
}
