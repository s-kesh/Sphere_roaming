#ifndef H5FILE_H_
#define H5FILE_H_

#include "simulation.h"

#include <H5Ipublic.h>
#include <H5public.h>
#include <hdf5.h>

herr_t write_1d_array(hid_t groupid, const char *name, hid_t type, void *data, hsize_t size);
hid_t initialize_file(const char *filename);
void save_timestep_to_hdf5(const hid_t file_id, const unsigned int step, const Particle *particles);
void close_file(const hid_t file_id);
#endif
