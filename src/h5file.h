#ifndef H5FILE_H_
#define H5FILE_H_

#include "simulation.h"

#include <H5Ipublic.h>
#include <H5public.h>
#include <hdf5.h>

/*
Function to initialize the HDF5 file.
@param filename: Name of the file to be created.
@return: File ID of the created file.
*/
hid_t initialize_file(const char *filename);

/*
Function to open the HDF5 file.
@param filename: Name of the file to be opened.
@return: File ID of the opened file.
*/
hid_t open_file(const char *filename);

/*
Function to save the particle properties to the HDF5 file.
@param file_id: File ID of the HDF5 file.
@param step: Timestep of the simulation.
@param particles: Pointer to the particles structure.
@return: Error code.
*/
herr_t save_timestep_to_hdf5(const hid_t file_id, const unsigned int step, const Particles *particles);

/*
Function to read the particle properties from the HDF5 file.
@param file_id: File ID of the HDF5 file.
@param step: Timestep of the simulation.
@param particles: Pointer to the particles structure.
@return: Error code.
*/
herr_t read_timestep_from_hdf5(const hid_t file_id, const unsigned int step, Particles *particles);

/*
Function to close the HDF5 file.
@param file_id: File ID of the HDF5 file.
*/
void close_file(const hid_t file_id);
#endif
