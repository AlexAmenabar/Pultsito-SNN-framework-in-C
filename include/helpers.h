#ifndef HELPERS_H
#define HELPERS_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>



/**
Helper functions not related only to SNN use cases
*/


/*#define print_matrix(x) _Generic((x), \
    int: print_matrix_int, \
    uint8: print_matrix_uint8, \
    double: print_matrix_double \
)(x)*/


/// @brief Print matrix of int
/// @param matrix Matrix to be printed
/// @param n_row Number of rows of the matrix
/// @param n_col Number of columns of the matrix
void print_matrix_int(int *matrix, int n_row, int n_col);


/// @brief Print array of int
/// @param array Array to be printed
/// @param n Number of elements of the array
void print_array(int *array, int n);

/// @brief Print matrix of uint8
/// @param matrix Matrix to be printed
/// @param n_row Number of rows of the matrix
/// @param n_col Number of columns of the matrix
void print_matrix_uint8(__uint8_t *matrix, int n_row, int n_col);

/// @brief Print an array of uint8
/// @param array Array to be printed
/// @param n Number of elements in the array
void print_array_uint8(__uint8_t *array, int n);


/// @brief Print matrix (of real numbers)
/// @param matrix Matrix to be printed
/// @param n_row Number of rows of the matrix
/// @param n_col Number of columns of the matrix
void print_matrix_f(double *matrix, int n_row, int n_col);


/// @brief Print an array (of real numbers)
/// @param array Array to be printed
/// @param n Number of elements of the array
void print_array_f(double *array, int n);


/// @brief Function to generate random spike trains of different spike densities
/// @param spike_trains Array of arrays (of different lengths) to store the spike trains
/// @param n_spikes Number of spikes for each spike train
/// @param n_input_spike_trains Number of spike trains to generate
/// @param time_steps Maximum number of time steps for each spike train 
/// @param prob Probability to generate a spike in a time step
void random_input_spike_train_generator(int **spike_trains, int *n_spikes, int n_input_spike_trains, int time_steps, int prob);

#endif