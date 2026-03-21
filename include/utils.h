#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/// @brief Functon to open a file and check if a error occurs
/// @param f FILE struct to store the opened file
/// @param file_name path of the file to open
int open_file(FILE **f, const char *file_name);

/// @brief Function to open a file in write mode and check if a error occurs
/// @param f FILE struct to store the opened file
/// @param file_name path of the file to open
int open_file_w(FILE **f, const char *file_name);

/// @brief Get maximum value in the array
/// @param arr Array
/// @param n Number of elements in the array
/// @return Maximum value in the array
int get_max_value(int *arr, size_t n);


#ifdef __cplusplus
}
#endif

#endif