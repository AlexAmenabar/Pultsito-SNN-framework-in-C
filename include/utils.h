#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

/// @brief Functon to open a file and check if a error occurs
/// @param f FILE struct to store the opened file
/// @param file_name path of the file to open
int open_file(FILE **f, const char *file_name);

/// @brief Function to open a file in write mode and check if a error occurs
/// @param f FILE struct to store the opened file
/// @param file_name path of the file to open
int open_file_w(FILE **f, const char *file_name);

#endif