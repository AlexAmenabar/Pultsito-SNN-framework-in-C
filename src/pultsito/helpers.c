#include "helpers.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


void print_matrix( int *matrix, int n_row, int n_col){
    for(int i = 0; i<n_row; i++){
        for(int j = 0; j<n_col; j++){
            printf("%d ", matrix[i*n_col + j]);
        }
        printf("\n");
    }
}

void print_array(int *array, int length){
    int i;
    for(i=0; i<length; i++){
        printf("%d ", array[i]);
    }
    printf("\n");
}

void print_matrix_uint8(__uint8_t *matrix, int n_row, int n_col){
    for(int i = 0; i<n_row; i++){
        for(int j = 0; j<n_col; j++){
            printf("%d ", matrix[i*n_col + j]);
        }
        printf("\n");
    }
}

void print_array_uint8(__uint8_t *array, int length){
    int i;
    for(i=0; i<length; i++){
        printf("%d ", array[i]);
    }
    printf("\n");
}


void print_matrix_double(double *matrix, int n_row, int n_col){
    for(int i = 0; i<n_row; i++){
        for(int j = 0; j<n_col; j++){
            printf("%f ", matrix[i*n_col + j]);
        }
        printf("\n");
    }
}

void print_array_double(double *array, int length){
    int i;
    for(i=0; i<length; i++){
        printf("%f ", array[i]);
    }
    printf("\n");
}