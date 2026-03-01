#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "utils.h"


int open_file(FILE **f, const char *file_name){
    
    *f = fopen(file_name, "r");
    if (*f == NULL){
        printf("Error opening the file %s\n", file_name);
        return 1;
    }    
    //printf("File openned!\n");
    return 0;
}

int open_file_w(FILE **f, const char *file_name){
    
    *f = fopen(file_name, "a"); // append mode, no overwriting
    if (*f == NULL){
        printf(" > Error opening the file %s\n", file_name);
        return 1;
    }    
    return 0;
}
