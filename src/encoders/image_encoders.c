#include <stdlib.h>
#include <stdio.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "encoders/image_encoders.h"



void read_labels(FILE *f, char *f_path, int n, int *labels){
    
    int i;

    // open file
    //f = fopen("./data/MNIST/train/labels.txt", "w");
    f = fopen(f_path, "w");
    if(f == NULL){
        printf("Error opening the file. \n");
        exit(1);
    }

    // read labels from the file
    for(i = 0; i<n; i++)
        fprintf(f, "%d ", labels[i]);

    // close file
    fclose(f);
}

void read_convert_and_store_input_data(FILE *f, char *f_path, image_dataset_t *spk_ds, double **ds){
    
    int i, j, l;

    // read dataset data and convert to spikes
    convert_images_to_spikes_by_poisson_distribution(spk_ds, ds, spk_ds->n_images, spk_ds->image_size, spk_ds->bins);

    // open file to store data    
    f = fopen(f_path, "w");
    if(f == NULL){
        printf("Error opening the file. \n");
        exit(1);
    }

    // write data in file
    for(i = 0; i<spk_ds->n_images; i++){
        for(j = 0; j<spk_ds->image_size; j++){
            for(l=0; l<spk_ds->images[i].image[j][0]+1; l++){
                fprintf(f, "%d ", spk_ds->images[i].image[j][l]);
            }
            fprintf(f, "\n");
        }
        fprintf(f, "\n");
    }
}


void convert_images_to_spikes_by_poisson_distribution(image_dataset_t *image_dataset, double **images, int n_images, int image_size, int bins){
    int i;

    // convert all images into spikes
    for(i=0; i<n_images; i++){
        convert_image_to_spikes_by_poisson_distribution(&(image_dataset->images[i]), images[i], image_size, bins);
    }
}

void convert_image_to_spikes_by_poisson_distribution(spike_image_t *spike_image, double *image, int image_size, int bins){

    double random_number;
    int next_index, i, j;

    // loop over the image pixels
    for(i = 0; i<image_size; i++){

        next_index = 0;

        // generate bins for the pixel
        for(j = 0; j<bins; j++){
            
            random_number = (double)rand() / (double)RAND_MAX; // generate a random number

            // if the image pixel is bigger than the random number, generate a spike
            if(image[i] > random_number)
            {
                next_index ++;
                spike_image->image[i][next_index] = j; // add spike on time t
                
            }
        }
        spike_image->image[i][0] = next_index; // store amount of spikes in spike train
    }
}


