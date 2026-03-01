#include <stdlib.h>
#include <stdio.h>

#include "config/config_loader.h"
#include "datasets/datasets.h"
#include "utils.h"

/* [PUBLIC] */
GPU_dataset_t* allocate_dataset_str(size_t n_samples, size_t n_features, size_t n_classes, size_t n_spikes){

    GPU_dataset_t *dataset = (GPU_dataset_t*)malloc(sizeof(GPU_dataset_t));

    // load general dataset information from configuration struct
    dataset->n_classes = n_classes;
    dataset->n_samples = n_samples;
    dataset->n_features = n_features;
    dataset->n_spikes = n_spikes;

    // allocate memory for arrays
    dataset->sample_offset        = n_samples > 0                   ? (size_t*)malloc(n_samples * sizeof(size_t))              : NULL;
    dataset->n_spikes_per_feature = n_samples > 0 && n_features > 0 ? (size_t*)malloc(n_samples * n_features * sizeof(size_t)) : NULL;
    dataset->feature_offset       = n_samples > 0 && n_features > 0 ? (size_t*)malloc(n_samples * n_features * sizeof(size_t)) : NULL;
    dataset->freq                 = n_samples > 0 && n_features > 0 ? (size_t*)malloc(n_samples * n_features * sizeof(size_t)) : NULL;
    dataset->first_spk            = n_samples > 0 && n_features > 0 ? (size_t*)malloc(n_samples * n_features * sizeof(size_t)) : NULL;
    dataset->spikes               = n_spikes > 0                    ? (size_t*)malloc(n_spikes * sizeof(size_t))               : NULL;

    // return allocated structure
    return dataset;
}

void deallocate_dataset_str(GPU_dataset_t *dataset){
    
    // deallocate internal arrays
    if(dataset->n_spikes_per_feature) free(dataset->n_spikes_per_feature);
    if(dataset->sample_offset)        free(dataset->sample_offset);
    if(dataset->feature_offset)       free(dataset->feature_offset);
    if(dataset->spikes)               free(dataset->spikes);
    if(dataset->freq)                 free(dataset->freq);
    if(dataset->first_spk)            free(dataset->first_spk);

    // deallocate struct
    if(dataset)                        free(dataset);
}

GPU_dataset_t* load_dataset_from_file_cpu(const char *file_name, const char *labels_file_name, size_t n_samples, simulation_configuration_t *conf){

    size_t i, j, l;
    FILE *f = NULL;

    open_file(&f, file_name);

    // allocate memory for dataset (we don't know the number of spikes in the dataset)
    GPU_dataset_t *dataset = allocate_dataset_str(conf->n_samples, conf->input_size, conf->n_classes, 0);

    // count number of spikes in the dataset
    size_t **tmp_spikes = (size_t **)malloc(dataset->n_samples * dataset->n_features * sizeof(size_t*));
    size_t offset = 0;
    size_t n_spikes;
    for(i = 0; i<dataset->n_samples; i++){
        
        // set sample offset
        dataset->sample_offset[i] = offset;

        // loop over features
        for(j = 0; j<dataset->n_features; j++){

            // set feature offset
            dataset->feature_offset[i * dataset->n_features + j] = offset;

            // scan number of spikes of the feature
            fscanf(f, "%zu", &(n_spikes));

            // update offset
            offset += n_spikes;
            dataset->n_spikes += n_spikes;

            // store number of spikes
            dataset->n_spikes_per_feature[i * dataset->n_features + j] = n_spikes;

            // store the spike times of the feature in tmp_spikes
            tmp_spikes[i * dataset->n_features + j] = (size_t*)malloc(n_spikes * sizeof(size_t));
            for(l = 0; l < n_spikes; l++){
                fscanf(f, "%zu", &(tmp_spikes[i * dataset->n_features + j][l]));
            }
        }
    }

    // copy spikes to the dataset struct, compute frequencies and store first spike time
    dataset->spikes = (size_t*)malloc(dataset->n_spikes * sizeof(size_t));

    size_t next_spike = 0;
    for(i = 0; i<dataset->n_samples; i++){

        // loop over features
        for(j = 0; j<dataset->n_features; j++){


            n_spikes = dataset->n_spikes_per_feature[i * dataset->n_features + j];

            for(l = 0; l < n_spikes; l++){

                dataset->spikes[next_spike] = tmp_spikes[i * dataset->n_features + j][l];
                next_spike ++;
            }

            // store first spike and frequency
            if(n_spikes > 0){
                dataset->freq[i * dataset->n_features + j] = conf->max_input_spikes / n_spikes; // spikes each freq time steps
                dataset->first_spk[i * dataset->n_features + j] = dataset->spikes[dataset->feature_offset[i * dataset->n_features + j]];
            }
            else{
                dataset->freq[i * dataset->n_features + j] = 0; // spikes each freq time steps
                dataset->first_spk[i * dataset->n_features + j] = 0;
            }

            // deallocate memory
            free(tmp_spikes[i * dataset->n_features + j]);
        }
    }

    // free temporaly allocated memory
    free(tmp_spikes);

    // return dataset
    return dataset;
}

double get_dataset_size(GPU_dataset_t *dataset){
    
    size_t nS = dataset->n_samples;
    size_t nF = dataset->n_features;
    size_t nSpks = dataset->n_spikes;

    return (

        // general and scalars
        (sizeof(GPU_dataset_t) +
        sizeof(int) + 
        sizeof(size_t) * 4 +

        // arrays
        sizeof(size_t) * (nS + nS * nF * 4 + nSpks)) / 8.0
    );
}

void print_dataset(GPU_dataset_t *dataset){

    size_t i, j, l, next = 0;

    printf(" > Printing dataset: \n");
    for(i = 0; i<dataset->n_samples; i++){

        printf(" > > Sample %zu (offset = %zu)\n", i, dataset->sample_offset[i]);

        // loop over features
        for(j = 0; j<dataset->n_features; j++){

            printf(" > >> Feature %zu (offset = %zu): [", j, dataset->feature_offset[i * dataset->n_features + j]);

            // print spikes in feature
            for(l = 0; l<dataset->n_spikes_per_feature[i * dataset->n_features + j]-1; l++){

                printf("%zu, ", dataset->spikes[next]);
                next++;
            }
            printf("%zu]\n", dataset->spikes[next]);
            next++;
        }
    }
}