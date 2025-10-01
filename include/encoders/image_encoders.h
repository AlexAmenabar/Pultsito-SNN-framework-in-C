#ifndef IMAGE_ENCODERS_H
#define IMAGE_ENDOERS_H



/** 
This file contains the structures and functions to convert different types of input into spike trains and store them.
*/


// TODO: this should be generalized for different data types ??? 

/// @brief Image of spikes. Each image has a list of spikes for each pixel
typedef struct {
    int **image; // first dimension is a list of pixels, and the second dimension the spike train for each pixel
} spike_image_t;

/// @brief Image dataset converted into spikes
typedef struct {
    spike_image_t *images; // array of images (spike trains)
    int *labels; // array of labels
    int n_classes; // number of classes in the dataset
    int image_size; // image size
    int n_images; // number of images in the dataset
    int bins; // number of bins for each pixel
} image_dataset_t;

/// @brief Read the labels of a dataset
/// @param f File to write labels on
/// @param f_path file path to write labels on
/// @param n Number of samples in the dataset
/// @param labels pointer to read the labels from
void read_labels(FILE *f, char *f_path, int n, int *labels);


/// @brief Read input data, convert it to spikes and store in a file
/// @param f File to store the dataset in spikes format
/// @param f_path File path of the output file
/// @param spk_ds Dataset in spikes
/// @param ds Original dataset
void read_convert_and_store_input_data(FILE *f, char *f_path, image_dataset_t *spk_ds, double **ds);


/// @brief Convert dataset of images into spike images
/// @param spike_images Struct to store the image converted
/// @param images Images to be converted into spikes
/// @param n_images Number of images in the dataset to be converted into spikes
/// @param image_size Number of pixels into the image
/// @param bins number of spikes (max) into each pixel spike train
void convert_images_to_spikes_by_poisson_distribution(image_dataset_t *spike_images, double **images, int n_images, int image_size, int bins);


/// @brief Convert a image into spike trains
/// @param spike_image Struct to store the image converted
/// @param image Image to be converted into spikes
/// @param image_size Number of pixels into the image
/// @param bins number of spikes (max) into each pixel spike train
void convert_image_to_spikes_by_poisson_distribution(spike_image_t *spike_image, double *image, int image_size, int bins);

#endif