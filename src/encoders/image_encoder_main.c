#include <time.h>

#include "MNIST_for_C/mnist.h"
#include "encoders/image_encoders.h"


/* main.c */
int main(int argc, char *argv[]) {
    
    // randomize execution
    srand(time(NULL));

    int i, j;
    int *train_labels, *test_labels;
    int n_train_images, n_test_images;
    FILE *train_set_data_file, *test_set_data_file; // file to store spike trains
    FILE *train_set_labels_file, *test_set_labels_file; // file to store labels
    image_dataset_t train_spike_images, test_spike_images;
    int bins, image_size;


    // load dataset // TOdo generalize?
    load_mnist(); // This should be generalized????


    /* Initialize parameters // TODO: generalize? */
    // initialize number of samples in each dataset // TODO: this should be read from a meta.txt file
    n_train_images = 60000; 
    n_test_images = 10000;
    // compute bins
    //int pT = 255;
    //int pDT = 1;
    bins = 255; // number of time steps for the bin
    image_size = 28*28;
    // set parameters for datasets
    train_spike_images.bins = bins;
    train_spike_images.image_size = image_size; // this should be loaded in a function from a meta file
    train_spike_images.n_images = n_train_images;

    test_spike_images.bins = bins;
    test_spike_images.image_size = image_size;
    test_spike_images.n_images = n_test_images;

    /* Allocate memory */
    // allocate memory for labels
    train_labels = (int *)malloc(n_train_images * sizeof(int));
    test_labels = (int *)malloc(n_test_images * sizeof(int));

    // allocate memory for samples
    train_spike_images.images = malloc(train_spike_images.n_images * sizeof(spike_image_t)); // alloc memory for train images
    test_spike_images.images = malloc(test_spike_images.n_images * sizeof(spike_image_t)); // alloc memory for test images
    for(i=0; i<train_spike_images.n_images; i++){
        // alloc memory for image pixels
        train_spike_images.images[i].image = malloc(train_spike_images.image_size * sizeof(int *)); // each image has 728 pixels (this must be generalized)

        for(j=0; j<train_spike_images.image_size; j++){ // alloc bins integers for each pixel spike train --> this must be corrected, as bins is used for 0,1 spike trains, but I use spike times, so less positions are needed
            train_spike_images.images[i].image[j] = malloc(bins * sizeof(int));
        }
    }
    for(i=0; i<test_spike_images.n_images; i++){
        // alloc memory for image pixels
        test_spike_images.images[i].image = malloc(test_spike_images.image_size * sizeof(int *)); // each image has 728 pixels (this must be generalized)

        for(j=0; j<test_spike_images.image_size; j++){ // alloc bins integers for each pixel spike train --> this must be corrected, as bins is used for 0,1 spike trains, but I use spike times, so less positions are needed
            test_spike_images.images[i].image[j] = malloc(bins * sizeof(int));
        }
    }


    // read and convert input data
    read_convert_and_store_input_data(train_set_data_file, "./data/MNIST/train_data.txt", &train_spike_images, train_image);
    read_convert_and_store_input_data(test_set_data_file, "./data/MNIST/test_data.txt", &test_spike_images, test_image);

    // read labels
    read_labels(train_set_labels_file, "./data/MNIST/train_labels.txt", n_train_images, train_label);
    read_labels(test_set_labels_file, "./data/MNIST/test_labels.txt", n_test_images, test_label);

    return 0;
}