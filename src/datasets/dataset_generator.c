#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


// Function to generate random spike trains
int random_spike_train_generator(int *spike_train, int time_steps, int prob){
    
    int p, t, n_spks = 0;

    // generate spikes 
    for(t=0; t<time_steps; t++){
        
        p = rand() % 100;

        // generate spike on time t if p > prob
        if(p > prob){
            
            spike_train[n_spks] = t;
            n_spks ++;
        }
    }

    // return number of spikes
    return n_spks;
}


int main(int argc, char *argv[]) {

    int i, j;
    int n_samples, n_input_spike_trains, time_steps, prob;
    int *spk_train, n_spikes;
    FILE *f;

    
    // randomize
    srand(time(NULL)); 

    // read the number of samples in the dataset
    n_samples = strtoul(argv[1], NULL, 10);

    // read the number of spike trains to generate (the number of features of each sample)
    n_input_spike_trains = strtoul(argv[2], NULL, 10);

    // read the time steps for each spike train
    time_steps = strtoul(argv[3], NULL, 10);

    // variable to store the spike train
    spk_train = (int *)malloc(time_steps * sizeof(int));

    // open file
    f = fopen(argv[4], "w"); //f = fopen("input_spikes.txt", "w");
    if(f == NULL){
        printf(" > Error opening the file!\n");
        exit(1);
    }


    for(size_t i = 0; i<n_samples; i++){

        for(size_t j = 0; j<n_input_spike_trains; j++){
            
            int prob = rand() % 100; // compute a random probability for the spike train
            
            // generate random spike train
            n_spikes = random_spike_train_generator(spk_train, time_steps, prob);

            // store spike train
            fprintf(f, "%d ", n_spikes);

            for(size_t s = 0; s<n_spikes; s++){
                fprintf(f, "%d ", spk_train[s]);
            }
            fprintf(f, "\n");
        }
    }

    free(spk_train);

    // close file
    fclose(f);
}